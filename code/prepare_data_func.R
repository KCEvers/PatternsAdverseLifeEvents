##### PREPARE ADVERSE LIFE EVENTS DATA FROM SHP AND HILDA #####

# Master function to prepare SHP/HILDA data
prepare_data = function(name_dataset, filepath_base, event_dict, rerun = FALSE){

  # Set-up file paths and create necessary directories
  filepath_dataset = file.path(filepath_base, "data", name_dataset)
  filepath_deriv = file.path(filepath_base, "derivatives", name_dataset)
  dir.create(filepath_deriv, recursive = TRUE, showWarnings = FALSE)

  filepath_figs = file.path(filepath_base, "figs")
  dir.create(filepath_figs, showWarnings = FALSE)

  filepath_figs_dataset = file.path(filepath_figs, name_dataset)
  dir.create(filepath_figs_dataset, showWarnings = FALSE)

  filepath_df_binary = file.path(filepath_deriv, "df_binary.RDS")

  # Prepare either SHP or HILDA data
  if (name_dataset == "SHP"){

    prepare_SHP(filepath_dataset, filepath_df_binary, event_dict, rerun)

  } else if (name_dataset == "HILDA"){

    prepare_HILDA(filepath_dataset, filepath_df_binary, event_dict, rerun)

  }

  # Read cleaned data
  df_binary = readRDS(file.path(filepath_deriv, "df_binary.RDS"))

  # Get dataframes for analysis
  out = get_analysis_df(name_dataset, df_binary, event_dict, filepath_deriv, rerun)
  out$df_binary = df_binary
  out$event_dict = event_dict

  # Add filepaths to output
  out$filepath_deriv = filepath_deriv
  out$filepath_figs = filepath_figs
  out$filepath_figs_dataset = filepath_figs_dataset
  return(out)
}

# Prepare SHP data
prepare_SHP = function(filepath_dataset, filepath_df_binary, event_dict, rerun){

  # Create dictionary of subjective variables (not used)
  subj_dict = matrix(c(
    c("PC02", "Satisfaction with health status", "health_satisfaction", "positive"),
    c("PC44", "Satisfaction with life in general", "life_satisfaction", "positive"),
    c("PQL04", "Satisfaction with personal relationships", "relationship_satisfaction", "positive"),
    c("PC17",	"Depression, blues, anxiety: Frequency", "depression_anxiety_freq", "negative"),
    c("PC18",	"Frequency of energy and optimism", "energy_optimism_freq", "positive"),
    c("PC47",	"Emotions: Joy: Frequency", "joy_freq", "positive"),
    c("PC48",	"Emotions: Anger: Frequency", "anger_freq", "negative"),
    c("PC49",	"Emotions: Sadness: Frequency", "sadness_freq", "negative"),
    c("PC50",	"Emotions: Worry: Frequency", "worry_freq", "negative")
  ),
  ncol = 4, byrow = TRUE ) %>%
    magrittr::set_colnames(c("var_name", "description", "recode_var", "valence")) %>% as.data.frame() %>%
    mutate_all(~ stringr::str_replace(.x, ": yes, no", ""))

  # Define different types of variables
  id_vars = c("IDPERS", "IDHOUS", "IDSPOU", "YEAR", "STATUS", "SEX", "AGE", "CIVSTA", "COHAST",
              "PD29", # Partner
              "WSTAT", subj_dict$var_name
  )
  le_all_var = c("PL01", "PL01A", "PL01B", "PL06", "PL11", "PL16", "PL21", "PL26", "PL35")
  le_suff_var = c("PL05A", "PL05B", "PL05", "PL10", "PL15", "PL20", "PL25", "PL30")
  le_main_var = c("PL01R", "PL06", "PL11", "PL16", "PL21", "PL26", "PL36") # Main variables with life events

  if (rerun | !file.exists(filepath_df_binary)){

    # SHP: Long formatted file of all waves
    df_p_raw = haven::read_sav(file.path(filepath_dataset, "/Data_SPSS/SHP-Data-Longfile-SPSS/SHPLONG_P_USER.sav"))
    df_h_raw = haven::read_sav(file.path(filepath_dataset, "/Data_SPSS/SHP-Data-Longfile-SPSS/SHPLONG_H_USER.sav"))

    # Combine person with household file
    df_raw = merge(df_p_raw, df_h_raw %>% dplyr::select(IDHOUS, YEAR,
                                                        HLDFFS), all.x = TRUE) %>%
      # Following instructions to only select STATUS 0 and 1 from https://www.suz.uzh.ch/dataforstat/baa/data.html
      filter(STATUS %in% c(0, 1)) %>% arrange(IDPERS, YEAR)

    # Save memory
    rm(df_h_raw)
    rm(df_p_raw)

    # Recode column with attributes
    recode_col <- function(x, col_name, df) {
      recode_vec <- attributes(df[[col_name]])$labels

      # Reverse the name-value pairs
      dplyr::recode(x, !!!tibble::deframe(tibble::enframe(recode_vec) %>% dplyr::select(value, name)))
    }

    # Clean and prepare dataframe with occurence per person per household per year per event
    df_binary = df_raw %>%
      # Remove year 1999, no life event data
      dplyr::filter(YEAR != 1999) %>%
      # Mutate age to numeric, haven-labelled otherwise
      mutate(AGE = as.numeric(as.character(factor(AGE))),
             IDHOUS = as.character(IDHOUS),
             IDPERS = as.character(IDPERS)) %>%
      dplyr::select(all_of(c("YEAR", "IDPERS", "IDHOUS", "SEX", "AGE", le_suff_var, union(le_all_var, le_main_var),
                             subj_dict$var_name))) %>%

      # Rename subjective variables
      dplyr::rename_with(~ tibble::deframe(subj_dict %>% dplyr::select(var_name, recode_var)),
                         matches(sprintf("^%s", subj_dict$var_name))
      ) %>%
      mutate(across(matches(subj_dict$recode_var),
                    ~ as.numeric(as.character(.x)))
      ) %>%

      # Recode life event variables to 0 (did not occur) and 1 (occurred)
      mutate(across(
        all_of(c("SEX", le_all_var)),
        ~ recode_col(as.numeric(as.character(.x)), cur_column(), df_raw)
      )
      ) %>%
      # Make occurrence a binary 0, 1, or NA
      mutate(across(
        all_of(le_all_var),
        ~ ifelse(.x == "yes", 1, ifelse(.x == "no", 0, NA))
      )
      ) %>%
      # PL01 (illness or accident) was asked separately in 2021/2022 -> combine and fill into PL01
      mutate(PL01 = ifelse(YEAR >= 2021, PL01A + PL01B > 0, PL01)) %>%
      # Suffering from illness or accident - use median if both illness and accident happened
      mutate(PL05 = ifelse(YEAR >= 2021, median(PL05A, PL05B, na.rm = TRUE), PL05)) %>%
      dplyr::select(-c(PL01A, PL01B, PL05A, PL05B)) %>%

      # PL01R = Type of illness or accident
      # Make sure that if PL01 = 1, but PL01R = NA, it is coded as PL01R = "other"
      mutate(PL01R = ifelse(PL01 == 1 & is.na(PL01R),
                            tibble::enframe(attributes(df_raw$PL01R)$labels) %>%
                              filter(name == "other reason") %>% pull(value),
                            ifelse(PL01 == 0 & (is.na(PL01R) | PL01R == 0), 0,
                                   PL01R))) %>%


      # PL36 = Type of other life event
      mutate(PL36 = ifelse(
        # Make sure that if PL35 = 1, but PL36 = NA, it is coded as PL36 = "other"
        PL35 == 1 & is.na(PL36),
        tibble::enframe(attributes(df_raw$PL36)$labels) %>%
          filter(name == "other") %>% pull(value),
        # If PL35 = 0 but there is an event in PL36, code as the event in PL36
        ifelse(PL35 == 0 & (is.na(PL36) | PL36 == 0), 0,
               PL36))) %>%

      ## PL01R
      mutate(PL01R_char = as.character(PL01R)) %>%          # Treat NA as a valid category
      mutate(row_id = row_number()) %>%                     # Add a row identifier to preserve order
      tidyr::pivot_wider(names_from = PL01R_char,                  # Create indicator columns
                         values_from = PL01R_char,
                         values_fn = length,
                         names_prefix = "PL01R_") %>%
      select(-row_id) %>%  # Drop helper column
      mutate(across(starts_with("PL01R_"),
                    # If PL01R did occur and the current column isn't the event that occurred, set to zero
                    ~ ifelse(!is.na(PL01R) & cur_column() != paste0("PL01R_", PL01R), 0, .))) %>%
      # Remove original
      dplyr::select(-all_of(c("PL01", "PL01R", "PL01R_NA", "PL01R_0"))) %>%

      ## PL36
      mutate(PL36_char = as.character(PL36)) %>%          # Treat NA as a valid category
      mutate(row_id = row_number()) %>%                     # Add a row identifier to preserve order
      tidyr::pivot_wider(names_from = PL36_char,                  # Create indicator columns
                         values_from = PL36_char,
                         values_fn = length,
                         names_prefix = "PL36_") %>%
      select(-row_id) %>%  # Drop helper column
      mutate(across(starts_with("PL36_"),
                    # If PL01R did occur and the current column isn't the event that occurred, set to zero
                    ~ ifelse(!is.na(PL36) & cur_column() != paste0("PL36_", PL36), 0, .))) %>%
      # Remove original
      dplyr::select(-any_of(c("PL35", "PL36", "PL36_NA", "PL36_0"))) %>%

      # In 2021/2022, Illness and accident are asked separately. Make into one category to match previous years
      rename_with(
        .fn = ~ paste0(., "_occurred"),          # Append "_occurred" to matched column names
        .cols = matches("^PL01R|^PL06|^PL11|^PL16|^PL21|^PL26|^PL36") # Match prefixes
      ) %>%
      # Add suffering variables (not used)
      dplyr::rename(
        PL06_suffering = PL10,
        PL11_suffering = PL15,
        PL16_suffering = PL20,
        PL21_suffering = PL25,
        PL26_suffering = PL30) %>%
      # Add suffering variable for illness/accident (not used)
      mutate(PL01R_1_suffering = PL05,
             PL01R_2_suffering = PL05,
             PL01R_3_suffering = PL05,
             PL01R_4_suffering = PL05,
             PL01R_5_suffering = PL05,
             PL01R_6_suffering = PL05,
             PL01R_7_suffering = PL05) %>%
      dplyr::select(-PL05) %>%
      tidyr::pivot_longer(cols = matches("^PL01R|^PL06|^PL11|^PL16|^PL21|^PL26|^PL36"),
                          # names_pattern = c("(.*)_suffering(.*)"), names_to = c(".value", "event")
                          names_to = c("event_code", ".value"),            # Split column names into "event" and value type
                          # names_sep = "_"
                          names_pattern = "^(.*)_(suffering|occurred)$"                 # Match event and suffix

      ) %>%

      mutate_at(c("occurred", "suffering"), ~ as.numeric(as.character(.))) %>%

      # Remove person-year if all life events are NA
      group_by(IDPERS, YEAR) %>%
      mutate(all_events_missing = sum(is.na(occurred)) == length(occurred)) %>%
      ungroup() %>%
      filter(!all_events_missing) %>% select(-all_events_missing) %>%

      # Change event names
      mutate(event = dplyr::recode_factor(event_code, !!!tibble::deframe(event_dict %>% dplyr::select(var_name, recode_var)))) %>%
      # For correspondence with HILDA dataset, rename variables
      dplyr::rename(p_id = IDPERS, crosswave_h_id = IDHOUS, sex = SEX, age = AGE) %>%
      mutate(wave_nr = as.numeric(as.character(YEAR)) - min(as.numeric(as.character(YEAR)), na.rm = TRUE) + 1,
             ses = NA) %>%
      arrange(p_id, wave_nr, event_code)

    # Save dataframe
    saveRDS(df_binary, filepath_df_binary)
  }
  return(NULL)
}

# Prepare HILDA data
prepare_HILDA = function(filepath_dataset, filepath_df_binary, event_dict, rerun){

  # If interested: Create dictionary of subjective variables (not used)
  subj_dict = matrix(c(
    c("lssupac", "I don't have anyone that I can confide in", "no_confidant", "negative"),
    c("lssupcd", "There is someone who can always cheer me up when I'm down", "always_someone_to_cheer_me_up", "positive"),
    c("lssuplf", "I seem to have a lot of friends", "lot_of_friends", "positive"),
    c("lssuplt", "I have no one to lean on in times of trouble", "no_one_to_lean_on", "negative"),
    c("lssupnh", "I often need help from other people but can't get it", "need_help_from_people_but_cant_get_it", "negative"),
    c("lssuppi", "I enjoy the time I spend with the people who are important to me", "enjoy_time_with_important_people", "positive"),
    c("lssuppv", "People dont come to visit me as often as I would like", "people_dont_visit_often", "negative"),
    c("lssupsh", "When I need someone to help me out, I can usually find someone", "people_help_when_I_need_it", "positive"),
    c("lssuptp", "When somethings on my mind, just talking with the people I know can make me feel better", "talking_with_people_I_know_helps", "positive"),
    c("lssupvl", "I often feel very lonely", "feel_lonely","negative"),
    c("lsrelpc", "Satisfaction with: Partners relationship with children", "satisfaction_partners_relationship_with_children", "positive"),
    c("lsrelrp", "Satisfaction with: Relationship with parents", "satisfaction_relationship_with_parents", "positive"),
    c("lsrelrs", "Satisfaction with: Relationship with step parents", "satisfaction_relationship_with_step_parents", "positive"),
    c("lsrelsc", "Satisfaction with: Children", "satisfaction_with_children", "positive"),
    c("lsrelsp", "Satisfaction with: Partner", "satisfaction_with_partner", "positive"),
    c("lsrelst", "Satisfaction with: Relationship with step children", "satisfaction_relationship_with_stepchildren", "positive"),
    c("lssocal", "How often get together socially with friends/relatives
  not living with you", "frequency_get_together_friends_relatives", "positive"),
    c("losat", "How satisfied are you with your life", "life_satisfaction", "positive"),
    c("losateo", "Satisfaction - Your employment opportunities", "satisfaction_employment_opportunities", "positive"),
    c("losatfs","Satisfaction - Your financial situation", "satisfaction_financial_situation", "positive"),
    c("losatft", "Satisfaction - The amount of free time you have", "satisfaction_amount_free_time", "positive"),
    c("losathl", "Satisfaction - The home in which you live", "satisfaction_home", "positive"),
    c("losatlc","Satisfaction - Feeling part of your local community", "satisfaction_part_of_local_community", "positive"),
    c("losatnl", "Satisfaction - The neighbourhood in which you live", "satisfaction_neighbourhood", "positive"),
    c("losatsf", "Satisfaction - How safe you feel", "satisfaction_feeling_safe", "positive"),
    c("losatyh", "Satisfaction - Your health", "health_satisfaction", "positive")
  ),
  ncol = 4, byrow = TRUE ) %>%
    magrittr::set_colnames(c("var_name", "description", "recode_var", "valence")) %>% as.data.frame() %>%
    mutate_all(~ stringr::str_replace(.x, ": yes, no", ""))

  # Re-run creation of dataframe over all waves
  if (rerun | !file.exists(filepath_df_binary)){

    df_list = list()

    # Loop over waves
    for (wave in letters[2:23]){ #  Exclude wave 1 which didn't measure life events

      print(wave)
      folder = ifelse(wave %in% letters[1:11],
                      "2. SPSS 230c (Zip file 1 of 4 - Combined Data Files a-k)",
                      "2. SPSS 230c (Zip file 2 of 4 - Combined Data Files l-w)")
      df = haven::read_sav(file.path(filepath_dataset, sprintf("%s/Combined %s230c.sav", folder, wave)))

      # Recode entries with attributes of that column
      recode_col <- function(x, col_name, df) {

        recode_vec <- attributes(df[[col_name]])$labels

        # Reverse the name-value pairs
        dplyr::recode(x, !!!tibble::deframe(tibble::enframe(recode_vec) %>% dplyr::select(value, name)))
      }

      # Get events
      df1 = df %>% dplyr::select(xwaveid,
                                 all_of(sprintf("%shhhqivw", wave)),
                                 matches(sprintf("^%shhs3add", wave)),
                                 matches(sprintf("^%shgsex$", wave)),
                                 matches(sprintf("^%shgage$", wave)),
                                 matches(sprintf("^%shhrhid$", wave)),
                                 matches(sprintf("^%sle", wave))) %>%

        # Recode entries with attributes
        mutate(across(matches(sprintf("^%sle|^%shgsex$", wave, wave)),
                      ~ recode_col(as.numeric(as.character(.x)), cur_column(), df ))
        ) %>%

        # Remove wave prefix from all variables
        dplyr::rename_with(~ stringr::str_replace_all(.x, sprintf("^%s", wave), "")
        ) %>%
        dplyr::rename(sex = hgsex) %>%
        dplyr::rename(age = hgage) %>%
        dplyr::rename(h_id = hhrhid) %>%
        dplyr::rename(ses = hhs3add) %>%
        mutate(ses = as.numeric(as.character(ses))) %>%

        # Prepare column names for pivot_longer
        dplyr::rename_with(~ stringr::str_replace_all(.x, "na$",
                                                      "_na") ) %>%
        dplyr::rename_with(~ stringr::str_replace_all(.x, "q([0-9])$",
                                                      "_q\\1") ) %>%
        dplyr::rename_with(~ ifelse(!grepl("_", .x), paste0(.x, "_occurred"), .x), matches("^le") ) %>%

        # For each event, create row with columns indicating timing of event
        tidyr::pivot_longer(cols = matches("^le"),
                            names_pattern = c("(.*)_(.*)"),
                            names_to = c("event_code", ".value")) %>%
        as.data.frame() %>%
        # Make occurred binary 0, 1, or NA; note that the column "na" refers to no answer regarding when it happened
        dplyr::rename(na_when = na) %>%
        mutate(across(c(occurred, q1, q2, q3, q4, na_when), ~ ifelse(tolower(.x) == "no", 0, ifelse(tolower(.x) == "yes", 1, .x)))) %>%
        # Remove participant in this wave in case ALL events happened (e.g. subject "0101523" at wave "g" had all events happen in the first quarter, unlikely)
        dplyr::filter(sum(occurred == 1, na.rm = TRUE) != length(occurred), .by = xwaveid) %>%
        # Remove participant in this wave if all life events are NA
        filter(sum(!is.na(occurred)) > 0, .by = xwaveid)

      # Get subjective variables (not used)
      df2 = df %>% dplyr::select(xwaveid,
                                 matches(sprintf("^%s%s", wave, subj_dict$var_name))
      ) %>%
        # Remove wave prefix from all variables
        dplyr::rename_with(~ stringr::str_replace_all(.x, sprintf("^%s", wave), "")
        ) %>%
        dplyr::rename_with(~ tibble::deframe(subj_dict %>% dplyr::select(var_name, recode_var)),
                           matches(sprintf("^%s", subj_dict$var_name))
        ) %>%
        mutate(across(matches(subj_dict$recode_var),
                      ~ as.numeric(as.character(.x)))
        )

      # Merge subjective with event occurrence
      df3 = merge(df1, df2, all.x = TRUE) %>% mutate(wave = !!wave)

      df_list[[wave]] = df3
      rm(df)
      rm(df1)
      rm(df2)
      rm(df3)
    }
    rm(wave) # Remove variable wave

    # Check dataframes
    purrr::map(df_list, ncol)
    purrr::map(df_list, nrow)
    purrr::map(df_list, head, n = 1)
    purrr::map(df_list, function(df){df %>% pull(wave) %>% unique})

    # Check whether there is more than 1 entry per person per wave
    check_func = function(df){
      df %>% group_by(xwaveid, wave) %>%
        dplyr::summarise(nr_entries = n(), l = length(unique(life_satisfaction)), .groups = 'drop') %>%
        filter(l != 1) %>% as.data.frame()
    }
    check = purrr::map(df_list, check_func)
    stopifnot(sum(purrr::map(check, nrow) %>% unlist()) == 0) # Should all be dataframes with 0 rows

    # Merge files
    df_binary_ = do.call(dplyr::bind_rows, df_list) %>% as.data.frame() %>%
      arrange(xwaveid, wave) %>%
      # Recode event names
      dplyr::mutate(event = recode_factor(event_code,
                                          !!!tibble::deframe(event_dict %>% dplyr::select(var_name, recode_var) )  )) %>%
      # Interview date
      mutate(hhhqivw = as.Date(hhhqivw, format="%d/%m/%Y")) %>%
      mutate(YEAR = as.numeric(format(hhhqivw, "%Y")),
             interview_YEAR = as.numeric(format(hhhqivw, "%Y")),
             interview_month = as.numeric(format(hhhqivw, "%m")),
             interview_day = as.numeric(format(hhhqivw, "%d"))) %>%
      mutate(across(matches("^q[0-9]"), ~ as.numeric(as.character(.x)))) %>%
      # Get number of repetitions per event and the beginning and end quarter
      mutate(nr_rep_event = dplyr::select(., matches("^q[0-9]")) %>% rowSums(na.rm = TRUE),
             begin_q = ifelse(occurred == 0, NA, ifelse(q1 == 1, 1, ifelse(q2 == 1, 2, ifelse(q3 == 1, 3, ifelse(q4 == 1, 4, NA))))),
             end_q = ifelse(occurred == 0, NA, ifelse(q1 == 1 & q2 != 1, 1,
                                                      ifelse(q2 == 1 & q3 != 1, 2,
                                                             ifelse(q3 == 1 & q4 != 1, 3,
                                                                    ifelse(q4 == 1 & q3 == 0 & q2 == 0 & q1 == 0, 4,
                                                                           ifelse(nr_rep_event == 4, 4, NA )))))  )
      ) %>%
      # Transform year-quarter to window of months, i.e. Q1 = 0-3 months; Q2 = 4-6 months; Q3 = 7-9 months; Q4 = 10-12 months ago (not used)
      mutate(begin_q = ifelse(is.na(begin_q), NA, ((begin_q - 1) * 3 + 1) %>%
                                ifelse(. == 1, 0, .) ), # Adjust because Q1 = 0-3 months, not 1-3 months
             end_q = ifelse(is.na(end_q), NA, (end_q - 1) * 3 + 3),
             begin_window_event = lubridate::`%m-%`(hhhqivw, months(end_q)),
             end_window_event = lubridate::`%m-%`(hhhqivw, months(begin_q))
      ) %>%
      # For correspondence with SHP data, add suffering variable and rename subject ID
      mutate(suffering = NA) %>%
      dplyr::rename(IDPERS = xwaveid) %>%
      # Convert wave letter to numeric
      mutate(wave_nr = match(wave, letters)) %>%
      dplyr::rename(p_id = IDPERS)

    head(df_binary_) %>% as.data.frame()

    # Create crosswave household identifier: As an end product, we want for each xwaveid and wave the crosswave household identifier
    # Read Master file
    df_master = haven::read_sav(file.path(filepath_dataset, sprintf("2. SPSS 230c (Zip file 4 of 4 - Eperson and Other Data Files)/Master w230c.sav")))

    df_h_ids = df_master %>%
      dplyr::rename(p_id = xwaveid) %>%
      dplyr::select(p_id, matches("hhrhid$")) %>%
      mutate(across(!p_id, as.character)) %>%
      tidyr::pivot_longer(cols = !p_id,
                          names_pattern = c("(.*)hhrhid"),
                          names_to = c("wave"), values_to = "h_id") %>%
      # Delete empty rows
      filter(h_id != "") %>%
      # Household IDs may be the same across waves, but these are not the same indivduals - ensure unique household ID by appending wave letter
      mutate(h_id = paste0(wave, "_", h_id)) %>%
      group_by(h_id, wave) %>%
      dplyr::summarise(unique_p_ids = sort(unique(p_id)) %>% paste0(collapse=","), .groups = 'drop') %>%
      group_by(unique_p_ids) %>%
      mutate(crosswave_h_id = cur_group_id() %>% as.character()) %>% ungroup() %>%
      arrange(crosswave_h_id) %>%
      # Split the unique_p_ids into separate rows
      tidyr::separate_rows(unique_p_ids, sep = ",") %>%
      dplyr::select(crosswave_h_id, wave, unique_p_ids) %>%
      dplyr::rename(p_id = unique_p_ids)

    df_h_ids

    # Merge dataframes
    df_binary = df_binary_ %>% merge(df_h_ids, all.x = TRUE) %>%
      # Mutate age to numeric, haven-labelled otherwise
      mutate(age = as.numeric(as.character(factor(age))))
    df_binary %>% head

    saveRDS(df_binary, filepath_df_binary)
    rm(df_list)
  }
  return(NULL)
}


# Create dataframes for analyses
get_analysis_df = function(name_dataset, df_binary, event_dict, filepath_deriv, rerun){

  # Define filepaths
  filepath_df_per_event_pp_py = file.path(filepath_deriv, sprintf("df_per_event_pp_py.RDS"))
  filepath_df_per_event = file.path(filepath_deriv, sprintf("df_per_event.RDS"))
  filepath_df_nr_neg_events = file.path(filepath_deriv, sprintf("df_nr_neg_events.RDS"))
  filepath_df_nr_negevents_pp_py = file.path(filepath_deriv, sprintf("df_nr_negevents_pp_py.RDS"))
  filepath_df_nr_negevents_pp = file.path(filepath_deriv, sprintf("df_nr_negevents_pp.RDS"))

  # Rerun creation of analysis dataframes
  if (rerun | any(!sapply(c(filepath_df_per_event_pp_py,
                           filepath_df_per_event,
                           filepath_df_nr_neg_events,
                           filepath_df_nr_negevents_pp_py,
                           filepath_df_nr_negevents_pp
                           ), file.exists) ) ){

    # Remove the "other" event category in SHP as it was an open question, categorized post-hoc, and the categories are not very populated
    if (name_dataset == "SHP"){
      df_binary = df_binary %>%
        filter(!grepl("^PL36", event_code))
    }

    # Dataframe of event occurrences per event per person per year
    df_per_event_pp_py = df_binary %>%
      dplyr::select(p_id, age, sex, crosswave_h_id, wave_nr, event, event_code, occurred) %>%
      # Only adults
      dplyr::filter(age >= 18) %>%
      group_by(p_id, age, sex, crosswave_h_id, wave_nr, event, event_code) %>%

      # To be sure, we use summarise here, as there may be multiple occurrence in the dataset in case of HILDA
      dplyr::summarise(nr_occur = sum(occurred == 1, na.rm = TRUE),
                       nr_nooccur = sum(occurred == 0, na.rm = TRUE),
                       nr_missing_occur = sum(is.na(occurred)),
                       .groups = 'drop') %>%
      dplyr::arrange(p_id, wave_nr, event_code) %>%
      group_by(p_id) %>%

      # Summarise when they were observed
      dplyr::mutate(
        first_obs_wave_nr_pp = min(wave_nr, na.rm  = TRUE),
        last_obs_wave_nr_pp = max(wave_nr, na.rm  = TRUE),
        timespan_obs_pp = last_obs_wave_nr_pp - first_obs_wave_nr_pp + 1,
        nr_years_obs = length(unique(wave_nr)),
        nr_years_missing_in_timespan = timespan_obs_pp - nr_years_obs,
      ) %>%
      ungroup() %>%
      mutate(valence = dplyr::recode(event,
                                     !!!tibble::deframe(event_dict %>%
                                                          dplyr::select(recode_var, valence)))) %>%
      mutate(dependence = dplyr::recode(event,
                                        !!!tibble::deframe(event_dict %>%
                                                             dplyr::select(recode_var, in_dependent)))) %>%
      arrange(p_id, wave_nr) %>%
      # Remove years where ALL negative events were missing
      group_by(p_id, wave_nr) %>%
      dplyr::filter(!(sum(nr_occur[valence == "negative"]) == 0 &
                        sum(nr_nooccur[valence == "negative"]) == 0)) %>%
      dplyr::ungroup()

    # Should have only one row per person, wave, event combination
    stopifnot(nrow(df_per_event_pp_py) > 0)
    stopifnot(nrow(df_per_event_pp_py %>% group_by(p_id, wave_nr, event) %>%
                     dplyr::summarise(n = n(), .groups = 'drop') %>% filter(n > 1)) == 0)

    test = df_per_event_pp_py %>% select(p_id, sex) %>%
      group_by(p_id) %>%
      # select first sex
      slice(1) %>% ungroup() %>%
      # distinct() %>%
      pull(sex) %>% table()
    test / sum(test) * 100

    # Dataframe of event occurrences per event
    df_per_event = df_per_event_pp_py %>%
      group_by(event_code, event, valence, dependence) %>%
      dplyr::summarise(nr_occur = sum(nr_occur),
                       nr_nooccur = sum(nr_nooccur),
                       nr_missing_occur = sum(nr_missing_occur),
                       .groups ='drop')
    df_per_event %>% arrange(desc(nr_occur)) %>% as.data.frame()

    # Number of individuals, households, and person-years
    test = df_per_event_pp_py %>% filter(valence == "negative")
    test %>% dplyr::select(p_id, wave_nr) %>% distinct() %>% nrow()

    # Dataframe of number of NEGATIVE event occurrences per person per year
    df_nr_negevents_pp_py = df_per_event_pp_py %>%
      filter(valence == "negative") %>%
      group_by(p_id, crosswave_h_id, age, sex, wave_nr) %>%

      # Use summarise to collapse across event dependence
      dplyr::summarise(nr_occur = sum(nr_occur),
                       nr_nooccur = sum(nr_nooccur),
                       nr_missing_occur = sum(nr_missing_occur),
                       .groups = 'drop') %>%

      # If you want to remove years with ANY missing events, use:
      # dplyr::filter(nr_missing_occur == 0) %>%

      # Add t_id -> time id for OBSERVED years, should always be consecutive
      dplyr::arrange(p_id, wave_nr) %>%
      group_by(p_id) %>%
      dplyr::mutate(t_id = seq.int(n())) %>%

      # Cumulative number of events
      dplyr::mutate(nr_occur_cumsum = cumsum(nr_occur)) %>%

      # Summarise when they were observed
      dplyr::mutate(
        first_obs_wave_nr_pp = min(wave_nr, na.rm  = TRUE),
        last_obs_wave_nr_pp = max(wave_nr, na.rm  = TRUE),
        timespan_obs_pp = last_obs_wave_nr_pp - first_obs_wave_nr_pp + 1,
        nr_years_obs = length(unique(wave_nr)),
        nr_years_missing_in_timespan = timespan_obs_pp - nr_years_obs,
      ) %>%
      ungroup() %>% arrange(p_id, wave_nr)

    # Should be an empty dataframe, only want one entry per person per year
    stopifnot(nrow(df_nr_negevents_pp_py %>% group_by(p_id, wave_nr) %>%
                     dplyr::summarise(n = n(), .groups = 'drop') %>%
                     filter(n > 1)) == 0)

    # Dataframe of overall number of negative event occurrences in the dataset
    df_nr_neg_events = df_nr_negevents_pp_py %>%
      filter(!is.na(t_id)) %>%
      group_by(nr_occur) %>%
      dplyr::summarise(n = n(), .groups = 'drop') %>%
      dplyr::mutate(n_norm = n / sum(n))

    # Dataframe of number of negative event occurrences per person
    df_nr_negevents_pp = df_nr_negevents_pp_py %>%
      group_by(p_id, first_obs_wave_nr_pp, last_obs_wave_nr_pp,
               timespan_obs_pp, nr_years_obs, nr_years_missing_in_timespan) %>%
      dplyr::summarise(nr_occur = sum(nr_occur, na.rm = TRUE),
                       nr_nooccur = sum(nr_nooccur, na.rm = TRUE),
                       nr_missing_occur = sum(nr_missing_occur, na.rm = TRUE),
                       .groups = 'drop'
      )
    df_nr_negevents_pp %>% head(n=100) %>% as.data.frame()

    # Save all dataframes
    saveRDS(df_per_event_pp_py, filepath_df_per_event_pp_py)
    saveRDS(df_per_event, filepath_df_per_event)
    saveRDS(df_nr_neg_events, filepath_df_nr_neg_events)
    saveRDS(df_nr_negevents_pp_py, filepath_df_nr_negevents_pp_py)
    saveRDS(df_nr_negevents_pp, filepath_df_nr_negevents_pp)
  }

  df_per_event_pp_py = readRDS(filepath_df_per_event_pp_py)
  df_per_event = readRDS(filepath_df_per_event)
  df_nr_neg_events = readRDS(filepath_df_nr_neg_events)
  df_nr_negevents_pp_py = readRDS(filepath_df_nr_negevents_pp_py)
  df_nr_negevents_pp = readRDS(filepath_df_nr_negevents_pp)

  return(list(df_per_event_pp_py = df_per_event_pp_py,
              df_per_event = df_per_event,
              df_nr_neg_events = df_nr_neg_events,
              df_nr_negevents_pp_py = df_nr_negevents_pp_py,
              df_nr_negevents_pp = df_nr_negevents_pp))
}
