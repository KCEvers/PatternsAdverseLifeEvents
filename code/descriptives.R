##### DESCRIPTIVE ANALYSIS OF LIFE EVENTS - REFACTORED #####

#' Main wrapper function for descriptive analysis
#' @param name_dataset Dataset name (e.g., "SHP", "HILDA")
#' @param filepath_base Base file path
#' @param P Parameters list
#' @param rerun Whether to rerun analysis
run_descriptives <- function(name_dataset, filepath_base, P, rerun = FALSE) {

  cat(sprintf("\n=== Running descriptive analysis for %s ===\n", name_dataset))

  # 1. Prepare data
  cat("\n=== Preparing data ===\n")
  event_dict <- P$event_dicts[[name_dataset]]
  datalist <- prepare_data(name_dataset, filepath_base, event_dict, rerun = rerun)

  # 2. Calculate frequency statistics
  cat("\n=== Calculating frequencies ===\n")
  df_freq <- calculate_event_frequencies(datalist, name_dataset, event_dict, P)

  # 3. Calculate numerical descriptives
  cat("\n=== Calculating numerical descriptives ===\n")
  num_desc <- calculate_numerical_descriptives(datalist, name_dataset, event_dict)

  # 4. Create plots
  cat("\n=== Creating plots ===\n")
  plots <- create_all_plots(datalist, df_freq, name_dataset, P)

  # 5. Save plots
  cat("\n=== Saving plots ===\n")
  save_all_plots(plots, datalist, name_dataset)

  return(list(
    datalist = datalist,
    df_freq = df_freq,
    num_desc = num_desc,
    plots = plots
  ))
}


# ============================================================================
# DATA PREPARATION AND FREQUENCIES
# ============================================================================

#' Calculate event frequencies with ordering
#' @param datalist Prepared data list
#' @param name_dataset Dataset name
#' @param event_dict Event dictionary
#' @param P Parameters list
calculate_event_frequencies <- function(datalist, name_dataset, event_dict, P) {

  # Calculate basic frequencies
  df_freq <- datalist$df_per_event %>%
    dplyr::filter(valence == "negative") %>%
    dplyr::mutate(
      dataset = name_dataset,
      perc_occur = nr_occur / (nr_occur + nr_nooccur) * 100,
      perc_occur_incl_missing = nr_occur / (nr_occur + nr_nooccur + nr_missing_occur) * 100,
      sum_n = nr_occur + nr_nooccur + nr_missing_occur
    )

  # Order events by frequency
  df_freq <- df_freq %>%
    dplyr::mutate_at(
      "event",
      ~ dplyr::recode_factor(
        .x,
        !!!tibble::deframe(
          event_dict %>%
            dplyr::select(recode_var, description) %>%
            dplyr::slice(match(
              df_freq %>%
                dplyr::arrange(perc_occur_incl_missing) %>%
                dplyr::pull(event) %>%
                as.character(),
              recode_var
            ))
        ),
        .ordered = TRUE
      )
    )

  return(df_freq)
}


#' Calculate numerical descriptive statistics
#' @param datalist Prepared data list
#' @param name_dataset Dataset name
#' @param event_dict Event dictionary
calculate_numerical_descriptives <- function(datalist, name_dataset, event_dict) {

  res <- data.frame(dataset = name_dataset)

  # Top two most frequent events
  df_events <- datalist$df_per_event %>%
    dplyr::filter(valence == "negative") %>%
    dplyr::arrange(desc(nr_occur)) %>%
    dplyr::mutate(
      event = as.character(event),
      nr_occur_norm = nr_occur / sum(nr_occur)
    )

  res$top1_event <- df_events$event[1]
  res$top2_event <- df_events$event[2]
  res$total_nr_neg_events <- sum(df_events$nr_occur)
  res$perc_top_two_neg_events <- sprintf(
    "The top two events were responsible for %.2f%% of all events",
    sum(df_events %>% dplyr::slice(1:2) %>% dplyr::pull(nr_occur)) / res$total_nr_neg_events * 100
  )

  # Number of person-years < 18
  res$nr_personyears_underage <- datalist$df_binary %>%
    dplyr::filter(age < 18) %>%
    dplyr::select(p_id, wave_nr) %>%
    dplyr::distinct() %>%
    nrow()

  # Number of person-years missing all adverse events
  res$nr_all_negevents_missing <- datalist$df_binary %>%
    dplyr::select(p_id, age, sex, crosswave_h_id, wave_nr, event, event_code, occurred) %>%
    dplyr::filter(age >= 18) %>%
    dplyr::group_by(p_id, age, sex, crosswave_h_id, wave_nr, event, event_code) %>%
    dplyr::summarise(
      nr_occur = sum(occurred == 1, na.rm = TRUE),
      nr_nooccur = sum(occurred == 0, na.rm = TRUE),
      nr_missing_occur = sum(is.na(occurred)),
      .groups = 'drop'
    ) %>%
    dplyr::arrange(p_id, wave_nr, event_code) %>%
    dplyr::mutate(
      valence = dplyr::recode(event, !!!tibble::deframe(event_dict %>% dplyr::select(recode_var, valence))),
      dependence = dplyr::recode(event, !!!tibble::deframe(event_dict %>% dplyr::select(recode_var, in_dependent)))
    ) %>%
    dplyr::arrange(p_id, wave_nr) %>%
    dplyr::group_by(p_id, wave_nr) %>%
    dplyr::filter(
      sum(nr_occur[valence == "negative"]) == 0 &
        sum(nr_nooccur[valence == "negative"]) == 0
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(p_id, wave_nr) %>%
    dplyr::distinct() %>%
    nrow()

  # Sample demographics
  temp <- datalist$df_per_event_pp_py %>%
    dplyr::select(p_id, crosswave_h_id, age) %>%
    dplyr::distinct()

  res$nr_p_id_overall <- nrow(temp)
  res$nr_h_id_overall <- length(unique(temp$crosswave_h_id))

  # Co-occurrence statistics
  df_cooccur <- datalist$df_nr_negevents_pp_py %>%
    dplyr::filter(!is.na(t_id)) %>%
    dplyr::group_by(nr_occur) %>%
    dplyr::summarise(n = n(), .groups = 'drop') %>%
    dplyr::mutate(n_norm = n / sum(n))

  res$perc_no_event <- df_cooccur %>% dplyr::filter(nr_occur == 0) %>% dplyr::pull(n_norm)
  res$perc_one_event <- df_cooccur %>% dplyr::filter(nr_occur == 1) %>% dplyr::pull(n_norm)
  res$perc_cooccur <- df_cooccur %>% dplyr::filter(nr_occur > 1) %>% dplyr::pull(n_norm) %>% sum()

  # Co-occurrence conditional on having at least one event
  df_cooccur_cond <- df_cooccur %>%
    dplyr::filter(nr_occur > 0) %>%
    dplyr::mutate(n_norm = n / sum(n))

  res$perc_one_event_if_event <- df_cooccur_cond %>% dplyr::filter(nr_occur == 1) %>% dplyr::pull(n_norm)
  res$perc_cooccur_if_event <- df_cooccur_cond %>% dplyr::filter(nr_occur > 1) %>% dplyr::pull(n_norm) %>% sum()

  return(res)
}


# ============================================================================
# PLOT CREATION
# ============================================================================

#' Create all descriptive plots
#' @param datalist Prepared data list
#' @param df_freq Frequency dataframe
#' @param name_dataset Dataset name
#' @param P Parameters list
create_all_plots <- function(datalist, df_freq, name_dataset, P) {

  # Set plot parameters
  plot_params <- get_plot_parameters(name_dataset)

  plots <- list(
    pl_negevent_abs_freq = plot_absolute_event_frequency(df_freq, name_dataset, P, plot_params),
    pl_negevent_freq = plot_relative_event_frequency(df_freq, name_dataset, P, plot_params),
    pl_nr_events = plot_yearly_event_distribution(datalist, P),
    pl_nr_events_py = plot_event_distribution_per_year(datalist, name_dataset, P),
    pl_event_distr_over_time = plot_event_distribution_over_time(datalist, name_dataset, P),
    pl_cum_distr = plot_cumulative_distribution(datalist, name_dataset, P),
    # pl_age = plot_age_distribution(datalist, name_dataset, P),
    pl_demo = plot_demo_distribution(datalist, name_dataset, P),
    pl_nr_years_obs = plot_observation_years(datalist, name_dataset, P)
  )

  return(plots)
}


#' Get plot-specific parameters based on dataset
#' @param name_dataset Dataset name
get_plot_parameters <- function(name_dataset) {

  if (name_dataset == "HILDA") {
    return(list(
      height = 160,
      max_char = 25,
      text_size = 10,
      freq_offset = 4000,
      perc_offset = 0.02 * 100
    ))
  } else if (name_dataset == "SHP") {
    return(list(
      height = 160,
      max_char = 25,
      text_size = 10,
      freq_offset = 8000,
      perc_offset = 0.03 * 100
    ))
  }
}


#' Plot absolute event frequency
plot_absolute_event_frequency <- function(df_freq, name_dataset, P, plot_params) {

  pl <- df_freq %>%
    ggplot() +
    geom_bar(
      aes(x = nr_occur, y = event),
      col = P$col_data,
      fill = P$col_data,
      stat = 'identity',
      alpha = 0.5
    ) +
    geom_text(
      aes(x = nr_occur + plot_params$freq_offset, y = event, label = nr_occur),
      col = 'grey30',
      family = P$font_family,
      size = 4
    ) +
    ggh4x::facet_grid2(. ~ dataset) +
    P$own_theme +
    theme(
      axis.text = element_text(size = plot_params$text_size + 2),
      plot.margin = unit(c(0.01, 0.9, 0.01, 0.01), 'cm'),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    scale_y_discrete(
      labels = function(x) sapply(x, wrap_equally, max_chars = 31),
      expand = c(0, 0)
    ) +
    scale_x_continuous(
      n.breaks = 4,
      expand = expansion(mult = c(0, 0.075))
    ) +
    labs(
      x = "Frequency in person-years",
      y = "",
      title = "Event frequency\nacross people and years"
    )

  return(pl)
}


#' Plot relative event frequency
plot_relative_event_frequency <- function(df_freq, name_dataset, P, plot_params) {

  # Calculate reference lines for 10,000 instances
  x_vals <- seq(10000, ifelse(name_dataset == "SHP", 50000, 40000), by = 10000)
  perc_10000_instances <- data.frame(
    x = x_vals,
    perc = x_vals / median(df_freq$sum_n) * 100
  )

  pl <- df_freq %>%
    ggplot() +
    geom_vline(
      data = perc_10000_instances,
      aes(xintercept = perc),
      linetype = 'dashed',
      color = 'grey80'
    ) +
    geom_bar(
      aes(x = perc_occur_incl_missing, y = event),
      col = P$col_data,
      fill = P$col_data,
      stat = 'identity',
      alpha = 0.6
    ) +
    geom_text(
      aes(x = perc_occur_incl_missing + plot_params$perc_offset, y = event,
          label = sprintf("%.2f%%", perc_occur_incl_missing)),
      col = 'grey30',
      family = P$font_family,
      size = 4
    ) +
    ggh4x::facet_grid2(. ~ dataset) +
    P$own_theme +
    theme(
      axis.text = element_text(size = plot_params$text_size + 2),
      title = element_text(size = plot_params$text_size + 6),
      plot.margin = unit(c(0.01, 0.5, 0.01, 0.01), 'cm'),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    scale_y_discrete(
      labels = function(x) sapply(x, wrap_equally, max_chars = plot_params$max_char),
      expand = c(0, 0)
    ) +
    scale_x_continuous(
      n.breaks = 4,
      expand = expansion(mult = c(0, 0.1))
    ) +
    labs(
      x = "Frequency (%) of events across all person-years",
      y = "",
      title = "Event frequency\nacross people and years"
    )

  return(pl)
}


#' Plot yearly event distribution (overall)
plot_yearly_event_distribution <- function(datalist, P) {

  pl <- datalist$df_nr_negevents_pp_py %>%
    dplyr::filter(!is.na(t_id)) %>%
    dplyr::group_by(nr_occur) %>%
    dplyr::summarise(n = n(), .groups = 'drop') %>%
    dplyr::mutate(n_norm = n / sum(n)) %>%
    ggplot() +
    geom_bar(aes(x = nr_occur, y = n), stat = 'identity') +
    P$own_theme +
    theme(panel.grid.minor = element_blank()) +
    scale_y_continuous(
      n.breaks = 5,
      expand = expansion(mult = c(0, 0.01))
    ) +
    scale_x_continuous(
      n.breaks = max(datalist$df_nr_negevents_pp_py$nr_occur),
      expand = c(0, 0)
    ) +
    labs(
      x = "Number of events per year",
      y = "Frequency in person-years"
    )

  return(pl)
}


#' Plot event distribution per observation year
plot_event_distribution_per_year <- function(datalist, name_dataset, P) {

  pl <- datalist$df_nr_negevents_pp_py %>%
    dplyr::filter(!is.na(t_id)) %>%
    dplyr::group_by(t_id, nr_occur) %>%
    dplyr::summarise(n = n(), .groups = 'drop') %>%
    dplyr::group_by(t_id) %>%
    dplyr::mutate(
      n_norm = n / sum(n),
      t_id_label = sprintf("Year %02d, n = %d", t_id, sum(n)) %>%
        stringr::str_wrap(width = 9)
    ) %>%
    dplyr::ungroup() %>%
    ggplot() +
    geom_bar(aes(x = nr_occur, y = n_norm), stat = 'identity') +
    P$own_theme +
    theme(
      axis.text.y = element_text(size = 7),
      plot.title = element_text(hjust = 0.5),
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), 'cm'),
      text = element_text(family = "serif"),
      title = element_text(family = "serif", size = 12),
      plot.subtitle = element_text(hjust = 0.5),
      panel.spacing.y = unit(0.4, "lines"),
      strip.background = element_rect(colour = P$col_facets, fill = P$fill_facets),
      strip.text.y = element_text(angle = 0)
    ) +
    scale_x_continuous(n.breaks = 4, expand = c(0, 0)) +
    labs(
      title = "Empirical distribution of\nnumber of adverse life events per year",
      subtitle = name_dataset,
      x = "Number of events per year",
      y = "Frequency (normalized per year)"
    ) +
    ggh4x::facet_wrap2(
      t_id_label ~ .,
      ncol = 2,
      dir = 'v',
      strip.position = 'right'
    )

  return(pl)
}


#' Plot event distribution over time (area plot)
plot_event_distribution_over_time <- function(datalist, name_dataset, P) {

  pl <- datalist$df_nr_negevents_pp_py %>%
    dplyr::filter(!is.na(t_id)) %>%
    dplyr::group_by(t_id, nr_occur) %>%
    dplyr::summarise(n = n(), .groups = 'drop') %>%
    dplyr::group_by(t_id) %>%
    dplyr::mutate(n_norm = n / sum(n)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      nr_occur = factor(nr_occur %>% as.character(), levels = as.character(0:13)),
      dataset = name_dataset
    ) %>%
    ggplot() +
    geom_area(
      aes(x = t_id, y = n_norm, fill = nr_occur),
      col = 'grey30',
      linewidth = 0.8,
      position = 'stack'
    ) +
    ggh4x::facet_grid2(. ~ dataset) +
    viridis::scale_fill_viridis(
      discrete = TRUE,
      direction = -1,
      name = "Number of events",
      drop = FALSE
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = expansion(mult = c(0.002, 0.002))) +
    labs(x = "Observation year", y = "Percentage (normalized per year)") +
    P$own_theme +
    theme(
      legend.position = 'bottom',
      plot.margin = unit(c(1, 1, 1, 1), "cm"),
      axis.text = element_text(size = 10),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 12),
      axis.title = element_text(size = 14),
      strip.text = element_text(size = 12)
    ) +
    guides(fill = guide_legend(nrow = 2, direction = "horizontal", title.position = "top"))

  return(pl)
}


#' Plot cumulative distribution
plot_cumulative_distribution <- function(datalist, name_dataset, P) {

  pl <- datalist$df_nr_negevents_pp_py %>%
    dplyr::group_by(t_id, nr_occur_cumsum) %>%
    dplyr::summarise(frequency = n(), .groups = 'drop') %>%
    dplyr::group_by(t_id) %>%
    dplyr::mutate(
      frequency_norm = frequency / sum(frequency),
      nr_obs = sum(frequency),
      time_label = sprintf("Observed year %02d (n = %d)", t_id, nr_obs) %>%
        stringr::str_wrap(width = 9)
    ) %>%
    ggplot() +
    geom_bar(
      aes(x = nr_occur_cumsum, y = frequency_norm),
      linewidth = 0.5,
      stat = 'identity',
      position = 'identity'
    ) +
    P$own_theme +
    ggh4x::facet_wrap2(
      time_label ~ .,
      scales = "free_y",
      ncol = 2,
      strip.position = 'left',
      dir = 'v',
      labeller = label_wrap_gen(18)
    ) +
    scale_y_continuous(
      n.breaks = 3,
      expand = expansion(add = c(0, 0)),
      position = "right"
    ) +
    scale_x_continuous(expand = expansion(mult = c(0.01, 0))) +
    theme(
      panel.spacing.y = unit(0.45, "lines"),
      title = element_text(size = 10),
      axis.text.y = element_text(size = 9),
      strip.text.y.left = element_text(
        size = 9,
        angle = 0,
        margin = margin(t = 1, r = 1, b = 1, l = 1)
      ),
      legend.position = 'right'
    ) +
    labs(
      title = sprintf("Empirical distribution of\ncumulative number of life events (%s)", name_dataset),
      y = "Frequency (normalized)",
      x = "Number of events per person"
    )

  return(pl)
}


#' Plot age distribution
plot_age_distribution <- function(datalist, name_dataset, P) {

  pl <- datalist$df_per_event_pp_py %>%
    dplyr::select(p_id, age) %>%
    dplyr::distinct() %>%
    dplyr::mutate(dataset = name_dataset) %>%
    ggplot() +
    geom_histogram(
      aes(x = age, y = after_stat(count) / sum(after_stat(count)) * 100),
      linewidth = 0.3,
      col = 'grey30',
      binwidth = 1
    ) +
    ggh4x::facet_grid2(. ~ dataset) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(x = "Age (years)", y = "Percentage") +
    P$own_theme +
    theme(
      plot.margin = unit(c(1, 1, 1, 1), "cm"),
      axis.text = element_text(size = 10),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 12),
      axis.title = element_text(size = 14),
      strip.text = element_text(size = 12)
    )

  return(pl)
}


#' Plot demographics distribution
plot_demo_distribution <- function(datalist, name_dataset, P,
                                   nr_obs_years = c(10, 15, 20)) {

  # Create dataframt with summary demographics per subsample
  df <- datalist$df_nr_negevents_pp_py %>%
    select(p_id, t_id,nr_years_obs,
           nr_years_missing_in_timespan,
           age, sex, ses, civsta, edyear, wstat)

  head(df)

  df_demo <- lapply(c(0, nr_obs_years), function(nr_years){
    temp <- df

    if (nr_years > 0){
      temp <- temp %>%
        dplyr::filter(
          nr_years_missing_in_timespan == 0,
          nr_years_obs >= !!nr_years,
          t_id <= !!nr_years
        )
    }

    temp <- temp %>%
      select(p_id, age, sex, ses, civsta, edyear, wstat) %>%
      dplyr::distinct() %>%
      tidyr::gather(variable, value, -p_id) %>%
      dplyr::filter(!is.na(variable), !is.na(value)) %>%
      group_by(variable, value) %>%
      dplyr::summarise(
        count = n(),
        .groups = 'drop'
      ) %>%
        # Normalize
        group_by(variable) %>%
        dplyr::mutate(frequency_norm = count / sum(count)) %>%
        ungroup() %>%

        dplyr::mutate(subsample = ifelse(nr_years > 0,
                                         sprintf("Observed for %d consecutive years", nr_years),
                                         "Full sample"))
    # temp <- temp %>%
    #   dplyr::mutate(value = ifelse(variable == "age",
    #                                as.numeric(value), value))
    #
    # temp <- temp %>%
    #   dplyr::mutate(value = dplyr::case_when(
    #     variable == "age" ~ as.character(as.numeric(value)),
    #     TRUE ~ value
    #   ))

    if (name_dataset == "SHP"){
      temp <- temp %>%
        dplyr::filter(variable %in% c("age", "sex", "wstat", "edyear")) %>%
        dplyr::mutate(variable = dplyr::recode(variable,
                                                 "age" = "Age",
                                                 "sex" = "Sex",
                                                 "wstat" = "Working\nStatus",
                                                 "edyear" = "Education\nYears"),
                        variable = factor(variable,
                                        levels = c("Age", "Sex",
                                                   "Working\nStatus",
                                                   "Education\nYears")))
    } else if (name_dataset == "HILDA"){
      temp <- temp %>%
        dplyr::filter(variable %in% c("age", "sex", "ses"))  %>%
        dplyr::mutate(variable = dplyr::recode(variable,
                                               "age" = "Age",
                                               "sex" = "Sex",
                                               "ses" = "Socio-Economic\nStatus"),
                      variable = factor(variable,
                                        levels = c("Age", "Sex",
                                                   "Socio-Economic\nStatus")))
    }


    return(temp)
  }) %>% do.call(rbind, .) %>% as.data.frame() %>%
    dplyr::mutate(name_dataset = name_dataset)

  range_age <- range(as.numeric(df_demo[df_demo[["variable"]] == "Age", "value"]))

  if (name_dataset == "SHP"){
    range_edyears <- range(as.numeric(df_demo[df_demo[["variable"]] == "Education\nYears", "value"]))
    scales <- list(
      # scale_x_continuous(n.breaks = 4,
                         # expand = expansion(mult = c(0, 0.075))),
      scale_x_discrete(breaks = c("20", "40", "60", "80"),
                       limits = as.character(seq(range_age[1], range_age[2]))),
      scale_x_discrete(), scale_x_discrete(),
      scale_x_discrete(breaks = c("5", "10", "15", "20"),
                       limits = as.character(seq(range_edyears[1], range_edyears[2])))
    )
  } else if (name_dataset == "HILDA"){
    range_ses <- range(as.numeric(df_demo[df_demo[["variable"]] == "Socio-Economic\nStatus", "value"]))

    scales <- list(
      # scale_x_continuous(n.breaks = 4,
      # expand = expansion(mult = c(0, 0.075))),
      scale_x_discrete(breaks = c("20", "40", "60", "80"),
                       limits = as.character(seq(range_age[1], range_age[2]))),
      scale_x_discrete(),
      scale_x_discrete(breaks = c("2", "4", "6", "8", "10"),
                                           limits = as.character(seq(range_ses[1], range_ses[2])))
    )
  }

  pl <- df_demo %>%
    ggplot() +
    geom_bar(
      aes(x = value, y = frequency_norm),
      linewidth = 0.5,
      stat = 'identity',
      position = 'identity'
    ) +
    P$own_theme +
    ggh4x::facet_grid2(
      subsample ~ variable,
      scales = "free",
      # ncol = 2,
      switch = 'y',
      independent = "y",
      # dir = 'v',
      labeller = label_wrap_gen(12)
    ) +
    scale_y_continuous(
      n.breaks = 3,
      expand = expansion(add = c(0, 0)),
      position = "right"
    ) +
    # scale_x_continuous(expand = expansion(mult = c(0.01, 0))) +
    theme(
      panel.spacing.y = unit(0.45, "lines"),
      title = element_text(size = 10),
      axis.text.x = element_text(angle = 90, size = 9),
      axis.text.y = element_text(size = 9),
      strip.text.x.top = element_text(
        size = 10),
      strip.text.y.left = element_text(
        size = 8,
        angle = 0,
        margin = margin(t = 1, r = 1, b = 1, l = 1)
      ),
      legend.position = 'right'
    ) +
    labs(
      title = sprintf("Demographics (%s)", name_dataset),
      y = "Frequency (normalized)",
      x = ""
    ) +
    ggh4x::facetted_pos_scales(x = scales)

  pl

  return(pl)
}



#' Plot observation years distribution
plot_observation_years <- function(datalist, name_dataset, P) {

  pl <- datalist$df_nr_negevents_pp %>%
    dplyr::select(p_id, nr_years_obs) %>%
    dplyr::distinct() %>%
    dplyr::mutate(dataset = name_dataset) %>%
    ggplot() +
    geom_histogram(
      aes(x = nr_years_obs, y = after_stat(count) / sum(after_stat(count)) * 100),
      linewidth = 0.3,
      col = 'grey30',
      binwidth = 1
    ) +
    ggh4x::facet_grid2(. ~ dataset) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(x = "Number of observed years", y = "Percentage") +
    P$own_theme +
    theme(
      plot.margin = unit(c(1, 1, 1, 1), "cm"),
      axis.text = element_text(size = 10),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 12),
      axis.title = element_text(size = 14),
      strip.text = element_text(size = 12)
    )

  return(pl)
}


# ============================================================================
# PLOT SAVING
# ============================================================================

#' Save all individual plots
#' @param plots List of plots
#' @param datalist Data list with file paths
#' @param name_dataset Dataset name
save_all_plots <- function(plots, datalist, name_dataset) {

  plot_specs <- list(
    pl_negevent_abs_freq = list(
      name = sprintf("%s_absolute_frequency_per_negevent.pdf", name_dataset),
      height = 160
    ),
    pl_negevent_freq = list(
      name = sprintf("%s_frequency_per_negevent.pdf", name_dataset),
      height = 160
    ),
    pl_nr_events = list(
      name = sprintf("%s_distribution_nr_events_across_years.pdf", name_dataset),
      height = 150
    ),
    pl_nr_events_py = list(
      name = sprintf("%s_distribution_nr_events_per_year.pdf", name_dataset),
      height = 200
    ),
    pl_event_distr_over_time = list(
      name = sprintf("%s_areaplot_nr_events_per_year.pdf", name_dataset),
      height = 160
    ),
    pl_cum_distr = list(
      name = sprintf("%s_emp_cumulative_nr_events_per_year.pdf", name_dataset),
      height = 220
    ),
    # pl_age = list(
    #   name = sprintf("%s_age.pdf", name_dataset),
    #   height = 160
    # ),
    pl_demo = list(
      name = sprintf("%s_demo.pdf", name_dataset),
      height = 150
    ),
    pl_nr_years_obs = list(
      name = sprintf("%s_nr_years_obs.pdf", name_dataset),
      height = 160
    )
  )

  for (plot_name in names(plot_specs)) {
    filepath <- file.path(datalist$filepath_figs_dataset, plot_specs[[plot_name]]$name)
    save_plot(plots[[plot_name]], filepath, height = plot_specs[[plot_name]]$height)
    cat(sprintf("  Saved: %s\n", plot_specs[[plot_name]]$name))
  }
}


# ============================================================================
# COMBINED PLOTS
# ============================================================================

#' Create combined event frequency plot for SHP and HILDA
#' @param results_SHP SHP results
#' @param results_HILDA HILDA results
#' @param filepath_base Base file path
plot_combined_event_frequency <- function(results_SHP, results_HILDA, filepath_base) {

  # Adjust text size in plots
  results_SHP$plots$pl_negevent_freq$layers[[3]]$aes_params$size <- 2.25
  results_HILDA$plots$pl_negevent_freq$layers[[3]]$aes_params$size <- 2.25

  pl_combined <- (
    results_SHP$plots$pl_negevent_freq +
      labs(title = "") +
      theme(
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 14),
        axis.title.x = element_text(margin = margin(t = 10)),
        plot.margin = unit(c(0, 0, 0, 0), 'cm')
      )
  ) + (
    results_HILDA$plots$pl_negevent_freq +
      labs(title = "") +
      theme(
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 14),
        axis.title.x = element_text(margin = margin(t = 10)),
        plot.margin = unit(c(0, 0, 0, 0), 'cm')
      )
  ) +
    plot_layout(widths = c(1.05, 1)) +
    plot_layout(axis_titles = "collect")

  filepath <- file.path(filepath_base, "figs", "frequency_per_negevent.pdf")
  save_plot(pl_combined, filepath, height = 115)
  cat(sprintf("  Saved combined plot: frequency_per_negevent.pdf\n"))

  return(pl_combined)
}


#' Create combined event distribution over time plot
#' @param results_SHP SHP results
#' @param results_HILDA HILDA results
#' @param filepath_base Base file path
plot_combined_event_over_time <- function(results_SHP, results_HILDA, filepath_base) {

  pl_combined <- (
    (results_SHP$plots$pl_event_distr_over_time + theme(legend.position = "none")) +
      (results_HILDA$plots$pl_event_distr_over_time +
         theme(legend.margin = margin(t = 10, r = 20, b = 10, l = 300, unit = "pt"))) +
      plot_layout(axis_titles = "collect")
  ) +
    guide_area() +
    plot_layout(guides = 'collect') +
    plot_layout(heights = c(1, 0.3))

  filepath <- file.path(filepath_base, "figs", "areaplot_nr_events_per_year.pdf")
  save_plot(pl_combined, filepath, height = 120)
  cat(sprintf("  Saved combined plot: areaplot_nr_events_per_year.pdf\n"))

  return(pl_combined)
}


#' #' Create combined age distribution plot
#' #' @param results_SHP SHP results
#' #' @param results_HILDA HILDA results
#' #' @param filepath_base Base file path
#' plot_combined_age <- function(results_SHP, results_HILDA, filepath_base) {
#'
#'   pl_combined <- (
#'     results_SHP$plots$pl_age +
#'       results_HILDA$plots$pl_age +
#'       plot_layout(axis_titles = "collect")
#'   ) +
#'     guide_area() +
#'     plot_layout(heights = c(1, 0.3))
#'
#'   filepath <- file.path(filepath_base, "figs", "age.pdf")
#'   save_plot(pl_combined, filepath, height = 100)
#'   cat(sprintf("  Saved combined plot: age.pdf\n"))
#'
#'   return(pl_combined)
#' }


#' Create combined observation years plot
#' @param results_SHP SHP results
#' @param results_HILDA HILDA results
#' @param filepath_base Base file path
plot_combined_observation_years <- function(results_SHP, results_HILDA, filepath_base) {

  pl_combined <- (
    results_SHP$plots$pl_nr_years_obs +
      results_HILDA$plots$pl_nr_years_obs +
      plot_layout(axis_titles = "collect")
  ) +
    guide_area() +
    plot_layout(heights = c(1, 0.3))

  filepath <- file.path(filepath_base, "figs", "nr_years_obs.pdf")
  save_plot(pl_combined, filepath, height = 100)
  cat(sprintf("  Saved combined plot: nr_years_obs.pdf\n"))

  return(pl_combined)
}


# ============================================================================
# COMPARISON TABLE
# ============================================================================

#' Create comparison table of descriptive statistics
#' @param num_desc_SHP Numerical descriptives for SHP
#' @param num_desc_HILDA Numerical descriptives for HILDA
create_descriptive_comparison_table <- function(num_desc_SHP, num_desc_HILDA) {

  df_comparison <- rbind(num_desc_SHP, num_desc_HILDA)

  # Format table
  table <- df_comparison %>%
    kableExtra::kbl(
      format = "latex",
      escape = FALSE,
      caption = "Descriptive statistics comparison between SHP and HILDA datasets",
      booktabs = TRUE,
      label = "descriptive_comparison"
    ) %>%
    kableExtra::kable_styling(latex_options = c("striped", "scale_down"))

  return(table)
}
