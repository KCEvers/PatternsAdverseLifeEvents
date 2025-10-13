##### DESCRIPTIVE ANALYSIS OF LIFE EVENTS #####

### Set-up file with libraries, functions, and parameters
filepath_base = dirname(dirname(rstudioapi::getActiveDocumentContext()$path))
stopifnot(nzchar(filepath_base))

source(file.path(filepath_base, "code/setup.R"))
source(file.path(filepath_base, "code/prepare_data_func.R"))

# Get parameters
P = get_parameters()

# Function to plot descriptives
plot_descriptives = function(name_dataset, filepath_base, P, rerun = FALSE){

  # Load data
  datalist = prepare_data(name_dataset, filepath_base, P$event_dicts[[name_dataset]], rerun = rerun)

  # Set height and width of figure
  if (name_dataset == "HILDA"){
    height = 160
    max_char = 25
  } else if ("SHP" %in% name_dataset){
    height = 160
    max_char = 25
  }
  text_size = 10

  # Frequency dataframe
  df_freq = datalist$df_per_event %>%
    filter(valence == "negative") %>%
    dplyr::mutate(dataset = !!name_dataset,
                  perc_occur = nr_occur / (nr_occur + nr_nooccur) * 100,
                  perc_occur_incl_missing = nr_occur / (nr_occur + nr_nooccur + nr_missing_occur) * 100,
                  sum_n = nr_occur + nr_nooccur + nr_missing_occur) #%>%
    # recode_events(., "event", P$event_dicts[[name_dataset]], df_per_event = datalist$df_per_event)

  # Order events by frequency
  df_freq = df_freq %>%
    dplyr::mutate_at("event",
                     ~ dplyr::recode_factor(.x, !!!tibble::deframe(P$event_dicts[[name_dataset]] %>%
                                                                     dplyr::select(recode_var, description) %>%
                                                                     # Order events in event_dict according to frequency
                                                                     dplyr::slice(match(df_freq %>% dplyr::arrange(perc_occur_incl_missing) %>% dplyr::pull(event) %>% as.character(), recode_var))), .ordered = T))

  # Absolute event frequency
  pl_negevent_abs_freq = df_freq %>%
    ggplot() +
    geom_bar(aes(x = nr_occur, y = event), col = P$col_data, fill = P$col_data, stat = 'identity', alpha = .5) +
    geom_text(aes(x = nr_occur + ifelse(name_dataset == "SHP", 8000, 4000), y = event, label = nr_occur), col = 'grey30', family = P$font_family,
              size = 4) +
    ggh4x::facet_grid2(. ~ dataset) +
    P$own_theme +
    theme(
      axis.text = element_text(size = text_size+2),
      plot.margin = unit(c(.01,.9,.01,.01), 'cm'),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    scale_y_discrete(labels = function(x) sapply(x, wrap_equally, max_chars = 31),
                     expand = c(0,0)) +
    scale_x_continuous(n.breaks = 4, expand = expansion(mult = c(0, .075))) +
    labs(x = "Frequency in person-years", y = "", title = "Event frequency\nacross people and years")
  pl_negevent_abs_freq

  save_plot(pl_negevent_abs_freq, file.path(datalist$filepath_figs_dataset, sprintf("%s_absolute_frequency_per_negevent.pdf", name_dataset)), height = height)

  # Relative event frequency
  # Find at what percentages 10000 instances are
  x = seq(10000, ifelse(name_dataset == "SHP", 50000, 40000), by = 10000)
  perc_10000_instances = data.frame(x = x, perc = x / median(df_freq$sum_n) * 100)

  pl_negevent_freq = df_freq %>%
    ggplot() +
    geom_vline(data = perc_10000_instances,
               aes(xintercept = perc), linetype = 'dashed', color = 'grey80') +
    geom_bar(aes(x = perc_occur_incl_missing, y = event), col = P$col_data, fill = P$col_data, stat = 'identity', alpha = .6) +
    geom_text(aes(x = perc_occur_incl_missing + ifelse(name_dataset == "SHP", .03 * 100, .02 * 100), y = event, label = sprintf("%.2f%%", perc_occur_incl_missing)), col = 'grey30', family = P$font_family,
              # fontface = 'bold',
              size = 4) +
    ggh4x::facet_grid2(. ~ dataset) +
    P$own_theme +
    theme(
      axis.text = element_text(size = text_size+2),
      title = element_text(size = text_size+6),
      plot.margin = unit(c(.01,.5,.01,.01), 'cm'),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
    ) +
    scale_y_discrete(labels = function(x) sapply(x, wrap_equally, max_chars = max_char),
                     expand = c(0,0)) +
    scale_x_continuous(n.breaks = 4, expand = expansion(mult = c(0, .1))) +
    labs(x = "Frequency (%) of events across all person-years", y = "", title = "Event frequency\nacross people and years")
  pl_negevent_freq

  save_plot(pl_negevent_freq, file.path(datalist$filepath_figs_dataset, sprintf("%s_frequency_per_negevent.pdf", name_dataset)), height = height)


  # Number of events per year

  # Plot number of events across years
  pl_nr_events = datalist$df_nr_negevents_pp_py %>%
    filter(!is.na(t_id)) %>%
    group_by(nr_occur) %>%
    dplyr::summarise(n = n(), .groups = 'drop') %>%
    mutate(n_norm = n / sum(n)) %>%
    ggplot() +
    geom_bar(aes(x = nr_occur, y = n), stat = 'identity')+
    P$own_theme +
    theme(
      panel.grid.minor = element_blank(),
    ) +
    scale_y_continuous(n.breaks = 5, expand = expansion(mult = c(0, .01))) +
    scale_x_continuous(n.breaks = max(datalist$df_nr_negevents_pp_py$nr_occur), expand = c(0,0)) +
    labs(title = "Empirical distribution of number of adverse life events per year",
         subtitle = name_dataset, x = "Number of events per year", y = "Frequency in person-years") +
    labs(title ="", subtitle = "")
  pl_nr_events
  save_plot(pl_nr_events, file.path(datalist$filepath_figs_dataset, sprintf("%s_distribution_nr_events_across_years.pdf", name_dataset)), height = 150)


  # Plot number of events per observation year
  pl_nr_events_py = datalist$df_nr_negevents_pp_py %>%
    filter(!is.na(t_id)) %>%
    group_by(t_id, nr_occur) %>%
    dplyr::summarise(n = n(), .groups = 'drop') %>%
    group_by(t_id) %>%
    mutate(n_norm = n / sum(n)) %>%
    mutate(t_id_label = sprintf("Year %02d, n = %d", t_id, sum(n)) %>% stringr::str_wrap(width = 9)) %>%
    ungroup() %>%
    ggplot() +
    geom_bar(aes(x = nr_occur, y = n_norm), stat = 'identity')+
    P$own_theme +
    theme(axis.text.y = element_text(size=7),
          plot.title = element_text(hjust = .5),
          plot.margin = unit(c(.5,.5,.5,.5), 'cm'),
          text = element_text(family = "serif"),
          title = element_text(family = "serif", size = 12),
          plot.subtitle = element_text(hjust = .5)) +
    scale_x_continuous(n.breaks = 4, expand = c(0,0)) +
    labs(title = "Empirical distribution of\nnumber of adverse life events per year",
         subtitle = name_dataset, x = "Number of events per year", y = "Frequency (normalized per year)") +
    ggh4x::facet_wrap2(t_id_label ~ ., ncol = 2, dir = 'v', strip.position = 'right') +
    # Spacing between facets
    theme(panel.spacing.y = unit(0.4, "lines"),
          strip.background=element_rect(colour=P$col_facets,
                                        fill=P$fill_facets),
          strip.text.y = element_text(angle = 0))
  pl_nr_events_py
  save_plot(pl_nr_events_py, file.path(datalist$filepath_figs_dataset, sprintf("%s_distribution_nr_events_per_year.pdf", name_dataset)), height = 200)

  # Distribution of number of events across years in area plot
  pl_event_distr_over_time = datalist$df_nr_negevents_pp_py %>%
    filter(!is.na(t_id)) %>%
    group_by(t_id, nr_occur) %>%
    dplyr::summarise(n = n(), .groups = 'drop') %>%
    group_by(t_id) %>%
    mutate(n_norm = n / sum(n)) %>% ungroup() %>%
    dplyr::mutate(nr_occur = factor(nr_occur %>% as.character(), levels = as.character(0:13)),
                  dataset = !!name_dataset) %>%
    ggplot() +
    geom_area(aes(x = t_id, y = n_norm, fill = nr_occur), col = 'grey30',
    linewidth = .8,
    position = 'stack'
    ) +
    ggh4x::facet_grid2(.~dataset) +
    viridis::scale_fill_viridis(discrete = T, direction = -1, name = "Number of events",
                                drop=F) +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = expansion(mult = c(0.002, .002))
    ) +
    labs(x = "Observation year", y = "Percentage (normalized per year)") +
    P$own_theme +
    theme(legend.position = 'bottom') +
    theme(
      plot.margin=unit(c(1, 1, 1, 1),"cm"),
      axis.text = element_text(size = 10),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 12),
      axis.title = element_text(size = 14),
      strip.text = element_text(size = 12)
    ) + guides(fill = guide_legend(nrow =2, direction = "horizontal",
                                   title.position = "top"))
  pl_event_distr_over_time

  save_plot(pl_event_distr_over_time,
            file.path(datalist$filepath_figs_dataset, sprintf("%s_areaplot_nr_events_per_year.pdf", name_dataset)), height = 160)

  # Cumulative number of events
  pl_cum_distr = datalist$df_nr_negevents_pp_py %>%
    group_by(t_id, nr_occur_cumsum) %>%
    dplyr::summarise(frequency = n(), .groups = 'drop') %>%
    group_by(t_id) %>%
    dplyr::mutate(frequency_norm = frequency / sum(frequency),
                  nr_obs = sum(frequency)) %>%
    mutate(time_label = sprintf("Observed year %02d (n = %d)", t_id, nr_obs) %>% stringr::str_wrap(width = 9)) %>%
    ggplot() +
    geom_bar(aes(x = nr_occur_cumsum, y = frequency_norm), linewidth = .5, stat = 'identity', position = 'identity') +
    P$own_theme +
    ggh4x::facet_wrap2(time_label ~ ., scales = "free_y",
                       # switch = "y",
                       ncol = 2,
                       strip.position = 'left',
                       dir = 'v',
                       labeller = label_wrap_gen(18)
    ) +
    scale_y_continuous( n.breaks = 3, expand = expansion(add = c(0, 0)), position = "right") +
    scale_x_continuous(expand = expansion(mult = c(0.01, 0))) +
    # Spacing between facets
    theme(panel.spacing.y = unit(0.45, "lines")) +
    theme(
      title = element_text(size = 10),
      axis.text.y = element_text(size=9),
      strip.text.y.left = element_text(size = 9, angle = 0, margin = margin(t = 1, r = 1, b = 1, l = 1)),
      legend.position = 'right') +
    labs(title = sprintf("Empirical distribution of\ncumulative number of life events (%s)", name_dataset),
         y = "Frequency (normalized)", x = "Number of events per person")
  pl_cum_distr

  save_plot(pl_cum_distr, file.path(datalist$filepath_figs_dataset, sprintf("%s_emp_cumulative_nr_events_per_year.pdf", name_dataset)), height = 220)

  # Plot age distribution across person-years
  pl_age = datalist$df_per_event_pp_py %>%
    dplyr::select(p_id, age) %>% distinct() %>%
    dplyr::mutate(dataset = !!name_dataset) %>%
    ggplot() +
    geom_histogram(aes(x = age, y = after_stat(count) / sum(after_stat(count)) * 100),
                   linewidth = .3,
                   col = 'grey30',
                   binwidth = 1) +
    ggh4x::facet_grid2(.~dataset) +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = expansion(mult = c(0, .05))
    ) +
    labs(x = "Age (years)", y = "Percentage") +
    P$own_theme +
    theme(
      plot.margin=unit(c(1, 1, 1, 1),"cm"),
      axis.text = element_text(size = 10),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 12),
      axis.title = element_text(size = 14),
      strip.text = element_text(size = 12)
    )
  pl_age

  save_plot(pl_age,
            file.path(datalist$filepath_figs_dataset, sprintf("%s_age.pdf", name_dataset)), height = 160)

  # Number of available years per participant
  temp = datalist$df_nr_negevents_pp %>%
    dplyr::select(p_id, nr_years_obs) %>% distinct()
  nrow(temp)

  pl_nr_years_obs = datalist$df_nr_negevents_pp %>%
    dplyr::select(p_id, nr_years_obs) %>% distinct() %>%
    dplyr::mutate(dataset = !!name_dataset) %>%
    ggplot() +
    geom_histogram(aes(x = nr_years_obs, y = after_stat(count) / sum(after_stat(count)) * 100),
                   linewidth = .3,
                   col = 'grey30',
                   binwidth = 1) +
    ggh4x::facet_grid2(.~dataset) +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = expansion(mult = c(0, .05))
    ) +
    labs(x = "Number of observed years", y = "Percentage") +
    P$own_theme +
    theme(
      plot.margin=unit(c(1, 1, 1, 1),"cm"),
      axis.text = element_text(size = 10),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 12),
      axis.title = element_text(size = 14),
      strip.text = element_text(size = 12)
    )
  pl_nr_years_obs
  save_plot(pl_nr_years_obs,
            file.path(datalist$filepath_figs_dataset, sprintf("%s_nr_years_obs.pdf", name_dataset)),
            height = 160)


  return(list(
    datalist = datalist,
    df_freq = df_freq,
    pl_negevent_abs_freq = pl_negevent_abs_freq,
    pl_negevent_freq = pl_negevent_freq,
    pl_nr_events = pl_nr_events,
    pl_nr_events_py = pl_nr_events_py,
    pl_event_distr_over_time = pl_event_distr_over_time,
    pl_cum_distr = pl_cum_distr,
    pl_age = pl_age,
    pl_nr_years_obs = pl_nr_years_obs
  ))
}

# Function to get numerical descriptives
num_descriptives = function(name_dataset, filepath_base, P, rerun = FALSE){

  # Load data
  datalist = prepare_data(name_dataset, filepath_base, P$event_dicts[[name_dataset]], rerun = rerun)

  # Set up storage
  res <- data.frame(dataset = name_dataset)

  # Top two most frequent events
  df_ = datalist$df_per_event %>% filter(valence == "negative") %>% arrange(desc(nr_occur)) %>%
    mutate(event = as.character(event)) %>%
    mutate(nr_occur_norm = nr_occur / sum(nr_occur))

  res$top1_event = df_$event[1]
  res$top2_event = df_$event[2]
  res$total_nr_neg_events = sum(df_$nr_occur)
  res$perc_top_two_neg_events = sprintf("The top two events were responsible for %.2f%% of all events", sum(df_ %>% slice(1:2) %>% pull(nr_occur)) / res$total_nr_neg_events * 100)

  # Number of person-years < 18
  res$nr_personyears_underage = datalist$df_binary %>%
    # Only adults
    dplyr::filter(age < 18) %>%
    dplyr::select(p_id, wave_nr) %>%
    distinct() %>%
    nrow()

  # Number of person-years missing all adverse events
  res$nr_all_negevents_missing = df_binary %>%
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
    mutate(valence = dplyr::recode(event, !!!tibble::deframe(event_dict %>% dplyr::select(recode_var, valence)))  ) %>%
    mutate(dependence = dplyr::recode(event, !!!tibble::deframe(event_dict %>% dplyr::select(recode_var, in_dependent)))  ) %>%
    arrange(p_id, wave_nr) %>%
    group_by(p_id, wave_nr) %>%
    # Only keep person-years with all negative events missing
    dplyr::filter((sum(nr_occur[valence == "negative"]) == 0 & sum(nr_nooccur[valence == "negative"]) == 0)) %>%
    dplyr::ungroup() %>%
    dplyr::select(p_id, wave_nr) %>%
    distinct() %>%
    nrow()

  # Sex of sample
  temp = datalist$df_per_event_pp_py %>%
    dplyr::select(p_id, sex) %>%
    group_by(p_id) %>% slice(1) %>% pull(sex) %>%
    table()
  round(temp / sum(temp) * 100, 3)

  # Number of person-years and households overall sample
  temp = datalist$df_per_event_pp_py %>%
    dplyr::select(p_id, crosswave_h_id, age) %>% distinct()
  res$nr_p_id_overall = nrow(temp)
  res$nr_h_id_overall = length(unique(temp$crosswave_h_id))


  # Co-occurrence
  df_ = datalist$df_nr_negevents_pp_py %>%
    filter(!is.na(t_id)) %>%
    group_by(nr_occur) %>%
    dplyr::summarise(n = n(), .groups = 'drop') %>%
    mutate(n_norm = n / sum(n))
  res$perc_no_event = df_ %>% filter(nr_occur == 0) %>% pull(n_norm)
  res$perc_one_event = df_ %>% filter(nr_occur == 1) %>% pull(n_norm)
  res$perc_cooccur = df_ %>% filter(nr_occur > 1) %>% pull(n_norm) %>% sum()

  # Viewed another way: when an event occurs, how often does it co-occur?
  df_ = df_ %>% filter(nr_occur > 0) %>% mutate(n_norm = n / sum(n))
  res$perc_one_event_if_event = df_ %>% filter(nr_occur == 1) %>% pull(n_norm)
  res$perc_cooccur_if_event = df_ %>% filter(nr_occur > 1) %>% pull(n_norm) %>% sum()

  return(res)

}



# Plot descriptives of SHP and HILDA
plots_SHP = plot_descriptives("SHP", filepath_base, P)
plots_HILDA = plot_descriptives("HILDA", filepath_base, P)

# Numerical descriptives of SHP and HILDA
num_SHP = num_descriptives("SHP", filepath_base, P)
num_HILDA = num_descriptives("HILDA", filepath_base, P)


# Create combined plots of SHP and HILDA: Frequency of events

# Change size of geom_text
plots_SHP$pl_negevent_freq$layers[[3]]$aes_params$size = 2.25
plots_HILDA$pl_negevent_freq$layers[[3]]$aes_params$size = 2.25

pl_negevent_freq = (plots_SHP$pl_negevent_freq + labs(title = "") +
                      theme(axis.text = element_text(size = 7),
                            axis.title = element_text(size = 14),
                            axis.title.x = element_text(margin = margin(t = 10)),
                            plot.margin = unit(c(0,0,0,0), 'cm'))) + (plots_HILDA$pl_negevent_freq + labs(title = "") +
                                                                        theme(axis.text = element_text(size = 7),
                                                                              axis.title = element_text(size = 14),
                                                                              axis.title.x = element_text(margin = margin(t = 10)),
                                                                              plot.margin = unit(c(0,0,0,0), 'cm'))) +   plot_layout(widths = c(1.05, 1))  +
  plot_layout(axis_titles = "collect")
pl_negevent_freq

save_plot(pl_negevent_freq, file.path(plots_SHP$datalist$filepath_figs,
                                      sprintf("frequency_per_negevent.pdf")), height = 115)


# Create combined plots of SHP and HILDA: Consistency distribution over time
pl_event_distr_over_time = ((plots_SHP$pl_event_distr_over_time + theme(legend.position = "none")) +
                              (plots_HILDA$pl_event_distr_over_time + theme(legend.margin = margin(t = 10, r = 20, b = 10, l = 300, unit = "pt"))) +
                              plot_layout(axis_titles = "collect")) + guide_area() +
  plot_layout(guides = 'collect') + plot_layout(heights = c(1, .3))
pl_event_distr_over_time

save_plot(pl_event_distr_over_time, file.path(plots_SHP$datalist$filepath_figs,
                                              "areaplot_nr_events_per_year.pdf"), height = 120)



# Create combined plots of SHP and HILDA: Age
pl_age = ((plots_SHP$pl_age) +
                              (plots_HILDA$pl_age) +
                              plot_layout(axis_titles = "collect")) +
  guide_area() +
  # plot_layout(guides = 'collect') +
  plot_layout(heights = c(1, .3))
pl_age

save_plot(pl_age, file.path(plots_SHP$datalist$filepath_figs,
                                              "age.pdf"), height = 100)


# Number of observed years
pl_nr_years_obs = ((plots_SHP$pl_nr_years_obs) +
            (plots_HILDA$pl_nr_years_obs) +
            plot_layout(axis_titles = "collect")) +
  guide_area() +
  # plot_layout(guides = 'collect') +
  plot_layout(heights = c(1, .3))
pl_nr_years_obs

save_plot(pl_nr_years_obs, file.path(plots_SHP$datalist$filepath_figs,
                            "nr_years_obs.pdf"), height = 100)
