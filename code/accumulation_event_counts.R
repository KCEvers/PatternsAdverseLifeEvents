##### ACCUMULATION OF ADVERSE LIFE EVENT COUNTS - REFACTORED #####

#' Main wrapper function for accumulation analysis
#' @param name_dataset Dataset name (e.g., "SHP", "HILDA")
#' @param filepath_base Base file path
#' @param P Parameters list
#' @param rerun Whether to rerun analysis
run_accumulation <- function(name_dataset, filepath_base, P,
                             nr_years = 20, rerun = FALSE) {

  # 1. Prepare data
  cat("\n=== Preparing data ===\n")
  event_dict <- P$event_dicts[[name_dataset]]
  datalist <- prepare_data(name_dataset, filepath_base, event_dict)

  filepath <- file.path(
    datalist$filepath_deriv,
    sprintf("accemulation_%dyears.RDS", nr_years)
  )

  if (!file.exists(filepath) || rerun){
    # 2. Calculate number of trials for Polya urn
    cat("\n=== Calculating trial parameters ===\n")
    n_trials <- calculate_n_trials(name_dataset, event_dict)

    # 3. Create accumulation dataset
    cat("\n=== Creating accumulation dataset ===\n")
    df_accum <- create_accumulation_dataset(datalist, n_trials, nr_years = nr_years)

    # Assign to global environment (required for insight package functions)
    assign("df_accum", df_accum, envir = .GlobalEnv)

    # 4. Fit models
    cat("\n=== Fitting models ===\n")
    models <- fit_accumulation_models(df_accum)

    # 5. Extract model parameters
    cat("\n=== Extracting parameters ===\n")
    polya_params <- extract_polya_parameters(models$mod_polya)
    sample_sizes <- extract_sample_sizes(models$mod_frailty, nr_years)

    # 6. Compare model fit
    cat("\n=== Comparing model fit ===\n")
    fit_comparison <- compare_model_fit(models)

    # 7. Summarize models
    cat("\n=== Summarizing models ===\n")
    model_summaries <- list(
        mod_summ_poisson = summarise_model(models$mod_poisson),
        mod_summ_frailty = summarise_model(models$mod_frailty),
        mod_summ_polya = summarise_model(models$mod_polya)
      )

    # 8. Return consolidated results
    results <- list(
      datalist = datalist,
      df_accum = df_accum,
      nr_people = sample_sizes$nr_people,
      nr_households = sample_sizes$nr_households,
      mod_poisson = models$mod_poisson,
      mod_frailty = models$mod_frailty,
      mod_polya = models$mod_polya,
      mu = polya_params$mu,
      phi = polya_params$phi,
      alpha = polya_params$alpha,
      beta = polya_params$beta,
      mod_summ_poisson = model_summaries$mod_summ_poisson,
      mod_summ_frailty = model_summaries$mod_summ_frailty,
      mod_summ_polya = model_summaries$mod_summ_polya,
      df_ABIC = fit_comparison$df_ABIC,
      df_comp_AIC = fit_comparison$df_comp_AIC,
      df_comp_BIC = fit_comparison$df_comp_BIC
    )

    saveRDS(results, filepath)
    rm(df_accum, envir = .GlobalEnv)
  }
  results <- readRDS(filepath)

  return(results)
}


# ============================================================================
# DATA PREPARATION
# ============================================================================

#' Calculate number of trials for Polya urn model
calculate_n_trials <- function(name_dataset, event_dict) {

  if (name_dataset == "SHP") {
    # Remove personal accidents/illnesses (only one can occur) and "other" events
    n_trials <- event_dict %>%
      dplyr::filter(
        valence == "negative",
        !grepl("PL01", var_name),
        !grepl("PL36", var_name)
      ) %>%
      nrow() + 1  # Add back 1 for personal accident or illness

  } else if (name_dataset == "HILDA") {
    n_trials <- event_dict %>%
      dplyr::filter(valence == "negative") %>%
      nrow()
  }

  cat(sprintf("  Number of trials: %d\n", n_trials))
  return(n_trials)
}


#' Create accumulation dataset with complete observations
create_accumulation_dataset <- function(datalist, n_trials, nr_years = 20) {

  df_accum <- datalist$df_nr_negevents_pp_py %>%
    dplyr::filter(
      nr_years_missing_in_timespan == 0,
      nr_years_obs >= !!nr_years,
      t_id <= !!nr_years
    ) %>%
    dplyr::select(p_id, crosswave_h_id, t_id, nr_occur, nr_nooccur, nr_missing_occur) %>%
    dplyr::rename(nr_nooccur_old = nr_nooccur) %>%
    dplyr::arrange(p_id, t_id) %>%
    # Recalculate no-occurrences excluding missing data
    dplyr::mutate(nr_nooccur = n_trials - nr_occur)

  # Report missing events
  cat("  Missing events distribution:\n")
  print(table(df_accum$nr_missing_occur))

  return(df_accum)
}


# ============================================================================
# MODEL FITTING
# ============================================================================

#' Fit all three accumulation models
fit_accumulation_models <- function(df_accum) {

  # 1. Poisson model
  cat("  Fitting Poisson model...\n")
  mod_poisson <- glmmTMB::glmmTMB(
    nr_occur ~ 1,
    data = df_accum,
    family = poisson(link = "log")
  )

  # 2. Frailty model (Poisson with random effects)
  cat("  Fitting Frailty model...\n")
  mod_frailty <- glmmTMB::glmmTMB(
    nr_occur ~ 1 + (1 | p_id) + (1 | crosswave_h_id),
    data = df_accum,
    family = poisson(link = "log")
  )

  # 3. Polya urn model (beta-binomial)
  cat("  Fitting Polya urn model...\n")
  mod_polya <- glmmTMB::glmmTMB(
    cbind(nr_occur, nr_nooccur) ~ 1 + (1 | p_id) + (1 | crosswave_h_id),
    data = df_accum,
    family = betabinomial(link = "logit")
  )

  return(list(
    mod_poisson = mod_poisson,
    mod_frailty = mod_frailty,
    mod_polya = mod_polya
  ))
}


# ============================================================================
# PARAMETER EXTRACTION
# ============================================================================

#' Extract Polya urn parameters (alpha and beta)
extract_polya_parameters <- function(mod_polya) {

  mu <- plogis(mod_polya$fit$par[["beta"]])  # mu = a / (a + b) = prob
  phi <- exp(mod_polya$fit$par[["betadisp"]])  # phi = a + b
  alpha <- mu * phi
  beta <- phi - alpha

  return(list(
    mu = mu,
    phi = phi,
    alpha = alpha,
    beta = beta
  ))
}


#' Extract sample sizes from frailty model
extract_sample_sizes <- function(model, nr_years) {

  nr_people <- nrow(ranef(model)$cond$p_id)
  nr_households <- nrow(ranef(model)$cond$crosswave_h_id)

  cat(sprintf("  Number of people: %d\n", nr_people))
  cat(sprintf("  Number of households: %d\n", nr_households))
  cat(sprintf("  Number of observations: %d\n", nr_people * nr_years))

  return(list(
    nr_people = nr_people,
    nr_households = nr_households
  ))
}


# ============================================================================
# MODEL COMPARISON
# ============================================================================

#' Compare model fit using AIC and BIC
compare_model_fit <- function(models) {

  # Extract AIC and BIC
  df_ABIC <- cbind(
    AIC(models$mod_poisson, models$mod_frailty, models$mod_polya),
    BIC(models$mod_poisson, models$mod_frailty, models$mod_polya)
  )

  # Create pairwise comparisons
  model_names <- c("poisson", "frailty", "polya")

  df_comp_AIC <- outer(
    df_ABIC$AIC %>% stats::setNames(model_names),
    df_ABIC$AIC %>% stats::setNames(model_names),
    FUN = "-"
  )

  df_comp_BIC <- outer(
    df_ABIC$BIC %>% stats::setNames(model_names),
    df_ABIC$BIC %>% stats::setNames(model_names),
    FUN = "-"
  )

  return(list(
    df_ABIC = df_ABIC,
    df_comp_AIC = df_comp_AIC,
    df_comp_BIC = df_comp_BIC
  ))
}



# ============================================================================
# FORMATTING AND TABLES
# ============================================================================

#' Format model estimates for display
#' @param df_est Dataframe with model estimates
#' @param effect Type of effect ("fixed" or "ran_pars")
#' @param group Group name for random effects (e.g., "p_id")
#' @param transform Transformation to apply ("", "exp", or "prob")
format_accumulation_estimates <- function(df_est, effect, group = NULL,
                                          transform = "") {

  df <- df_est %>%
    dplyr::filter(effect == !!effect) %>%
    dplyr::rename(lower = "2.5 %", upper = "97.5 %") %>%
    dplyr::mutate(term = stringr::str_replace_all(term,
                                                  "\\(Intercept\\)", "intercept"))

  # Filter by group for random effects
  if (effect == "ran_pars") {
    df <- df %>% dplyr::filter(group == !!group)
  }

  # Apply transformation
  if (transform == "exp") {
    df <- df %>%
      dplyr::mutate_at(c("estimate", "lower", "upper"), exp)
  } else if (transform == "prob") {
    df <- df %>%
      dplyr::mutate_at(c("estimate", "lower", "upper"), ~ exp(.x) / (1 + exp(.)))
  }

  # Format as string
  df %>%
    dplyr::mutate(label = sprintf("%.2f [%.2f, %.2f]", estimate, lower, upper)) %>%
    dplyr::pull(label) %>%
    return()
}


#' Create dataframe with accumulation estimates for one dataset
create_accumulation_df <- function(name_dataset, dataset) {

  data.frame(
    dataset = name_dataset,
    poisson_intercept = format_accumulation_estimates(
      dataset$mod_summ_poisson$df_est, "fixed", transform = "exp"
    ),
    frailty_intercept = format_accumulation_estimates(
      dataset$mod_summ_frailty$df_est, "fixed", transform = "exp"
    ),
    frailty_sd_person = format_accumulation_estimates(
      dataset$mod_summ_frailty$df_est, "ran_pars", group = "p_id"
    ),
    frailty_sd_household = format_accumulation_estimates(
      dataset$mod_summ_frailty$df_est, "ran_pars", group = "crosswave_h_id"
    ),
    delta_AIC_poisson_frailty = round(
      AIC(dataset$mod_poisson) - AIC(dataset$mod_frailty), 2
    ),
    delta_BIC_poisson_frailty = round(
      BIC(dataset$mod_poisson) - BIC(dataset$mod_frailty), 2
    ),
    polya_prob = format_accumulation_estimates(
      dataset$mod_summ_polya$df_est, "fixed", transform = "prob"
    ),
    polya_dispersion = round(exp(dataset$mod_polya$fit$par[["betadisp"]]), 2),
    polya_sd_person = format_accumulation_estimates(
      dataset$mod_summ_polya$df_est, "ran_pars", group = "crosswave_h_id"
    ),
    polya_sd_household = format_accumulation_estimates(
      dataset$mod_summ_polya$df_est, "ran_pars", group = "crosswave_h_id"
    ),
    delta_AIC_frailty_polya = round(
      AIC(dataset$mod_frailty) - AIC(dataset$mod_polya), 2
    ),
    delta_BIC_frailty_polya = round(
      BIC(dataset$mod_frailty) - BIC(dataset$mod_polya), 2
    )
  )
}


#' Create LaTeX table comparing accumulation models
create_accumulation_table <- function(SHP, HILDA, nr_years) {

  rbind(
    create_accumulation_df("SHP", SHP),
    create_accumulation_df("HILDA", HILDA)
  ) %>%
    kableExtra::kbl(
      format = "latex",
      escape = FALSE,
      caption = paste0("Estimates of the Poisson, frailty, and Polya urn models in the accumulation analysis across ", nr_years, " years (profile 95\\% profile confidence intervals in brackets). Estimates are on the natural scale. Fit comparisons subtract the fit of the model on the right from that of the model on the left, such that a positive number indicates the model on the left fits better. $\\lambda$ = Yearly rate of adverse life events; $p$ = Per-trial probability of adverse life events; $\\phi$ = Dispersion; $i$ = Individual; $h$ = Household; AIC = Akaike Information Criterion; BIC = Bayesian Information Criterion (Source: Swiss Household Panel, SHP and Household, Income and Labour Dynamics in Australia Survey, HILDA)."),
      label = paste0("accumulation_", nr_years, "years"),
      booktabs = TRUE,
      col.names = NULL
    ) %>%
    kableExtra::add_header_above(
      c(" ", "$\\\\lambda$" = 1, "$\\\\lambda$" = 1, "$i$" = 1, "$h$" = 1,
        "$\\\\Delta$AIC" = 1, "$\\\\Delta$BIC" = 1,
        "$p$" = 1, "$\\\\phi$" = 1, "$i$" = 1, "$h$" = 1,
        "$\\\\Delta$AIC" = 1, "$\\\\Delta$BIC" = 1),
      escape = FALSE
    ) %>%
    kableExtra::add_header_above(
      c(" ", "Fixed" = 1, "Fixed" = 1, "Random $\\\\sigma$" = 2,
        "Fit (Poisson-Frailty)" = 2, "Fixed" = 2, "Random $\\\\sigma$" = 2,
        "Fit (Frailty-Polya)" = 2),
      escape = FALSE
    ) %>%
    kableExtra::add_header_above(
      c(" ", "Poisson" = 1, "Frailty" = 5, "Polya urn" = 6)
    ) %>%
    kableExtra::column_spec(1, "0.4in") %>%
    kableExtra::column_spec(2:13, "0.27in") %>%
    kableExtra::kable_styling(latex_options = "scale_down")
}


# ============================================================================
# SIMULATION
# ============================================================================

#' Format simulation results
#' @param sim_df Matrix of simulated data
#' @param df_accum Original accumulation dataframe
#' @param sim_type Label for simulation type
format_simulation_results <- function(sim_df, df_accum, sim_type) {

  # Bind participant and time information to each simulation
  sim_df_merged <- purrr::map(1:ncol(sim_df), function(i) {
    cbind(
      df_accum[, c("p_id", "t_id")],
      nr_occur = sim_df[, i],
      i = i
    )
  }) %>%
    do.call(rbind, .) %>%
    as.data.frame()

  # Calculate cumulative sums and frequencies
  sim_df_merged %>%
    dplyr::arrange(i, p_id, t_id) %>%
    dplyr::group_by(i, p_id) %>%
    dplyr::mutate(cumsum_nr_occur = cumsum(nr_occur)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(t_id, cumsum_nr_occur) %>%
    dplyr::summarise(frequency = n(), .groups = 'drop') %>%
    dplyr::group_by(t_id) %>%
    dplyr::mutate(frequency_norm = frequency / sum(frequency)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(sim_type = !!sim_type) %>%
    return()
}


#' Simulate data from all three models
#' @param models List containing fitted models
#' @param n_sim Number of simulations
#' @param sim_types Vector of simulation type labels
simulate_accumulation_data <- function(models, n_sim, sim_types) {

  cat("\n=== Simulating data ===\n")
  start_t <- Sys.time()

  # Simulate from Poisson model
  cat("  Simulating from Poisson model...\n")
  sim_poisson <- simulate(models$mod_poisson, nsim = n_sim) %>%
    as.matrix() %>%
    format_simulation_results(models$df_accum, "Poisson")

  # Simulate from Frailty model
  cat("  Simulating from Frailty model...\n")
  sim_frailty <- simulate(models$mod_frailty, nsim = n_sim) %>%
    as.matrix() %>%
    format_simulation_results(models$df_accum, "Frailty")

  # Simulate from Polya urn model
  cat("  Simulating from Polya urn model...\n")
  sim_polya <- lapply(1:n_sim, function(i) {
    simulate(models$mod_polya)
  }) %>%
    # Extract only number of successes (first column)
    lapply(function(x) {
      as.matrix(as.matrix(x)[, 1], ncol = 1)
    }) %>%
    do.call(cbind, .) %>%
    format_simulation_results(models$df_accum, "Polya urn")

  # Prepare empirical data
  cat("  Preparing empirical data...\n")
  df_accum_plot <- models$df_accum %>%
    dplyr::arrange(p_id, t_id) %>%
    dplyr::mutate(cumsum_nr_occur = cumsum(nr_occur), .by = p_id) %>%
    dplyr::group_by(t_id, cumsum_nr_occur) %>%
    dplyr::summarise(frequency = n(), .groups = 'drop') %>%
    dplyr::group_by(t_id) %>%
    dplyr::mutate(
      frequency_norm = frequency / sum(frequency),
      sim_type = "Empirical"
    ) %>%
    dplyr::ungroup() %>%
    rbind(sim_poisson, sim_frailty, sim_polya) %>%
    dplyr::mutate(
      sim_type = factor(sim_type, levels = sim_types, ordered = FALSE),
      time_label = sprintf("Year %02d", t_id) %>%
        stringr::str_wrap(width = 9)
    )
  cat(sprintf("Time elapsed: %s\n", format(Sys.time() - start_t)))

  return(df_accum_plot)
}


# ============================================================================
# VISUALIZATION
# ============================================================================

#' Plot cumulative distribution of events
#' @param name_dataset Dataset name
#' @param df Dataframe with simulation and empirical data
#' @param P Parameters list
#' @param col_values Named vector of colors
#' @param fill_values Named vector of fill colors
#' @param size_text Text size for plot
#' @param chosen_t_ids Time points to display
#' @param clip_to_emp Whether to clip simulations to empirical range
plot_cum_distr <- function(name_dataset, df, P, col_values, fill_values,
                                         size_text = 10,
                                         chosen_t_ids = c(1, 5, 10, 15, 20),
                                         clip_to_emp = TRUE) {

  # Filter to chosen time points
  df <- df %>%
    dplyr::filter(t_id %in% chosen_t_ids) %>%
    dplyr::mutate(dataset = name_dataset)

  # Clip simulated data to empirical range
  if (clip_to_emp) {
    max_emp <- df %>%
      dplyr::filter(sim_type == "Empirical") %>%
      dplyr::group_by(t_id) %>%
      dplyr::summarise(max = max(cumsum_nr_occur), .groups = 'drop') %>%
      dplyr::pull(max, t_id)

    df <- df %>%
      dplyr::group_by(t_id) %>%
      dplyr::filter(
        cumsum_nr_occur <= max_emp[names(max_emp) == as.character(unique(t_id))]
      ) %>%
      dplyr::ungroup()
  }

  # Create plot
  pl_distr <- df %>%
    dplyr::filter(sim_type != "Empirical") %>%
    ggplot() +
    # Simulated data as lines
    geom_line(
      aes(x = cumsum_nr_occur, y = frequency_norm, col = sim_type),
      alpha = 1,
      linewidth = 1.1
    ) +
    # Empirical data as segments
    geom_segment(
      data = df %>% dplyr::filter(sim_type == "Empirical"),
      aes(x = cumsum_nr_occur, xend = cumsum_nr_occur,
          y = 0, yend = frequency_norm, col = sim_type),
      alpha = 0.7,
      linewidth = 0.5
    ) +
    # Empirical data as points
    geom_point(
      data = df %>% dplyr::filter(sim_type == "Empirical"),
      aes(x = cumsum_nr_occur, y = frequency_norm, col = sim_type, size = t_id),
      # size = 1.5,
      alpha = 0.85
    ) +
    scale_size(range = c(1.5, 1), guide = "none") +
    ggh4x::facet_grid2(
      time_label ~ dataset,
      axes = "all",
      switch = "y",
      independent = "x",
      scales = "free"
    ) +
    scale_color_manual(
      name = "",
      values = col_values,
      breaks = names(col_values)
    ) +
    P$own_theme +
    scale_y_continuous(
      position = "right",
      labels = function(x) sprintf("%.2f", x),
      expand = expansion(mult = c(0.03, 0.08)),
      n.breaks = 3
    ) +
    scale_x_continuous(
      n.breaks = 4,
      expand = expansion(mult = c(0.015, 0.015))
    ) +
    labs(
      y = "Frequency (normalized per year)",
      x = "Number of cumulative events"
    ) +
    theme(legend.position = "top") +
    theme(
      # panel.grid.major = element_line(color = "gray90", size = 0.3),
      panel.grid.major = element_blank(),
      panel.border = element_rect(colour = "grey30", fill = NA, linewidth = 1),
      panel.spacing.y = unit(0.75, "lines"),
      panel.spacing.x = unit(0.01, "lines"),
      axis.text.x = element_text(size = size_text),
      strip.text.y = element_text(size = size_text),
      axis.text.y = element_text(size = size_text, vjust = 0),
      axis.title = element_text(size = size_text + 4),
      legend.text = element_text(size = size_text + 4)
    ) +
    theme(plot.margin = unit(c(0, 0.5, 0, 0), 'cm'))

  return(pl_distr)
}


#' Plot range of cumulative distribution of events
plot_cum_distr_range <- function(name_dataset, df, P, col_values, fill_values,
                             size_text = 10,
                             chosen_t_ids = c(1, 5, 10, 15, 20)){

  # Filter to chosen time points
  df <- df %>%
    dplyr::filter(t_id %in% chosen_t_ids)

  # Find range
  df_range <- df %>%
      dplyr::group_by(sim_type, t_id, time_label) %>%
      dplyr::summarise(min = min(cumsum_nr_occur),
                       mean = mean(cumsum_nr_occur),
                       median = median(cumsum_nr_occur),
                       max = max(cumsum_nr_occur),
                       .groups = 'drop') %>%
    dplyr::mutate(dataset = name_dataset)

  # Create range plot
  pl_range <- df_range %>%
    dplyr::mutate(sim_type = factor(sim_type,
                                    levels = rev(c("Empirical", "Poisson", "Frailty", "Polya urn")))) %>%
    ggplot(aes(y = time_label, col = sim_type)) +
    # Error bar caps at min and max
    geom_errorbar(
      aes(xmin = min, xmax = max),
      width = 0.5,
      alpha = 1,
      linewidth = 0.75,
      position = position_dodge(width = 0.8)
    ) +
    # Median as point
    geom_point(
      aes(x = median),
      alpha = 1,
      size = 2,
      position = position_dodge(width = 0.8)
    ) +
    ggh4x::facet_grid2(
      time_label ~ dataset,
      axes = "all",
      switch = "y",
      independent = "x",
      scales = "free"
    ) +
    scale_color_manual(
      name = "",
      values = col_values,
      breaks = names(col_values)
    ) +
    P$own_theme +
    scale_x_continuous(
      n.breaks = 4,
      expand = expansion(mult = c(0.015, 0.015))
    ) +
    labs(
      y = "",
      x = "Number of cumulative events"
    ) +
    theme(legend.position = "top") +
    theme(
      panel.grid.major = element_line(color = "gray90", size = 0.3),
      panel.border = element_rect(colour = "grey30", fill = NA, linewidth = 1),
      panel.spacing.y = unit(0.75, "lines"),
      panel.spacing.x = unit(0.01, "lines"),
      axis.text.x = element_text(size = size_text),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      strip.text.y = element_text(size = size_text),
      axis.title = element_text(size = size_text + 4),
      legend.text = element_text(size = size_text + 4)
    ) +
    theme(plot.margin = unit(c(0, 0.5, 0, 0), 'cm'))

  return(pl_range)
}


#' Plot SHP and HILDA's cumulative distributions together
plot_cum_distr_combo <- function(pl_distr_SHP, pl_distr_HILDA, filepath_base){
  pl_combo <- ((pl_distr_SHP +
                 theme(legend.margin = margin(t = 10, r = 20, b = 10, l = 250, unit = "pt")) +
                 (pl_distr_HILDA + theme(strip.text.y.left = element_blank()) +
                    theme(legend.margin = margin(t = 10, r = 20, b = 10, l = 250, unit = "pt"))) +
                 plot_layout(axis_titles = "collect")) + guide_area() +
                 plot_layout(guides = 'collect') + plot_layout(heights = c(1, .1)) )
  pl_combo

  filepath <- file.path(filepath_base, "figs", "cumulative_nr_events.pdf")
  save_plot(pl_combo, filepath, height = 170)
}

#' Plot range of SHP and HILDA's cumulative distributions together
plot_cum_range_combo <- function(pl_range_SHP, pl_range_HILDA, filepath_base){
  pl_combo <- ((pl_range_SHP +
                  theme(legend.margin = margin(t = 10, r = 20, b = 10, l = 250, unit = "pt")) +
                  (pl_range_HILDA + theme(strip.text.y.left = element_blank()) +
                     theme(legend.margin = margin(t = 10, r = 20, b = 10, l = 250, unit = "pt"))) +
                  plot_layout(axis_titles = "collect")) + guide_area() +
                 plot_layout(guides = 'collect') + plot_layout(heights = c(1, .1)) )
  pl_combo

  filepath <- file.path(filepath_base, "figs", "cumulative_nr_events_range.pdf")
  save_plot(pl_combo, filepath, height = 170)
}


#' Create animated cumulative distribution plot
#' @param name_dataset Dataset name
#' @param filepath_base Base file path for saving animation
#' @param df Dataframe with simulation and empirical data
#' @param P Parameters list
#' @param col_values Named vector of colors
#' @param fill_values Named vector of fill colors
#' @param clip_to_emp Whether to clip simulations to empirical range
plot_cumulative_distribution_animated <- function(name_dataset, filepath_base, df, P,
                                                  col_values, fill_values,
                                                  clip_to_emp = TRUE) {

  size_text <- 40

  # Clip simulated data to empirical range
  if (clip_to_emp) {
    max_emp <- df %>%
      dplyr::filter(sim_type == "Empirical") %>%
      dplyr::group_by(t_id) %>%
      dplyr::summarise(max = max(cumsum_nr_occur), .groups = 'drop') %>%
      dplyr::pull(max, t_id)

    df <- df %>%
      dplyr::group_by(t_id) %>%
      dplyr::filter(
        cumsum_nr_occur <= max_emp[names(max_emp) == as.character(unique(t_id))]
      ) %>%
      dplyr::ungroup()
  }

  # Ensure proper factor ordering
  df <- df %>%
    dplyr::mutate(
      sim_type = factor(sim_type, levels = levels(.data$sim_type), ordered = TRUE)
    )

  # Create smoothed lines for simulated data
  df_smooth <- df %>%
    dplyr::filter(sim_type != "Empirical") %>%
    dplyr::arrange(sim_type, t_id, cumsum_nr_occur) %>%
    dplyr::group_by(t_id, sim_type) %>%
    do({
      # Create finer grid for interpolation
      x_new <- seq(min(.$cumsum_nr_occur), max(.$cumsum_nr_occur), length.out = 1000)
      # Spline interpolation
      spline_fit <- spline(.$cumsum_nr_occur, .$frequency_norm,
                           xout = x_new, method = "fmm")
      data.frame(
        cumsum_nr_occur = spline_fit$x,
        frequency_norm_smooth = spline_fit$y,
        sim_type = unique(.$sim_type)
      )
    }) %>%
    dplyr::ungroup()

  # Create base plot
  pl_distr <- df_smooth %>%
    ggplot() +
    # Smoothed lines for simulated data
    geom_line(
      aes(x = cumsum_nr_occur, y = frequency_norm_smooth, col = sim_type),
      alpha = 1,
      linewidth = 1.5
    ) +
    # Empirical data as bars
    geom_bar(
      data = df %>% dplyr::filter(sim_type == "Empirical"),
      aes(x = cumsum_nr_occur, y = frequency_norm),
      stat = 'identity',
      col = col_values[names(col_values) == "Empirical"],
      fill = fill_values[names(fill_values) == "Empirical"],
      position = "identity",
      alpha = 0.5,
      linewidth = 0.5
    ) +
    scale_color_manual(
      name = "",
      values = col_values[names(col_values) != "Empirical"]
    ) +
    scale_fill_manual(
      name = "",
      values = fill_values[names(fill_values) != "Empirical"]
    ) +
    P$own_theme +
    scale_y_continuous(
      expand = expansion(mult = c(0.0, 0.02)),
      n.breaks = 3
    ) +
    scale_x_continuous(
      n.breaks = 4,
      limits = c(-0.5, 75),
      expand = expansion(mult = c(0, 0))
    ) +
    labs(
      y = "Frequency (normalized)",
      x = "Number of cumulative adverse life events"
    ) +
    theme(legend.position = "bottom") +
    theme(
      panel.border = element_blank(),
      panel.spacing.y = unit(0.75, "lines"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(0.01, "lines"),
      axis.text.x = element_text(size = size_text),
      strip.text.y = element_text(size = size_text),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title = element_text(size = size_text + 4),
      title = element_text(size = size_text + 10),
      legend.text = element_text(size = size_text + 4)
    ) +
    theme(plot.margin = unit(c(0, 0.5, 0, 0), 'cm'))

  # Create animation
  anim <- pl_distr +
    gganimate::transition_time(t_id) +
    gganimate::view_follow(fixed_x = TRUE) +
    labs(title = 'Year: {frame_time}') +
    theme(title = element_text(size = size_text + 10))

  # Render animation
  gganimate::animate(
    anim,
    nframes = 100,
    height = 12,
    width = 18,
    units = "cm",
    res = 300
  )

  # Save animation
  filepath <- file.path(filepath_base, "figs", name_dataset, "cumulative_nr_events.gif")
  gganimate::anim_save(filepath)

  return(filepath)
}


# ============================================================================
# HEAVY-TAILED DISTRIBUTION ANALYSIS
# ============================================================================

#' Fit heavy-tailed distributions to cumulative event counts
#' @param counts Vector of cumulative event counts
#' @param set_xmin Optional fixed xmin value
fit_heavy_tailed_distributions <- function(counts, set_xmin = NULL) {

  # 1. Fit power-law distribution
  cat("  Fitting power-law distribution...\n")
  pl_model <- poweRlaw::displ$new(counts)

  if (is.null(set_xmin)) {
    est <- poweRlaw::estimate_xmin(pl_model)
  } else {
    est <- set_xmin
  }

  pl_model$setXmin(est)
  pl_model$setPars(poweRlaw::estimate_pars(pl_model))
  cat(sprintf("    Power-law: xmin = %d, alpha = %.2f\n",
              pl_model$getXmin(), pl_model$pars))

  # 2. Fit lognormal distribution
  cat("  Fitting lognormal distribution...\n")
  ln_model <- poweRlaw::dislnorm$new(counts)

  if (is.null(set_xmin)) {
    ln_est <- poweRlaw::estimate_xmin(ln_model)
  } else {
    ln_est <- set_xmin
  }

  ln_model$setXmin(ln_est)
  ln_model$setPars(poweRlaw::estimate_pars(ln_model))
  cat(sprintf("    Lognormal: xmin = %d, meanlog = %.2f, sdlog = %.2f\n",
              ln_model$getXmin(), ln_model$pars[1], ln_model$pars[2]))

  # 3. Fit exponential distribution
  cat("  Fitting exponential distribution...\n")
  exp_model <- poweRlaw::disexp$new(counts)

  if (is.null(set_xmin)) {
    exp_est <- poweRlaw::estimate_xmin(exp_model)
  } else {
    exp_est <- set_xmin
  }

  exp_model$setXmin(exp_est)
  exp_model$setPars(poweRlaw::estimate_pars(exp_model))
  cat(sprintf("    Exponential: xmin = %d, rate = %.2f\n",
              exp_model$getXmin(), exp_model$pars))

  # 4. Fit Poisson distribution
  cat("  Fitting Poisson distribution...\n")
  pois_model <- poweRlaw::dispois$new(counts)

  if (is.null(set_xmin)) {
    pois_est <- poweRlaw::estimate_xmin(pois_model)
  } else {
    pois_est <- set_xmin
  }

  pois_model$setXmin(pois_est)
  pois_model$setPars(poweRlaw::estimate_pars(pois_model))
  cat(sprintf("    Poisson: xmin = %d, rate = %.2f\n",
              pois_model$getXmin(), pois_model$pars))

  # 5. Create diagnostic plot
  plot(pl_model, main = "Fitted Distributions", xlab = "Value", ylab = "CDF", pch = 19)
  lines(pl_model, col = "red", lwd = 3, lty = 1)
  lines(ln_model, col = "blue", lwd = 3, lty = 2)
  lines(exp_model, col = "green", lwd = 3, lty = 3)
  lines(pois_model, col = "magenta", lwd = 3, lty = 4)

  legend(
    "bottomleft",
    legend = c("Power law", "Lognormal", "Exponential", "Poisson"),
    col = c("red", "blue", "green", "magenta"),
    lty = c(1, 2, 3, 4),
    lwd = 2
  )

  return(list(
    pl_model = pl_model,
    ln_model = ln_model,
    exp_model = exp_model,
    pois_model = pois_model
  ))
}


#' Compare two heavy-tailed distributions
#' @param model1_ First model
#' @param model2_ Second model
#' @param model1_name Name of first model
#' @param model2_name Name of second model
#' @param choose_xmin How to choose xmin ("model1", "model2", or "median")
compare_heavy_tailed_distributions <- function(model1_, model2_,
                                               model1_name, model2_name,
                                               choose_xmin = "median") {

  # Copy models to avoid overwriting (shallow copies don't work)
  model1 <- model1_$copy()
  model2 <- model2_$copy()

  # Set same xmin for both models
  if (choose_xmin == "median") {
    xmin_min <- round(median(c(model1$getXmin(), model2$getXmin())))
    model1$setXmin(xmin_min)
    model1$setPars(poweRlaw::estimate_pars(model1))
    model2$setXmin(xmin_min)
    model2$setPars(poweRlaw::estimate_pars(model2))

  } else if (choose_xmin == "model1") {
    model2$setXmin(model1$getXmin())
    model2$setPars(poweRlaw::estimate_pars(model2))

  } else if (choose_xmin == "model2") {
    model1$setXmin(model2$getXmin())
    model1$setPars(poweRlaw::estimate_pars(model1))
  }

  # Compare models
  comparison_1_to_2 <- poweRlaw::compare_distributions(model1, model2)
  comparison_2_to_1 <- poweRlaw::compare_distributions(model2, model1)

  # Format results
  res <- data.frame(
    loglik_ratio = comparison_1_to_2$test_statistic,
    p_two_sided = comparison_1_to_2$p_two_sided,
    p_1_better_than_2 = comparison_1_to_2$p_one_sided,
    p_2_better_than_1 = comparison_2_to_1$p_one_sided
  ) %>%
    tidyr::gather() %>%
    dplyr::mutate(key = sprintf("%s_vs_%s_%s", model1_name, model2_name, key))

  return(res)
}


#' Run complete heavy-tailed distribution analysis
#' @param name_dataset Dataset name
#' @param df_accum Accumulation dataframe
#' @param datalist Data list
#' @param P Parameters list
run_heavy_tails <- function(name_dataset, df_accum, datalist, P) {

  cat("\n=== Running heavy-tailed distribution analysis ===\n")

  # 1. Calculate cumulative counts at 20 years
  counts <- df_accum %>%
    dplyr::arrange(p_id, t_id) %>%
    dplyr::group_by(p_id) %>%
    dplyr::summarise(nr_occur_cumsum = sum(nr_occur), .groups = 'drop')

  cat(sprintf("  Number of people with zero events: %d\n",
              nrow(counts %>% dplyr::filter(nr_occur_cumsum == 0))))

  # Remove zero counts (required for power-law fitting)
  counts <- counts %>%
    dplyr::filter(nr_occur_cumsum != 0) %>%
    dplyr::pull(nr_occur_cumsum)

  # 2. Calculate summary statistics
  cat("\n  Summary statistics:\n")
  min_cumsum <- min(counts)
  max_cumsum <- max(counts)
  mean_cumsum <- mean(counts)
  median_cumsum <- median(counts)
  var_to_mean <- var(counts) / mean(counts)

  cat(sprintf("    Min: %d\n", min_cumsum))
  cat(sprintf("    Max: %d\n", max_cumsum))
  cat(sprintf("    Mean: %.2f\n", mean_cumsum))
  cat(sprintf("    Median: %.2f\n", median_cumsum))
  cat(sprintf("    Variance-to-mean ratio: %.2f\n", var_to_mean))

  # 3. Fit heavy-tailed distributions
  cat("\n  Fitting distributions:\n")
  fitted_models <- fit_heavy_tailed_distributions(counts)

  # 4. Create plot data
  ecdf_data <- ecdf(counts)
  x_vals <- unique(sort(counts))
  y_vals <- 1 - ecdf_data(x_vals)  # Complementary CDF (CCDF)

  # Empirical data
  df_emp <- data.frame(
    x = unique(sort(counts)),
    y = 1 - ecdf_data(x_vals)
  ) %>%
    dplyr::mutate(type = name_dataset)

  # Fitted line data
  pl_model_line <- lines(fitted_models$pl_model, draw = FALSE)
  ln_model_line <- lines(fitted_models$ln_model, draw = FALSE)
  exp_model_line <- lines(fitted_models$exp_model, draw = FALSE)
  pois_model_line <- lines(fitted_models$pois_model, draw = FALSE)

  df_fitted <- rbind(
    pl_model_line %>% dplyr::mutate(type = "Power-law"),
    ln_model_line %>% dplyr::mutate(type = "Log-normal"),
    exp_model_line %>% dplyr::mutate(type = "Exponential"),
    pois_model_line %>% dplyr::mutate(type = "Poisson")
  )

  # 5. Create visualization
  pl <- df_fitted %>%
    dplyr::mutate(dataset = name_dataset) %>%
    ggplot() +
    geom_point(
      data = df_emp,
      aes(x = x, y = y),
      size = 1,
      col = 'red', alpha = .9
    ) +
    geom_line(
      aes(x = x, y = y, col = type),
      linewidth = .65,
      # linetype = 'dashed',
      alpha = 1
    ) +

    ggh4x::facet_grid2(. ~ dataset) +
    scale_x_continuous(
      trans = "log10",
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x)),
      limits = range(df_emp$x)
    ) +
    scale_y_continuous(
      trans = "log10",
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    labs(x = "Cumulative count of events", y = "Probability") +
    viridis::scale_color_viridis(discrete = TRUE, name = "", end = 0.95) +
    P$own_theme +
    theme(
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.text = element_text(size = 10)
    )

  # 6. Compare distributions
  cat("\n  Comparing distributions:\n")
  comparison_results <- rbind(
    compare_heavy_tailed_distributions(
      fitted_models$ln_model, fitted_models$pl_model,
      "lognormal", "powerlaw"
    ),
    compare_heavy_tailed_distributions(
      fitted_models$ln_model, fitted_models$exp_model,
      "lognormal", "exponential"
    ),
    compare_heavy_tailed_distributions(
      fitted_models$ln_model, fitted_models$pois_model,
      "lognormal", "poisson"
    )
  ) %>%
    tidyr::spread(key, value)

  # 7. Format results for table
  df_est <- data.frame(
    dataset = name_dataset,
    pois_xmin = fitted_models$pois_model$xmin,
    pois_rate = fitted_models$pois_model$pars,
    exp_xmin = fitted_models$exp_model$xmin,
    exp_rate = fitted_models$exp_model$pars,
    lognormal_xmin = fitted_models$ln_model$xmin,
    lognormal_meanlog = fitted_models$ln_model$pars[1],
    lognormal_sdlog = fitted_models$ln_model$pars[2],
    powerlaw_xmin = fitted_models$pl_model$xmin,
    powerlaw_alpha = fitted_models$pl_model$pars
  ) %>%
    dplyr::mutate_at(
      setdiff(colnames(.), c("dataset", "pois_xmin", "exp_xmin",
                             "lognormal_xmin", "powerlaw_xmin")),
      ~ sprintf("%.2f", .)
    ) %>%
    dplyr::mutate(
      lognormal_vs_poisson = sprintf(
        "%.2f ($p$ = %.2f)",
        comparison_results[["lognormal_vs_poisson_loglik_ratio"]],
        comparison_results[["lognormal_vs_poisson_p_two_sided"]]
      ),
      lognormal_vs_exp = sprintf(
        "%.2f ($p$ = %.2f)",
        comparison_results[["lognormal_vs_exponential_loglik_ratio"]],
        comparison_results[["lognormal_vs_exponential_p_two_sided"]]
      ),
      lognormal_vs_powerlaw = sprintf(
        "%.2f ($p$ = %.2f)",
        comparison_results[["lognormal_vs_powerlaw_loglik_ratio"]],
        comparison_results[["lognormal_vs_powerlaw_p_two_sided"]]
      )
    ) %>%
    dplyr::mutate_at(
      c("lognormal_vs_poisson", "lognormal_vs_exp", "lognormal_vs_powerlaw"),
      ~ stringr::str_replace_all(., "= 0.00", "$<$ .01")
    )

  return(list(
    df_est = df_est,
    pl = pl
  ))
}


#' Create LaTeX table for heavy-tailed distribution analysis
#' @param SHP_tails Results from SHP analysis
#' @param HILDA_tails Results from HILDA analysis
create_heavy_tails_table <- function(SHP_tails, HILDA_tails) {

  rbind(
    SHP_tails$df_est,
    HILDA_tails$df_est
  ) %>%
    kableExtra::kbl(
      format = "latex",
      escape = FALSE,
      caption = "Model estimates of distributions fit to twenty-year cumulative adverse life event counts (Source: Swiss Household Panel, SHP and Household, Income and Labour Dynamics in Australia Survey, HILDA). The log-likelihood ratio compares the fit of the log-normal model to an alternative model (Poisson, exponential, power-law). Positive values indicate the log-normal model fits better, whereas negative values favour the alternative model. The two-sided p-value tests whether the difference in fit is statistically significant.",
      booktabs = TRUE,
      label = "heavytails",
      col.names = NULL
    ) %>%
    kableExtra::add_header_above(
      c(" ", "$x_{min}$" = 1, "$\\\\lambda$" = 1, "$x_{min}$" = 1,
        "$\\\\lambda$" = 1, "$x_{min}$" = 1, "$\\\\mu$" = 1, "$\\\\sigma$" = 1,
        "$x_{min}$", "$\\\\alpha$", "Poisson", "Exponential", "Power-law"),
      escape = FALSE
    ) %>%
    kableExtra::add_header_above(
      c(" ", "Poisson" = 2, "Exponential" = 2, "Log-normal" = 3,
        "Power-law" = 2, "Log-Likelihood Ratio cf. Log-normal" = 3),
      escape = FALSE
    ) %>%
    kableExtra::kable_styling(latex_options = "scale_down")
}

#' Plot heavy tails of SHP and HILDA together
plot_heavy_tails <- function(pl_SHP, pl_HILDA, filepath_base){

  pl <- ((pl_SHP + theme(legend.margin = margin(t = 10, r = 20, b = 10, l = 225, unit = "pt")) +
           (pl_HILDA + theme(legend.margin = margin(t = 10, r = 20, b = 10, l = 225, unit = "pt"))) +
           plot_layout(axis_titles = "collect")) + guide_area() +
          plot_layout(guides = 'collect') + plot_layout(heights = c(1, .1)) )


  filepath_image <- file.path(filepath_base, "figs", sprintf("heavy_tails.pdf"))
  save_plot(pl, filepath_image, height = 100)
}
