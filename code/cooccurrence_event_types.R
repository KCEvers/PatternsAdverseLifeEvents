##### CO-OCCURRENCE ANALYSIS - REFACTORED #####

#' Main wrapper function for co-occurrence analysis
#' @param name_dataset Dataset name (e.g., "SHP", "HILDA")
#' @param filepath_base Base file path
#' @param P Parameters list
#' @param chosen_leads Lag values to analyze (e.g., c(0, 1))
#' @param run_unadjusted Whether to run unadjusted sensitivity analysis
#' @param rerun Whether to rerun analysis
run_cooccurrence <- function(name_dataset, filepath_base, P,
                             chosen_leads = c(0, 1),
                             run_unadjusted = TRUE,
                             rerun = FALSE) {

  # 1. Prepare data
  cat("\n=== Preparing data ===\n")
  event_dict = P$event_dicts[[name_dataset]]
  datalist <- prepare_data(name_dataset, filepath_base, event_dict)

  # 2. Format event data
  cat("\n=== Formatting event data ===\n")
  df_events <- format_event_data(datalist)
  negevent_names <- get_event_names(df_events)

  # 3. Run models for each lag
  cat("\n=== Running models ===\n")
  model_dfs <- list()

  for (chosen_lead in chosen_leads) {
    cat(sprintf("\n--- Analyzing lag-%d ---\n", chosen_lead))

    # Run fully-adjusted models
    model_dfs[[as.character(chosen_lead)]] <- run_models_for_lead(
      df_events = df_events,
      negevent_names = negevent_names,
      event_dict = event_dict,
      chosen_lead = chosen_lead,
      name_dataset = name_dataset,
      datalist = datalist,
      run_unadjusted = run_unadjusted,
      rerun = rerun
    )

  }

  # 4. Return consolidated results
  return(list(
    datalist = datalist,
    df_events = df_events,
    negevent_names = negevent_names,
    model_dfs = model_dfs
  ))
}


# ============================================================================
# DATA PREPARATION
# ============================================================================

#' Format event data for analysis
format_event_data <- function(datalist, keep_h_id = TRUE) {
  df_events <- datalist$df_per_event_pp_py %>%
    dplyr::filter(valence == "negative") %>%
    dplyr::mutate(occur = ifelse(nr_missing_occur > 0, NA, nr_occur)) %>%
    dplyr::select(p_id, crosswave_h_id, wave_nr, event, occur)

  if (!keep_h_id){
    df_events <- dplyr::select(df_events, -crosswave_h_id)
  }

  df_events <- df_events %>%
    tidyr::spread(event, occur) %>%
    tidyr::complete(p_id, wave_nr) %>%
    dplyr::arrange(p_id, wave_nr)

  return(df_events)
}

#' Extract event names from formatted data
get_event_names <- function(df_events) {
  negevent_names <- sort(setdiff(
    colnames(df_events),
    c("p_id", "crosswave_h_id", "wave_nr")
  ))
  return(negevent_names)
}


# ============================================================================
# MODEL FITTING
# ============================================================================

#' Run all models for a specific lead/lag
run_models_for_lead <- function(df_events, negevent_names, event_dict, chosen_lead,
                                name_dataset, datalist, run_unadjusted, rerun) {

  # File paths
  filepath_models <- file.path(
    datalist$filepath_deriv,
    sprintf("cooccurrence_event_types_lead%d.RDS", chosen_lead)
  )

  # Fit models (or load)
  if (!file.exists(filepath_models) || rerun) {
    cat("Fitting models...\n")
    start_t <- Sys.time()

    model_list <- fit_full_models(
      df_events = df_events,
      negevent_names = negevent_names,
      chosen_lead = chosen_lead,
      name_dataset = name_dataset,
      event_dict = event_dict
    )

    cat(sprintf("Time elapsed: %s\n", format(Sys.time() - start_t)))
    saveRDS(model_list, filepath_models)
  } else {
    cat("Loading saved models...\n")
    model_list <- readRDS(filepath_models)
  }

  # Optionally run unadjusted models on same sample
  if (run_unadjusted) {
      cat("Running unadjusted sensitivity analysis...\n")
      run_unadjusted_sensitivity(
        model_list = model_list, negevent_names = negevent_names,
        chosen_lead = chosen_lead,
        name_dataset = name_dataset, datalist = datalist, P = P
      )
  }

  # Extract sample sizes and estimates
  model_df <- lapply(model_list, `[[`, "model_df") %>%
    do.call(rbind, .) %>% as.data.frame()

  return(model_df)
}


#' Fit fully-adjusted models (all events as predictors)
fit_full_models <- function(df_events, negevent_names, chosen_lead,
                            name_dataset, event_dict) {

  model_list <- foreach::foreach(
    i = 1:length(negevent_names),
    .packages = c("glmmTMB", "dplyr", "purrr"),
    .export = c("negevent_names", "df_events", "chosen_lead",
                "name_dataset", "event_dict",
                "get_valid_predictors", "prepare_model_data",
                "build_model_formulas", "fit_model_variants",
                "summarise_model", "get_plot_df",
                "extract_timepoints_per_person")
  ) %dopar% {

    outcome <- negevent_names[i]
    predictors <- get_valid_predictors(
      outcome = outcome,
      all_events = negevent_names,
      chosen_lead = chosen_lead,
      name_dataset = name_dataset,
      event_dict = event_dict
    )

    # Prepare data (lag if needed)
    df_model <- prepare_model_data(
      df_events = df_events,
      outcome = outcome,
      chosen_lead = chosen_lead
    )

    # Build formulas
    formulas <- build_model_formulas(
      outcome = outcome,
      predictors = predictors,
      chosen_lead = chosen_lead
    )

    # Fit three model specifications
    model <- fit_model_variants(formulas, df_model)

    # Summarise
    model_summ <- summarise_model(model)

    # Extract model estimates
    model_df <- get_plot_df(model_summ, outcome)

    # Extract sample sizes from models
    nr_obs <- purrr::map_vec(ranef(model)$cond, nrow)

    # Extract number of timepoints per person
    nr_timepoints <- extract_timepoints_per_person(model)

    return(list(outcome = outcome,
                model = model,
                model_summ = model_summ,
                model_df = model_df,
                nr_obs = nr_obs,
                nr_timepoints = nr_timepoints))
  }

  return(model_list)
}

# Get time points per person in mixed-effects model
extract_timepoints_per_person <- function(model) {
  if ("p_id" %in% colnames(model$frame)) {
    model$frame %>%
      dplyr::group_by(p_id) %>%
      dplyr::summarise(nr_waves = n(), .groups = 'drop') %>%
      dplyr::group_by(nr_waves) %>%
      dplyr::summarise(nr_people = n(), .groups = 'drop')
  }
}

# Create dataframe with model estimates
get_plot_df = function(model_summ, outcome){

  # Get fixed effects estimates
  model_est = model_summ$df_est %>%
    filter(effect == "fixed") %>%
    dplyr::rename(lower = "2.5 %", upper = "97.5 %") %>%
    dplyr::select(all_of(c("term", "estimate", "lower", "upper"))) %>%
    dplyr::rename(statistic = term) %>%
    dplyr::mutate(subplot = "fixed")

  # Get SD of random intercepts - note that not all models have a household intercept
  sd_intercept = model_summ$df_est %>%
    filter(effect == "ran_pars") %>%
    dplyr::rename(lower = "2.5 %", upper = "97.5 %") %>%
    dplyr::mutate(term = paste0(group, "_", term) %>%
                    stringr::str_replace_all(., "\\(Intercept\\)", "intercept") %>%
                    stringr::str_replace_all(., "__", "_")) %>%
    dplyr::select(all_of(c("term", "estimate", "lower", "upper"))) %>%
    dplyr::rename(statistic = term) %>%
    dplyr::mutate(subplot = "sd_intercept")

  # Get variance partitioning statistics
  var_part = rbind(
    model_summ[["ICC"]] %>% as.data.frame() %>% dplyr::select("ICC_adjusted") %>% t(),
    model_summ[["R2"]] %>% as.data.frame() %>% dplyr::select(-optional) %>% t()
  ) %>% magrittr::set_colnames("estimate") %>% as.data.frame() %>%
    tibble::rownames_to_column("statistic") %>%
    dplyr::mutate(subplot = "var_part")

  model_df = cbind(response = outcome,
                   dplyr::bind_rows(model_est, sd_intercept, var_part))

  return(model_df)
}


##### CO-OCCURRENCE VISUALIZATION - REFACTORED #####

# ============================================================================
# MAIN PLOTTING FUNCTIONS
# ============================================================================

#' Plot co-occurrence heatmap for all events
#' @param name_dataset Dataset name (e.g., "SHP", "HILDA")
#' @param datalist Data list containing paths and data
#' @param P Parameters list
#' @param model_df Model results dataframe
#' @param chosen_lead Lag value (0 or 1)
plot_cooccur <- function(name_dataset, datalist, P, model_df, chosen_lead,
                         color_lab = "Adjusted odds ratio",
                         suffix = "") {

  # 1. Prepare plot data
  plot_df <- prepare_cooccur_plot_data(model_df, chosen_lead, all_events = TRUE)
  ylab <- ifelse(chosen_lead == 0, " Outcome", "Outcome (next year)")

  # 2. Create fixed effects dataframe
  fixed_df <- create_fixed_effects_df(
    plot_df = plot_df,
    name_dataset = name_dataset,
    datalist = datalist,
    P = P
  )

  # 3. Apply dataset-specific adjustments
  if (name_dataset == "SHP") {
    fixed_df <- fixed_df %>%
      dplyr::mutate_at(
        c("response"),
        ~ forcats::fct_relevel(.x, "Other or unspecified illness or accident")
      )
  }

  # 4. Set plot parameters
  plot_params <- get_plot_parameters(
    name_dataset = name_dataset,
    selected = FALSE
  )

  # 5. Create plot
  pl_fixed <- create_odds_ratio_heatmap(
    fixed_df = fixed_df,
    plot_params = plot_params,
    P = P,
    ylab = ylab,
    color_lab = color_lab
  )

  # 6. Save plot
  save_cooccur_plot(
    plot = pl_fixed,
    datalist = datalist,
    name_dataset = name_dataset,
    chosen_lead = chosen_lead,
    suffix = suffix
  )

  return(pl_fixed)
}


#' Plot co-occurrence heatmap for selected events only
#' @param name_dataset Dataset name (e.g., "SHP", "HILDA")
#' @param datalist Data list containing paths and data
#' @param P Parameters list
#' @param model_df Model results dataframe
#' @param chosen_lead Lag value (0 or 1)
plot_selected_cooccur <- function(name_dataset, datalist, P, model_df, chosen_lead) {

  # Define selected events
  selected_events <- c(
    "jailed",
    "separated_from_spouse",
    "victim_physical_violence",
    "fired",
    "injury_illness_self"
  )

  # 1. Prepare plot data
  plot_df <- prepare_cooccur_plot_data(model_df, chosen_lead, all_events = FALSE)
  ylab <- ifelse(chosen_lead == 0, " Outcome", "Outcome (next year)")

  # 2. Create fixed effects dataframe
  fixed_df <- create_fixed_effects_df(
    plot_df = plot_df,
    name_dataset = name_dataset,
    datalist = datalist,
    P = P
  )

  # 3. Filter to selected events only
  fixed_df <- fixed_df %>%
    dplyr::filter(
      .data$response %in% selected_events,
      .data$statistic %in% selected_events
    )

  # 4. Set plot parameters
  plot_params <- get_plot_parameters(
    name_dataset = name_dataset,
    selected = TRUE
  )

  # 5. Create plot
  pl_fixed <- create_odds_ratio_heatmap(
    fixed_df = fixed_df,
    plot_params = plot_params,
    P = P,
    ylab = ylab
  )

  # 6. Save plot
  save_cooccur_plot(
    plot = pl_fixed,
    datalist = datalist,
    name_dataset = name_dataset,
    chosen_lead = chosen_lead,
    suffix = "_selected"
  )

  return(pl_fixed)
}


#' Plot difference fully adjusted versus unadjusted co-occurrence heatmap for all events
#' @param name_dataset Dataset name (e.g., "SHP", "HILDA")
#' @param datalist Data list containing paths and data
#' @param P Parameters list
#' @param unadjusted_df Model results dataframe
#' @param chosen_lead Lag value (0 or 1)
plot_cooccur_unadjusted <- function(name_dataset, datalist, P, unadjusted_df,
                                    chosen_lead) {

  # 1. Prepare plot data
  ylab <- ifelse(chosen_lead == 0, "Outcome", "Outcome (next year)")

  # Get odds ratios
  plot_df <- unadjusted_df %>%
    # Exponentiate to obtain odds ratio
    dplyr::mutate_at(c("estimate", "lower", "upper",
                       "full_estimate", "full_lower", "full_upper"),
                     ~ round(exp(.), 2)) %>%
    # Determine statistical significance
    dplyr::mutate(
      excludes_1 = (lower < 1 & upper < 1) | (lower > 1 & upper > 1),
      full_excludes_1 = (full_lower < 1 & full_upper < 1) | (full_lower > 1 & full_upper > 1),
      diff_OR = estimate - full_estimate) %>%
    dplyr::mutate(
      label = ifelse(is.na(diff_OR), "",
                                 sprintf("%.2f\n", diff_OR)),
                  label_diff = ifelse(is.na(diff_OR), "",
                                      sprintf("\n%.2f - %.2f",
                                       estimate, full_estimate
                                       ))
    ) %>%

  # plot_df <- unadjusted_df %>%
  #   dplyr::mutate(
  #     pct_change = 100 * (unadjusted_OR - full_OR) / full_OR,
  #     label = ifelse(is.na(unadjusted_OR), "",
  #                    sprintf("%.2f\n", pct_change)),
  #     label_diff = ifelse(is.na(unadjusted_OR), "",
  #                         sprintf("\n%.2f → %.2f",
  #                                 unadjusted_OR, full_OR
  #                         ))

      # label = ifelse(is.na(unadjusted_minus_full_OR), "",
      #                            sprintf("%.2f\n", unadjusted_minus_full_OR)),
      #             label_diff = ifelse(is.na(unadjusted_minus_full_OR), "",
      #                                 sprintf("\n%.2f - %.2f",
      #                                  unadjusted_OR, full_OR
      #                                  ))

      # ) %>%
    # Rename events with descriptions
    recode_events(
      .,
      c("response", "statistic"),
      P$event_dicts[[name_dataset]],
      datalist$df_per_event
    )


  # 3. Apply dataset-specific adjustments
  if (name_dataset == "SHP") {
    plot_df <- plot_df %>%
      dplyr::mutate_at(
        c("outcome"),
        ~ forcats::fct_relevel(.x, "Other or unspecified illness or accident")
      )

    # Remove from predictors
    plot_df <- plot_df %>%
      dplyr::filter(statistic != "Other or unspecified illness or accident")
  }

  # 4. Set plot parameters
  plot_params <- get_plot_parameters(
    name_dataset = name_dataset,
    selected = FALSE
  )

  # 5. Create plot
  # Build theme
  adapt_theme <- build_heatmap_theme(plot_params, P)

  # Create base plot with tiles
  pl_fixed <- plot_df %>%
    ggplot() +
    geom_tile(
      aes(x = statistic, y = response, col = diff_OR, fill = diff_OR),
      alpha = 0.95,
      linewidth = 0.5
    )

  # Add scales and labels
  color_lab <- "Odds ratio difference (unadjusted - fully adjusted)"
  pl_fixed <- pl_fixed +
    scale_fill_gradient2(
      name = color_lab,
      low = "blue",      # negative values
      mid = "white",     # zero/midpoint
      high = "red",      # positive values
      midpoint = 0,      # where the color transitions
      na.value = "white"
    ) +
    scale_color_gradient2(
      name = color_lab,
      low = "blue",      # negative values
      mid = "white",     # zero/midpoint
      high = "red",      # positive values
      midpoint = 0,      # where the color transitions
      na.value = "white"
    ) +
    guides(
      col = NULL,
      fill = guide_colourbar(
        title.position = "top",
        theme = theme(
          legend.key.width = unit(plot_params$width_legendbar, "lines"),
          legend.key.height = unit(plot_params$height_legendbar, "lines")
        )
      )
    ) +
    labs(y = ylab, x = "Predictor") +
    adapt_theme +
    scale_x_discrete(
      labels = function(x) sapply(x, wrap_equally, max_chars = plot_params$max_chars),
      limits = rev,
      position = 'top',
      expand = expansion(add = c(0, 0))
    ) +
    scale_y_discrete(
      labels = function(x) sapply(x, wrap_equally, max_chars = plot_params$max_chars),
      expand = expansion(add = c(0, 0))
    )

  pl_fixed <- pl_fixed +
    # Subscript
    geom_text(
      aes(x = statistic, y = response,
          label = label_diff),
      hjust = 0.5,
      vjust = plot_params$vjust,
      size = plot_params$size_text,
      family = P$font_family
    ) +
    # Point estimates (bold)
    geom_text(
      aes(x = statistic, y = response,
          label = label),
      fontface = "bold",
      hjust = 0.5,
      vjust = plot_params$vjust,
      size = plot_params$size_text + .8,
      family = P$font_family
    )

  # 6. Save plot
  save_cooccur_plot(
    plot = pl_fixed,
    datalist = datalist,
    name_dataset = name_dataset,
    chosen_lead = chosen_lead,
    suffix = "_unadjusted_diff"
  )

  return(pl_fixed)
}



# ============================================================================
# DATA PREPARATION
# ============================================================================

#' Prepare data for co-occurrence plotting
prepare_cooccur_plot_data <- function(model_df, chosen_lead, all_events = TRUE) {

  plot_df <- model_df

  if (chosen_lead == 0) {
    # Remove intercept and self-comparisons
    plot_df <- plot_df %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        statistic = ifelse(statistic == "(Intercept)", response, statistic)
      ) %>%
      dplyr::ungroup() %>%
      dplyr::filter(response != statistic)

  } else if (chosen_lead == 1) {
    # Remove intercept only
    plot_df <- plot_df %>%
      dplyr::filter(statistic != "(Intercept)")
  }

  return(plot_df)
}


#' Create fixed effects dataframe with odds ratios
create_fixed_effects_df <- function(plot_df, name_dataset, datalist, P) {

  fixed_df <- plot_df %>%
    dplyr::filter(subplot == "fixed") %>%
    # Exponentiate to obtain odds ratio
    dplyr::mutate_at(c("estimate", "lower", "upper"), ~ round(exp(.), 2)) %>%
    # Determine statistical significance
    dplyr::mutate(
      excludes_1 = (lower < 1 & upper < 1) | (lower > 1 & upper > 1)
    ) %>%
    # Create display labels
    dplyr::mutate(
      fill = ifelse(excludes_1, estimate, NA),
      label = sprintf("%.2f\n[%.2f, %.2f]", estimate, lower, upper),
      est = sprintf("%.2f\n", estimate),
      ci = sprintf("\n[%.2f, %.2f]", lower, upper)
    ) %>%
    # Rename events with descriptions
    recode_events(
      .,
      c("response", "statistic"),
      P$event_dicts[[name_dataset]],
      datalist$df_per_event
    )

  return(fixed_df)
}



# ============================================================================
# PLOT PARAMETERS
# ============================================================================

#' Get plot-specific parameters
get_plot_parameters <- function(name_dataset, selected = FALSE) {

  if (selected) {
    # Parameters for selected events plot
    size_text_theme <- ifelse(name_dataset == "SHP", 9.5 * 1.6, 9.5 * 1.6)
    size_text <- ifelse(name_dataset == "SHP", 1.8 * 2.6, 1.5 * 2.6)
    width_legendbar <- 25
    height_legendbar <- 1
    vjust <- 0.6
    max_chars <- 15
    axis_margin_top <- 20
    axis_margin_right <- 20

  } else {
    # Parameters for all events plot
    size_text_base <- 9.5
    size_text_theme <- size_text_base * 1.6
    size_text <- ifelse(name_dataset == "SHP", 1.8, 1.5) - 0.1
    width_legendbar <- 10
    height_legendbar <- 1
    vjust <- 0.6
    max_chars <- 31
    axis_margin_top <- 0
    axis_margin_right <- 0
  }

  return(list(
    size_text_theme = size_text_theme,
    size_text = size_text,
    width_legendbar = width_legendbar,
    height_legendbar = height_legendbar,
    vjust = vjust,
    max_chars = max_chars,
    axis_margin_top = axis_margin_top,
    axis_margin_right = axis_margin_right
  ))
}


# ============================================================================
# PLOT CREATION
# ============================================================================

#' Create odds ratio heatmap
create_odds_ratio_heatmap <- function(fixed_df, plot_params, P, ylab,
                                      color_lab = "Adjusted odds ratio") {

  # Build theme
  adapt_theme <- build_heatmap_theme(plot_params, P)

  # Create base plot with tiles
  pl_fixed <- fixed_df %>%
    ggplot() +
    geom_tile(
      aes(x = statistic, y = response, col = fill, fill = fill),
      alpha = 0.95,
      linewidth = 0.5
    )

  # Add text layers for significant estimates
  pl_fixed <- add_text_layers_significant(
    pl_fixed,
    fixed_df,
    plot_params,
    P
  )

  # Add text layers for non-significant estimates
  pl_fixed <- add_text_layers_nonsignificant(
    pl_fixed,
    fixed_df,
    plot_params,
    P
  )

  # Add scales and labels
  pl_fixed <- pl_fixed +
    scale_color_gradient(
      name = color_lab,
      low = "yellow",
      high = "red",
      na.value = "white"
    ) +
    scale_fill_gradient(
      name = color_lab,
      low = "yellow",
      high = "red",
      na.value = "white"
    ) +
    guides(
      col = NULL,
      fill = guide_colourbar(
        title.position = "top",
        theme = theme(
          legend.key.width = unit(plot_params$width_legendbar, "lines"),
          legend.key.height = unit(plot_params$height_legendbar, "lines")
        )
      )
    ) +
    labs(y = ylab, x = "Predictor") +
    adapt_theme +
    scale_x_discrete(
      labels = function(x) sapply(x, wrap_equally, max_chars = plot_params$max_chars),
      limits = rev,
      position = 'top',
      expand = expansion(add = c(0, 0))
    ) +
    scale_y_discrete(
      labels = function(x) sapply(x, wrap_equally, max_chars = plot_params$max_chars),
      expand = expansion(add = c(0, 0))
    )

  return(pl_fixed)
}


#' Build heatmap theme
build_heatmap_theme <- function(plot_params, P) {

  # Base theme from P
  base_size <- ifelse(
    plot_params$max_chars == 31,
    9.5,  # All events
    plot_params$size_text_theme  # Selected events
  )

  adapt_theme <- P$own_theme +
    theme(
      legend.position = 'bottom',
      axis.text.x = element_text(
        size = base_size,
        angle = 90,
        hjust = 0,
        vjust = 1
      ),
      axis.text.y = element_text(size = base_size),
      axis.title.x.top = element_text(
        size = plot_params$size_text_theme + 4,
        margin = margin(b = plot_params$axis_margin_top)
      ),
      axis.title.y = element_text(
        size = plot_params$size_text_theme + 4,
        margin = margin(r = plot_params$axis_margin_right)
      ),
      legend.text = element_text(
        size = ifelse(plot_params$max_chars == 31, base_size, base_size - 4)
      ),
      legend.title = element_text(
        size = ifelse(plot_params$max_chars == 31, base_size + 4, base_size)
      ),
      axis.title = element_text(size = base_size + 4),
      plot.margin = unit(c(.01, .15, .01, .01), 'cm'),
      panel.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.spacing = unit(0, "cm"),
      legend.margin = margin(0, 0, 0, 0)
    )

  return(adapt_theme)
}


#' Add text layers for significant estimates
add_text_layers_significant <- function(pl, fixed_df, plot_params, P) {

  sig_df <- fixed_df %>% dplyr::filter(excludes_1 == TRUE)

  pl +
    # Confidence intervals
    geom_text(
      data = sig_df,
      aes(x = statistic, y = response, label = ifelse(is.na(estimate), "", ci)),
      hjust = 0.5,
      vjust = plot_params$vjust,
      size = plot_params$size_text,
      family = P$font_family
    ) +
    # Point estimates (bold)
    geom_text(
      data = sig_df,
      aes(x = statistic, y = response, label = ifelse(is.na(estimate), "", est)),
      fontface = "bold",
      hjust = 0.5,
      vjust = plot_params$vjust,
      size = plot_params$size_text + 1.1,
      family = P$font_family
    )
}


#' Add text layers for non-significant estimates
add_text_layers_nonsignificant <- function(pl, fixed_df, plot_params, P) {

  nonsig_df <- fixed_df %>% dplyr::filter(excludes_1 == FALSE)

  pl +
    # Confidence intervals (grey)
    geom_text(
      data = nonsig_df,
      aes(x = statistic, y = response, label = ifelse(is.na(estimate), "", ci)),
      color = "grey70",
      hjust = 0.5,
      vjust = plot_params$vjust,
      size = plot_params$size_text,
      family = P$font_family
    ) +
    # Point estimates (bold grey)
    geom_text(
      data = nonsig_df,
      aes(x = statistic, y = response, label = ifelse(is.na(estimate), "", est)),
      fontface = "bold",
      color = "grey70",
      hjust = 0.5,
      vjust = plot_params$vjust,
      size = plot_params$size_text + 1.1,
      family = P$font_family
    )
}


# ============================================================================
# SAVING PLOTS
# ============================================================================

#' Save co-occurrence plot
save_cooccur_plot <- function(plot, datalist, name_dataset, chosen_lead,
                              suffix = "") {

  # Determine height
  height <- ifelse(name_dataset == "HILDA", 200, 200)

  # Build filename
  filename <-  sprintf("%s_heatmap_odds_ratio%s_lead%d.pdf",
                       name_dataset, suffix, chosen_lead)

  # Build filepath
  filepath_image <- file.path(datalist$filepath_figs_dataset, filename)

  # Save plot
  save_plot(plot, filepath_image, height = height)

  invisible(filepath_image)
}

# Create Supplementary Table
df_to_latex_cooccur = function(name_dataset, chosen_lead, model_df, datalist, P){

  latex_df = model_df %>% dplyr::filter(subplot != "fixed") %>%
    dplyr::mutate(estimate = ifelse(is.na(estimate), "",
                                    ifelse(subplot == "sd_intercept",
                                           sprintf("%.2f [%.2f, %.2f]", estimate, lower, upper),
                                           sprintf("%.2f", estimate)))) %>%
    dplyr::select(-c(lower, upper)) %>%
    recode_events(., "response", P$event_dicts[[name_dataset]], datalist$df_per_event) %>%
    tidyr::pivot_wider(names_from = c(subplot, statistic), values_from = estimate) %>%
    dplyr::mutate_at(setdiff(colnames(.), "response"), ~ ifelse(is.na(.), "", .)) %>%
    dplyr::select(response, sd_intercept_p_id_sd_intercept,
                  sd_intercept_crosswave_h_id_sd_intercept,
                  var_part_R2_marginal, var_part_R2_conditional,
                  var_part_ICC_adjusted)

  latex_df %>%
    kableExtra::kbl(format = "latex",
        col.names = NULL,
        align = 'lccccc',
        caption = sprintf("Model fit and estimates of %s associations between events (Source: %s). The household intercept had to be dropped for models with no indicated household %s. $i$ = Individual; $h$ = Household; Marg. $R^2$ = Marginal explained variance; Cond. $R^2$ = Conditional explained variance; Adj. ICC = Adjusted Intra-class Correlation Coefficient.",
                          ifelse(chosen_lead == 0, "contemporaneous", sprintf("lag-%d", chosen_lead)),
                          ifelse(name_dataset == "SHP", "Swiss Household Panel, SHP",
                                 "Household, Income and Labour Dynamics in Australia Survey, HILDA"),
                          # name_dataset,
                          "$\\sigma$"
        ),
        label = sprintf("%s_odds_ratio_lead%d", name_dataset, chosen_lead),
        booktabs = TRUE, escape = FALSE) %>%
    kableExtra::add_header_above(c(" ", "$i$" = 1, "$h$" = 1,
                                   "Marg. $R^2$", "Cond. $R^2$", "Adj. ICC"),
                                 escape = FALSE) %>%
    kableExtra::add_header_above(c(" ", "Random $\\\\sigma$" = 2,
                       "Variance Partitioning" = 3),
                     escape = FALSE) %>%
    kableExtra::column_spec(1, "1.5in") %>%
    kableExtra::column_spec(2:3, "1.1in") %>%
    kableExtra::column_spec(4:6, "0.2in")

}


#' Get valid predictors based on dataset rules
get_valid_predictors <- function(outcome, all_events, chosen_lead,
                                 name_dataset, event_dict) {

  predictors <- all_events

  if (name_dataset == "SHP") {
    # Drop "unspecified" as predictor to avoid collinearity (use as reference category)
    predictors <- setdiff(predictors, "unspecified_illness_accident")

    # For contemporaneous models, illness/accident events can't co-occur
    if (chosen_lead == 0) {
      illness_events <- event_dict[grepl("PL01R", event_dict[["var_name"]]),
                                   "recode_var"]

      if (outcome %in% illness_events) {
        predictors <- setdiff(predictors, illness_events)
      }
    }
  }

  # Remove outcome from predictors (except for lag models where we include it)
  if (chosen_lead == 0) {
    predictors <- setdiff(predictors, outcome)
  }

  return(predictors)
}


#' Prepare data for modeling (apply lag if needed)
prepare_model_data <- function(df_events, outcome, chosen_lead) {

  if (chosen_lead > 0) {
    new_name <- paste0("lead_", outcome)

    df_model <- df_events %>%
      group_by(p_id) %>%
      mutate(!!new_name := lead(!!sym(outcome), chosen_lead)) %>%
      ungroup()
  } else {
    df_model <- df_events
  }

  return(df_model)
}


#' Build model formulas for all three specifications
build_model_formulas <- function(outcome, predictors, chosen_lead) {

  # Determine outcome variable name
  outcome_var <- if (chosen_lead > 0) {
    paste0("lead_", outcome)
  } else {
    outcome
  }

  # Base formula
  formula_base <- paste0(
    outcome_var, " ~ ",
    paste(predictors, collapse = " + ")
  )

  # Add random effects
  formulas <- list(
    no_random = formula_base,
    p_id = paste0(formula_base, " + (1 | p_id)"),
    p_id_h_id = paste0(formula_base, " + (1 | p_id) + (1 | crosswave_h_id)")
  )

  return(formulas)
}


#' Fit all three model variants
fit_model_variants <- function(formulas, df_model) {

  model = glmmTMB(
    as.formula(formulas$p_id_h_id),
    data = df_model,
    family = binomial(link = "logit")
  )

  if (performance::check_singularity(model) | !performance::check_convergence(model)){
    model = glmmTMB(
      as.formula(formulas$p_id),
      data = df_model,
      family = binomial(link = "logit")
    )
  }

  if (performance::check_singularity(model) | !performance::check_convergence(model)){
    model = glmmTMB(
      as.formula(formulas$no_random),
      data = df_model,
      family = binomial(link = "logit")
    )
  }

  return(model)
}


# ============================================================================
# UNADJUSTED SENSITIVITY ANALYSIS
# ============================================================================

#' Run unadjusted models
run_unadjusted_sensitivity <- function(model_list, negevent_names, chosen_lead,
                                      name_dataset, datalist, P
                                      ) {

  # Fit unadjusted models
  filepath_unadjusted <- file.path(
    datalist$filepath_deriv,
    sprintf("cooccurrence_event_types_unadjusted_lead%d.RDS", chosen_lead)
  )

  cat("Fitting unadjusted models...\n")

  if (!file.exists(filepath_unadjusted) || rerun) {
    start_t = Sys.time()

    unadjusted_models <- fit_unadjusted_models(model_list, negevent_names,
                                             chosen_lead)
    saveRDS(unadjusted_models, filepath_unadjusted)

    cat(sprintf("Time elapsed: %s\n", format(Sys.time() - start_t)))
  }
  unadjusted_models <- readRDS(filepath_unadjusted)

  # Plot difference
  unadjusted_df <- lapply(unadjusted_models, function(x){
    data.frame(subplot = "fixed",
      response = x[["outcome"]],
      statistic = x[["predictor"]],
      estimate = x[["unadjusted_est"]],
      lower = x[["unadjusted_lower"]],
      upper = x[["unadjusted_upper"]],
      full_estimate = x[["full_est"]],
      full_lower = x[["full_lower"]],
      full_upper = x[["full_upper"]]
               )
  }) %>% do.call(rbind, .) %>% as.data.frame()

  plot_cooccur(name_dataset, datalist, P, unadjusted_df, chosen_lead,
               color_lab = "Unadjusted odds ratio",
               suffix = "_unadjusted")

  # Get odds ratios
  comp_df <- unadjusted_df %>%
    # Exponentiate to obtain odds ratio
    dplyr::mutate_at(c("estimate", "lower", "upper",
                     "full_estimate", "full_lower", "full_upper"),
                   ~ round(exp(.), 2)) %>%
    # Determine statistical significance
    dplyr::mutate(
      excludes_1 = (lower < 1 & upper < 1) | (lower > 1 & upper > 1),
      full_excludes_1 = (full_lower < 1 & full_upper < 1) | (full_lower > 1 & full_upper > 1)
    )

  # head(comp_df)

  # Which event combinations show different results for unadjusted vs. adjusted?
  cat("\nEvent combinations where the unadjusted model shows a significant relationship, but the adjusted model does not:\n")
  print(comp_df %>% dplyr::filter(excludes_1, !full_excludes_1) )

  cat("\nEvent combinations where the unadjusted model does not show a significant relationship, but the adjusted model does:\n")
  print(comp_df %>% dplyr::filter(!excludes_1, full_excludes_1) )

  # Different conclusions
  cat("\nEvent combinations where the unadjusted and adjusted model disagree on whether the association is positive (OR > 1) or negative (OR < 1), for models which show a significant difference in either the unadjusted and/or the adjusted model:\n")
  print(comp_df %>%
      # dplyr::filter(excludes_1 != full_excludes_1) %>%
      dplyr::filter(excludes_1 | full_excludes_1) %>%
      dplyr::filter(estimate > 1 & full_estimate < 1 | estimate < 1 & full_estimate > 1))

  # Different conclusions
  cat("\nEvent combinations where the adjusted model is significant, and the unadjusted and adjusted model disagree on whether the association is positive (OR > 1) or negative (OR < 1):\n")
  print(comp_df %>%
          dplyr::filter(full_excludes_1) %>%
          dplyr::filter(estimate > 1 & full_estimate < 1 | estimate < 1 & full_estimate > 1))

  cat("\nChange in size of association:\n")
  x <- comp_df %>%
          dplyr::mutate(pct_change = 100 * (estimate - full_estimate) / full_estimate) %>%
    dplyr::pull(pct_change)
  x <- x[!is.na(x)]
  print(sprintf("\nMean: %.2f; Median: %.2f; Variance: %.2f; SD: %.2f",
                mean(x), median(x), var(x), sd(x)))



  # unadjusted_df <- data.frame(
  #     outcome = unlist(lapply(unadjusted_models, `[[`, "outcome")),
  #     predictor = unlist(lapply(unadjusted_models, `[[`, "predictor")),
  #     full_OR = exp(unlist(lapply(unadjusted_models, `[[`, "full_est"))),
  #     unadjusted_OR = exp(unlist(lapply(unadjusted_models, `[[`, "unadjusted_est")))
  #   )
  # unadjusted_df[["unadjusted_minus_full_OR"]] = unadjusted_df[["unadjusted_OR"]] - unadjusted_df[["full_OR"]]
  #
  # # Check if any OR are different in direction
  # check <- unadjusted_df
  # check2 <- check[check[["unadjusted_OR"]] > 1 & check[["full_OR"]] < 1 | check[["unadjusted_OR"]] < 1 & check[["full_OR"]] > 1, ]
  # print(check2)


  # plot_cooccur_unadjusted(name_dataset, datalist, P, unadjusted_df,
  #                                     chosen_lead)

  return(NULL)
}



#' Fit unadjusted models using same sample as full models
fit_unadjusted_models <- function(model_list, negevent_names, chosen_lead) {

  n <- length(model_list)

  unadjusted_models <- foreach::foreach(
    i = rep(1:n, each = n), # outcome
    j = rep(1:n, n), # predictor
    .packages = c("glmmTMB", "dplyr", "purrr", "broom.mixed"),
    .export = c("model_list", "negevent_names", "chosen_lead")
  ) %dopar% {

    outcome <- model_list[[i]]$outcome
    predictor <- negevent_names[j]

    cat("\n")
    cat(sprintf("  Outcome: %s\n", outcome))
    cat(sprintf("  Predictor: %s\n", predictor))

    # Don't run
    if (i == j & chosen_lead == 0){
      cat("  Skip\n")
      return(list(
        outcome = outcome,
        predictor = predictor,
        model_df = NA,
        full_est = NA,
        unadjusted_est = NA,
        full_lower = NA,
        full_upper = NA,
        unadjusted_lower = NA,
        unadjusted_upper = NA
      ))
    }

    # Get dataframe with same observations
    df <- model_list[[i]]$model$frame

    # If the predictor was not a predictor in the original model, stop
    if (!predictor %in% colnames(df)){
      cat("  Skip\n")
      return(list(
        outcome = outcome,
        predictor = predictor,
        model_df = NA,
        full_est = NA,
        unadjusted_est = NA,
        full_lower = NA,
        full_upper = NA,
        unadjusted_lower = NA,
        unadjusted_upper = NA
      ))
    }

    # Find which random effects are needed
    nr_obs <- model_list[[i]]$nr_obs
    incl_p_id <- "p_id" %in% names(nr_obs)
    incl_h_id <- "crosswave_h_id" %in% names(nr_obs)

    # Build unadjusted formula
    outcome_var <- if (chosen_lead > 0) paste0("lead_", outcome) else outcome
    formula <- paste0(outcome_var, " ~ ", predictor,
                      ifelse(incl_p_id, " + (1 | p_id)", ""),
                      ifelse(incl_h_id, " + (1 | crosswave_h_id)", ""))

    # Fit model
    model <- glmmTMB(
        as.formula(formula),
        data = df,
        family = binomial(link = "logit")
      )

    # Ensure the number of observations is the same
    stopifnot(identical(purrr::map_vec(ranef(model)$cond, nrow), nr_obs))

    # # Summarise
    # model_summ <- summarise_model(model)
    # model_df <- get_plot_df(model_summ, outcome)

    # Extract model estimates
    model_df <- broom.mixed::tidy(model, effects = "fixed")

    ci_df <- as.data.frame(confint(model)) %>%
      tibble::rownames_to_column("term")
    unadjusted_lower <- ci_df[ci_df[["term"]] == predictor, ][["2.5 %"]]
    unadjusted_upper <- ci_df[ci_df[["term"]] == predictor, ][["97.5 %"]]


    # Get difference in estimate
    model_df_orig <- model_list[[i]]$model_df
    full_est <- model_df_orig[model_df_orig[["statistic"]] == predictor, ][["estimate"]]
    full_lower <- model_df_orig[model_df_orig[["statistic"]] == predictor, ][["lower"]]
    full_upper <- model_df_orig[model_df_orig[["statistic"]] == predictor, ][["upper"]]
    unadjusted_est <- model_df[model_df[["term"]] == predictor, ][["estimate"]]

    cat(sprintf("  Odds Ratio Difference: %.4f\n", exp(full_est) - exp(unadjusted_est)))

    return(list(
      outcome = outcome,
      predictor = predictor,
      model_df = model_df,
      full_est = full_est,
      unadjusted_est = unadjusted_est,
      full_lower = full_lower,
      full_upper = full_upper,
      unadjusted_lower = unadjusted_lower,
      unadjusted_upper = unadjusted_upper
    ))
  }

  return(unadjusted_models)
}




##### JOINT AND CONDITIONAL PROBABILITY ANALYSIS - REFACTORED #####

#' Main wrapper function for joint/conditional probability analysis
#' @param name_dataset Dataset name (e.g., "SHP", "HILDA")
#' @param filepath_base Base file path
#' @param P Parameters list
#' @param chosen_lags Lag values to analyze (e.g., c(0, 1))
#' @param rerun Whether to rerun analysis
run_joint_cond_prob <- function(name_dataset, filepath_base, P,
                                chosen_lags = c(0, 1)) {

  # 1. Prepare data
  cat("\n=== Preparing data ===\n")
  event_dict <- P$event_dicts[[name_dataset]]
  datalist <- prepare_data(name_dataset, filepath_base, event_dict)

  # 2. Format event data
  cat("\n=== Formatting event data ===\n")
  df_events <- format_event_data(datalist, keep_h_id = FALSE)
  negevent_names <- get_event_names(df_events)

  # 3. Generate event combinations
  cat("\n=== Generating event combinations ===\n")
  combo <- generate_event_combinations(negevent_names)

  # 4. Calculate probabilities for each lag
  cat("\n=== Calculating probabilities ===\n")
  prob_df <- calculate_probabilities_for_lags(
    df_events = df_events,
    combo = combo,
    chosen_lags = chosen_lags
  )

  # 5. Return consolidated results
  return(list(
    datalist = datalist,
    df_events = df_events,
    prob_df = prob_df
  ))
}


# ============================================================================
# DATA PREPARATION
# ============================================================================

#' Generate all event combinations
generate_event_combinations <- function(negevent_names) {
  combo <- expand.grid(
    eventA = negevent_names,
    eventB = negevent_names,
    stringsAsFactors = FALSE
  ) %>%
    mutate_at(c("eventA", "eventB"), ~ as.character(.)) %>%
    split(seq(nrow(.)))

  return(combo)
}


# ============================================================================
# PROBABILITY CALCULATIONS
# ============================================================================

#' Calculate probabilities for all specified lags
calculate_probabilities_for_lags <- function(df_events, combo, chosen_lags) {

  start_t <- Sys.time()
  prob_df <- lapply(chosen_lags, function(chosen_lag) {
    cat(sprintf("  Calculating for lag-%d\n", chosen_lag))

    cbind(
      lag_eventA = chosen_lag,
      combo %>% do.call(rbind, .),
      cond_prob_given_eventB = get_prob(
        df_events, combo, chosen_lag,
        lag_event = "A",
        type = "cond_prob_given_eventB"
      ),
      joint_prob = get_prob(
        df_events, combo, chosen_lag,
        lag_event = "A",
        type = "joint_prob"
      )
    )
  }) %>%
    do.call(rbind, .) %>%
    as.data.frame()
  cat(sprintf("Time elapsed: %s\n", format(Sys.time() - start_t)))

  return(prob_df)
}

#' Calculate joint or conditional probabilities
#' @param df_events Event data frame
#' @param combo List of event combinations
#' @param chosen_lag Lag value (0 for contemporaneous)
#' @param lag_event Which event to lag ("A" or "B")
#' @param type Type of probability to calculate
get_prob <- function(df_events, combo,
                     chosen_lag = 0,
                     lag_event = c("A", "B")[1],
                     type = c("joint_prob", "cond_prob_given_eventA",
                              "cond_prob_given_eventB")[1]) {

  cond <- foreach(
    x = combo,
    .packages = c("dplyr"),
    .combine = "c"
  ) %do% {

    # Select relevant columns
    df_ <- df_events %>%
      select_at(c("p_id", unname(unlist(x))))

    # Apply lag if needed
    if (chosen_lag > 0) {
      df_ <- apply_lag_to_events(df_, x, chosen_lag, lag_event)
    }

    # Remove NA and calculate probability
    df_ <- df_ %>%
      dplyr::select(-p_id) %>%
      dplyr::slice(which(complete.cases(.)))

    prob <- calculate_probability_metric(df_, x, chosen_lag, type)

    rm(df_)
    return(prob)
  }

  return(cond)
}

#' Apply lag to events
apply_lag_to_events <- function(df_, x, chosen_lag, lag_event) {
  if (lag_event == "A") {
    df_[["orig"]] <- df_[[x$eventA]]
    df_[[x$eventA]] <- unsplit(
      lapply(split(df_[[x$eventA]], df_$p_id),
             function(y) c(rep(NA, chosen_lag), head(y, -chosen_lag))),
      df_$p_id
    )
  } else if (lag_event == "B") {
    df_[["orig"]] <- df_[[x$eventB]]
    df_[[x$eventB]] <- unsplit(
      lapply(split(df_[[x$eventB]], df_$p_id),
             function(y) c(rep(NA, chosen_lag), head(y, -chosen_lag))),
      df_$p_id
    )
  }

  return(df_)
}

#' Calculate specific probability metric
calculate_probability_metric <- function(df_, x, chosen_lag, type) {

  # Handle same event combinations
  if (x$eventA == x$eventB) {
    if (chosen_lag == 0) {
      return(mean(df_[[x$eventA]]))
    } else {
      return(sum(rowSums(df_) == 2) / nrow(df_))
    }
  }

  # Remove original unlagged event
  if (chosen_lag > 0) {
    df_$orig <- NULL
  }

  # Calculate requested probability type
  if (type == "joint_prob") {
    prob <- sum(rowSums(df_) == 2) / nrow(df_)
  } else if (type == "cond_prob_given_eventA") {
    prob <- mean(rowSums(df_) == 2) / mean(df_[[x$eventA]])
  } else if (type == "cond_prob_given_eventB") {
    prob <- mean(rowSums(df_) == 2) / mean(df_[[x$eventB]])
  }

  return(prob)
}


# ============================================================================
# VISUALIZATION
# ============================================================================

#' Plot joint or conditional probability heatmap
plot_joint_cond <- function(name_dataset, datalist, prob_df, chosen_lag,
                            event_dict, P,
                            lag_event = c("A", "B")[1],
                            type = c("joint_prob", "cond_prob_given_eventA",
                                     "cond_prob_given_eventB")[1],
                            size_text = 9,
                            size_label = 2.8,
                            size_point = 8) {

  # 1. Prepare plot data
  plot_df <- prepare_plot_data(
    prob_df = prob_df,
    chosen_lag = chosen_lag,
    lag_event = lag_event,
    type = type,
    name_dataset = name_dataset,
    event_dict = event_dict
  )

  # 2. Get plot labels
  labels <- get_plot_labels(chosen_lag, lag_event, type)

  # 3. Recode events with descriptions
  plot_df <- recode_events(
    plot_df,
    c("eventA", "eventB"),
    P$event_dicts[[name_dataset]],
    datalist$df_per_event
  )

  # 4. Create plot
  pl_prob <- create_probability_plot(
    plot_df = plot_df,
    labels = labels,
    P = P,
    size_text = size_text,
    size_label = size_label,
    size_point = size_point
  )

  # 5. Save plot
  height <- ifelse(name_dataset == "HILDA", 190, 190)
  filepath_image <- file.path(
    datalist$filepath_figs_dataset,
    sprintf("%s_%s_lag%d.pdf", name_dataset, type, chosen_lag)
  )
  save_plot(pl_prob, filepath_image, height = height)

  return(pl_prob)
}

#' Prepare data for plotting
prepare_plot_data <- function(prob_df, chosen_lag, lag_event, type,
                              name_dataset, event_dict) {

  plot_df <- prob_df %>%
    filter(!!sym(sprintf("lag_event%s", lag_event)) == !!chosen_lag) %>%
    mutate(estimate = !!sym(type))

  # Set estimate to NA for events that cannot co-occur
  if (name_dataset == "SHP" && chosen_lag == 0) {
    illness_accident_events <- event_dict %>%
      filter(grepl("PL01R", var_name)) %>%
      pull(recode_var)

    plot_df <- plot_df %>%
      mutate(estimate = ifelse(
        eventA %in% illness_accident_events &
          eventB %in% illness_accident_events &
          eventA != eventB,
        NA,
        estimate
      ))
  }

  # Add label column based on type
  if (type == "joint_prob") {
    plot_df$label <- log10(plot_df$estimate)
  } else {
    plot_df$label <- plot_df$estimate
  }

  return(plot_df)
}

#' Get plot axis labels
get_plot_labels <- function(chosen_lag, lag_event, type) {

  # Determine axis labels
  if (chosen_lag == 0) {
    xlab <- "Event B"
    ylab <- "Event A"
    x <- "eventB"
    y <- "eventA"
  } else {
    lag_text <- sprintf("%d year%s later", chosen_lag,
                        ifelse(chosen_lag > 1, "s", ""))

    if (lag_event == "A") {
      ylab <- sprintf("Event A (%s)", lag_text)
      xlab <- "Event B"
      x <- "eventB"
      y <- "eventA"
    } else {
      ylab <- sprintf("Event B (%s)", lag_text)
      xlab <- "Event A"
      x <- "eventA"
      y <- "eventB"
    }
  }

  # Determine fill label
  if (type == "joint_prob") {
    fill_name <- "Joint probability"
  } else if (type == "cond_prob_given_eventA") {
    fill_name <- "Conditional probability of event B given event A"
  } else {
    fill_name <- "Conditional probability of event A given event B"
  }

  return(list(
    xlab = xlab,
    ylab = ylab,
    x = x,
    y = y,
    fill_name = fill_name
  ))
}

#' Create the probability plot
create_probability_plot <- function(plot_df, labels, P,
                                    size_text, size_label, size_point) {

  # Set theme
  adapt_theme <- P$own_theme +
    theme(
      legend.position = 'bottom',
      axis.text.x = element_text(size = size_text + 2, angle = 90,
                                 hjust = 0, vjust = 1),
      axis.text.y = element_text(size = size_text + 2),
      legend.text = element_text(size = size_text),
      legend.title = element_text(size = size_text + 4),
      axis.title = element_text(size = size_text + 4),
      plot.margin = unit(c(.01, .15, .01, .01), 'cm'),
      panel.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.spacing = unit(0, "cm"),
      legend.margin = margin(0, 0, 0, 0)
    )

  width_legendbar <- 10
  height_legendbar <- 1

  # Create plot
  pl_prob <- plot_df %>%
    ggplot() +
    geom_point(
      aes(x = .data[[labels$x]], y = .data[[labels$y]], col = estimate),
      shape = 19,
      size = size_point
    ) +
    geom_text(
      aes(x = .data[[labels$x]], y = .data[[labels$y]],
          label = ifelse(is.na(estimate), "", sprintf("%.3f", as.numeric(estimate)))),
      size = size_label,
      family = P$font_family
    ) +
    scale_color_gradient(
      name = labels$fill_name,
      low = "yellow",
      high = "red",
      na.value = "white"
    ) +
    guides(color = guide_colourbar(
      title.position = "top",
      theme = theme(
        legend.key.width = unit(width_legendbar, "lines"),
        legend.key.height = unit(height_legendbar, "lines")
      )
    )) +
    labs(y = labels$ylab, x = labels$xlab) +
    adapt_theme +
    scale_x_discrete(
      labels = function(x) sapply(x, wrap_equally, max_chars = 31),
      limits = rev,
      position = 'top',
      expand = expansion(add = c(0.6, 0.6))
    ) +
    scale_y_discrete(
      labels = function(x) sapply(x, wrap_equally, max_chars = 31),
      expand = expansion(add = c(0.6, 0.6))
    )

  return(pl_prob)
}
