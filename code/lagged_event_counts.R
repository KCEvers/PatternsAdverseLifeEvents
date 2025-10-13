##### LAGGED EVENT COUNTS ANALYSIS - REFACTORED #####

#' Main wrapper function for lagged event counts analysis
#' @param name_dataset Dataset name (e.g., "SHP", "HILDA")
#' @param filepath_base Base file path
#' @param P Parameters list
#' @param rerun Whether to rerun analysis
run_lagged <- function(name_dataset, filepath_base, P, rerun = FALSE) {

  # 1. Prepare data
  cat("\n=== Preparing data ===\n")
  event_dict <- P$event_dicts[[name_dataset]]
  datalist <- prepare_data(name_dataset, filepath_base, event_dict)

  # 2. Create lagged dataset
  cat("\n=== Creating lagged dataset ===\n")
  df_lag <- create_lagged_dataset(datalist)

  # 3. Define file path for saved models
  filepath_models <- file.path(
    datalist$filepath_deriv,
    "lagged_event_counts.RDS"
  )

  # 4. Fit models or load existing
  if (!file.exists(filepath_models) || rerun) {
    cat("\n=== Fitting models ===\n")
    results <- fit_lagged_models(df_lag, filepath_models)
  } else {
    cat("\n=== Loading saved models ===\n")
    results <- readRDS(filepath_models)
  }

  # 5. Return consolidated results
  return(results)
}


# ============================================================================
# DATA PREPARATION
# ============================================================================

#' Create dataset with lagged event counts
#' @param datalist Data list from prepare_data
create_lagged_dataset <- function(datalist) {

  df_lag <- datalist$df_nr_negevents_pp_py %>%
    dplyr::select(p_id, crosswave_h_id, age, wave_nr, nr_occur) %>%
    # Ensure each participant has a row for each wave (for correct lagging)
    tidyr::complete(p_id, wave_nr) %>%
    dplyr::arrange(p_id, wave_nr) %>%
    # Add lag-1 column to event counts
    dplyr::mutate(lag_nr_occur = lag(nr_occur), .by = p_id) %>%
    # Remove missing values created by lagging
    dplyr::filter(!is.na(lag_nr_occur), !is.na(nr_occur)) %>%
    # Scale age predictor (don't scale lag to interpret one-unit increase)
    dplyr::mutate(age_scaled = scale(age) %>% as.numeric()) %>%
    # Convert IDs to factors
    dplyr::mutate(
      p_id = as.factor(p_id),
      crosswave_h_id = as.factor(crosswave_h_id)
    )

  # Report dataset statistics
  cat(sprintf("  Number of unique participants: %d\n", length(unique(df_lag$p_id))))
  cat(sprintf("  Number of observations: %d\n", nrow(df_lag)))

  return(df_lag)
}


# ============================================================================
# MODEL FITTING
# ============================================================================

#' Fit lagged event count models
#' @param df_lag Lagged dataset
#' @param filepath_models Path to save models
fit_lagged_models <- function(df_lag, filepath_models) {

  # Define model formulas
  formulas <- build_lagged_formulas()

  # 1. Fit random intercept model
  cat("  Fitting random intercept model...\n")
  model_random_int <- glmmTMB::glmmTMB(
    as.formula(formulas$random_int),
    data = df_lag,
    family = poisson(link = "log")
  )

  # 2. Fit random intercept and slope model
  cat("  Fitting random intercept + slope model...\n")
  model_random_int_slope <- glmmTMB::glmmTMB(
    as.formula(formulas$random_int_slope),
    data = df_lag,
    family = poisson(link = "log")
  )

  # 3. Check convergence
  cat("  Checking convergence...\n")
  df_convergence <- check_model_convergence(model_random_int, model_random_int_slope)

  # 4. Compare model fit
  cat("  Comparing model fit...\n")
  fit_comparison <- compare_lagged_model_fit(model_random_int, model_random_int_slope)

  # 5. Compile results
  results <- list(
    df_lag = df_lag,
    model_random_int = model_random_int,
    model_random_int_slope = model_random_int_slope,
    df_convergence = df_convergence,
    df_ABIC = fit_comparison$df_ABIC,
    df_comp_AIC = fit_comparison$df_comp_AIC,
    df_comp_BIC = fit_comparison$df_comp_BIC
  )

  # 6. Save results
  saveRDS(results, filepath_models)

  return(results)
}


#' Build formulas for lagged models
build_lagged_formulas <- function() {

  # Base formula with fixed effects
  base_formula <- "nr_occur ~ 1 + age_scaled + I(age_scaled^2) + lag_nr_occur"

  # Random intercept only (partially crossed random effects)
  formula_random_int <- paste0(
    base_formula,
    " + (1 | p_id) + (1 | crosswave_h_id)"
  )

  # Random intercept and slope for lag
  formula_random_int_slope <- paste0(
    base_formula,
    " + (1 + lag_nr_occur | p_id) + (1 | crosswave_h_id)"
  )

  return(list(
    random_int = formula_random_int,
    random_int_slope = formula_random_int_slope
  ))
}


# ============================================================================
# MODEL DIAGNOSTICS
# ============================================================================

#' Check convergence of fitted models
#' @param model_random_int Random intercept model
#' @param model_random_int_slope Random intercept + slope model
check_model_convergence <- function(model_random_int, model_random_int_slope) {

  df_convergence <- data.frame(
    model_random_int = model_random_int$fit$convergence == 0,
    model_random_int_slope = model_random_int_slope$fit$convergence == 0
  )

  # Report convergence status
  if (all(df_convergence)) {
    cat("  ✓ All models converged successfully\n")
  } else {
    cat("  ⚠ Warning: Some models did not converge\n")
    print(df_convergence)
  }

  return(df_convergence)
}


#' Compare fit of lagged models using AIC and BIC
#' @param model_random_int Random intercept model
#' @param model_random_int_slope Random intercept + slope model
compare_lagged_model_fit <- function(model_random_int, model_random_int_slope) {

  # Extract AIC and BIC
  df_ABIC <- cbind(
    AIC(model_random_int, model_random_int_slope),
    BIC(model_random_int, model_random_int_slope)
  )

  # Create pairwise comparisons
  model_names <- c("model_random_int", "model_random_int_slope")

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

  # Report which model fits better
  delta_AIC <- df_comp_AIC[1, 2]
  if (abs(delta_AIC) > 2) {
    better_model <- ifelse(delta_AIC > 0, "Random intercept + slope", "Random intercept")
    cat(sprintf("  Better fit: %s (ΔAIC = %.2f)\n", better_model, abs(delta_AIC)))
  } else {
    cat("  Models have similar fit (ΔAIC < 2)\n")
  }

  return(list(
    df_ABIC = df_ABIC,
    df_comp_AIC = df_comp_AIC,
    df_comp_BIC = df_comp_BIC
  ))
}


# ============================================================================
# FORMATTING AND TABLES
# ============================================================================

#' Create dataframe with lagged analysis estimates for one dataset
#' @param name_dataset Dataset name
#' @param df_est Dataframe with model estimates
create_lagged_estimates_df <- function(name_dataset, df_est) {

  # Clean up estimate names
  df <- df_est %>%
    dplyr::rename(lower = "2.5 %", upper = "97.5 %") %>%
    dplyr::mutate(
      term = stringr::str_replace_all(term, "\\(Intercept\\)", "intercept")
    )

  # Extract and format fixed effects (exponentiate to get rate ratios)
  df_fixed <- df %>%
    dplyr::filter(effect == "fixed") %>%
    dplyr::mutate_at(c("estimate", "lower", "upper"), exp) %>%
    dplyr::mutate(label = sprintf("%.2f [%.2f, %.2f]", estimate, lower, upper)) %>%
    dplyr::select(term, label) %>%
    dplyr::pull(label, term)

  # Extract and format random effects (keep as SD)
  df_random <- df %>%
    dplyr::filter(effect == "ran_pars") %>%
    dplyr::mutate(label = sprintf("%.2f [%.2f, %.2f]", estimate, lower, upper)) %>%
    dplyr::pull(label, group)

  # Compile into single row
  data.frame(
    dataset = name_dataset,
    intercept = df_fixed[["intercept"]],
    lag_nr_occur = df_fixed[["lag_nr_occur"]],
    age = df_fixed[["age_scaled"]],
    age_squared = df_fixed[["I(age_scaled^2)"]],
    sd_intercept_p_id = df_random[["p_id"]],
    sd_intercept_h_id = df_random[["crosswave_h_id"]]
  )
}


#' Create LaTeX table for lagged analysis
#' @param summlist_SHP Summary list for SHP dataset
#' @param summlist_HILDA Summary list for HILDA dataset
create_lagged_table <- function(summlist_SHP, summlist_HILDA) {

  rbind(
    create_lagged_estimates_df("SHP", summlist_SHP$df_est),
    create_lagged_estimates_df("HILDA", summlist_HILDA$df_est)
  ) %>%
    kableExtra::kbl(
      format = "latex",
      escape = FALSE,
      caption = "Model estimates of the autocorrelation in event counts (profile 95\\% profile confidence intervals in brackets). Estimates are on the natural scale. $\\lambda$ = Yearly rate of adverse life events; $i$ = Individual; $h$ = Household (Source: Swiss Household Panel, SHP and Household, Income and Labour Dynamics in Australia Survey, HILDA).",
      booktabs = TRUE,
      label = "autocorrelation",
      col.names = NULL
    ) %>%
    kableExtra::add_header_above(
      c(" ", "$\\\\lambda$" = 1, "Lag-1 counts" = 1, "Age" = 1,
        "Age$^2$" = 1, "$i$" = 1, "$h$" = 1),
      escape = FALSE
    ) %>%
    kableExtra::add_header_above(
      c(" ", "Fixed" = 4, "Random $\\\\sigma$" = 2),
      escape = FALSE
    ) %>%
    kableExtra::kable_styling(latex_options = "scale_down")
}
