##### SUPPORTING FUNCTIONS FOR LAGGED EVENT COUNTS ANALYSIS #####
run_lagged = function(name_dataset, filepath_base, P, rerun = F){
  
  # Load data
  datalist = prepare_data(name_dataset, filepath_base, P$event_dicts, rerun = rerun)
  
  ## Prepare dataframe
  df_lag = datalist$df_nr_negevents_pp_py %>%
    dplyr::select(p_id, crosswave_h_id, age, wave_nr, nr_occur) %>%
    # Ensure each participant has a row for each wave so that the lag is computed correctly
    tidyr::complete(p_id, wave_nr) %>%
    arrange(p_id, wave_nr) %>%
    # Add lag column to event counts
    mutate(lag_nr_occur = lag(nr_occur), .by = p_id) %>%
    # Remove missing values because of lagging
    dplyr::filter(!is.na(lag_nr_occur), !is.na(nr_occur)) %>%
    # Don't scale lag predictor as we want to interpret the effect of a one-unit increase in lag
    # Scale age predictor
    mutate(age_scaled = scale(age) %>% as.numeric()) %>%
    mutate(p_id = as.factor(p_id), 
           crosswave_h_id = as.factor(crosswave_h_id))
  head(df_lag)
  length(unique(df_lag$p_id))
  nrow(df_lag)
  
  filepath = file.path(datalist$filepath_deriv, sprintf("lagged_event_counts.RDS"))
  if (!file.exists(filepath)){
    
    # Formulas (partially crossed random effects)
    formula_random_int = "nr_occur ~ 1 + age_scaled + I(age_scaled^2) + lag_nr_occur + (1 | p_id) + (1 | crosswave_h_id)"
    formula_random_int_slope = "nr_occur ~ 1 + age_scaled + I(age_scaled^2) + lag_nr_occur + (1 + lag_nr_occur| p_id) + (1 | crosswave_h_id)" 
    
    ## Fit mixed-effects models
    model_random_int <- glmmTMB(
      as.formula(formula_random_int),
      data = df_lag, 
      family = poisson(link = "log"))
    
    model_random_int_slope <- glmmTMB(
      as.formula(formula_random_int_slope),
      data = df_lag, 
      family = poisson(link = "log"))
    
    # Model convergence
    df_convergence = data.frame(
      model_random_int = !model_random_int$fit$convergence,
      model_random_int_slope = !model_random_int_slope$fit$convergence
    )
    
    # Compare model fit
    df_ABIC = cbind(AIC(model_random_int, model_random_int_slope),
                    BIC(model_random_int, model_random_int_slope))
    
    model_names <- c("model_random_int", "model_random_int_slope")  # Label models
    df_comp_AIC = outer(df_ABIC$AIC %>% stats::setNames(model_names), 
                        df_ABIC$AIC %>% stats::setNames(model_names), FUN = "-")
    df_comp_BIC = outer(df_ABIC$BIC %>% stats::setNames(model_names), 
                        df_ABIC$BIC %>% stats::setNames(model_names), FUN = "-")
    
    out = list(df_lag = df_lag,
               model_random_int = model_random_int, 
               model_random_int_slope = model_random_int_slope,
               df_convergence = df_convergence,
               df_ABIC = df_ABIC, df_comp_AIC = df_comp_AIC, df_comp_BIC = df_comp_BIC)
    saveRDS(out, filepath)
  }
  out = readRDS(filepath)
  
  return(out)
  
}

# Create dataframe with estimates of lagged analysis
lagged_create_df = function(name_dataset, df_est){
  
  df = df_est %>% 
    dplyr::rename(lower = "2.5 %", upper = "97.5 %") %>%
    dplyr::mutate(term = stringr::str_replace_all(term, "\\(Intercept\\)", "intercept"))
  
  df_fixed = df %>% filter(effect == "fixed") %>%
    dplyr::mutate_at(c("estimate", "lower", "upper"), exp)  %>% 
    mutate(label = sprintf("%.2f [%.2f, %.2f]", estimate, lower, upper)) %>%
    dplyr::select(term, label) %>%
    pull(label, term)
  
  df_random = df %>% filter(effect == "ran_pars")  %>% 
    mutate(label = sprintf("%.2f [%.2f, %.2f]", estimate, lower, upper)) %>%
    pull(label, group)
  
  data.frame(dataset = name_dataset,
             intercept = df_fixed[["intercept"]],
             lag_nr_occur = df_fixed[["lag_nr_occur"]],
             age = df_fixed[["age_scaled"]],
             age_squared = df_fixed[["I(age_scaled^2)"]],
             sd_intercept_p_id = df_random[["p_id"]],
             sd_intercept_h_id = df_random[["crosswave_h_id"]]
  )
  
}

# Create table combining SHP and HILDA with estimates of lagged analysis
lagged_table = function(summlist_SHP, summlist_HILDA){
  rbind(lagged_create_df("SHP", summlist_SHP$df_est), lagged_create_df("HILDA", summlist_HILDA$df_est)) %>% 
    kbl(format = "latex",
        escape=F,
        caption = "Model estimates of the autocorrelation in event counts (profile 95\\% profile confidence intervals in brackets). Estimates are on the natural scale. $\\lambda$ = Yearly rate of adverse life events; $i$ = Individual; $h$ = Household (Source: Swiss Household Panel, SHP and Household, Income and Labour Dynamics in Australia Survey, HILDA).",
        # label = "tab:autocorrelation",
        booktabs = TRUE,
        label = "autocorrelation",
        col.names = NULL) %>%
    add_header_above(c(" ", "$\\\\lambda$" = 1, "Lag-1 counts" = 1, "Age" = 1, "Age$^2$" = 1, "$i$" = 1, "$h$" = 1),
                     escape = F) %>%
    add_header_above(c(" ", "Fixed" = 4, "Random $\\\\sigma$" = 2),
                     escape = F)  %>%
    # column_spec(1, "0.4in") %>%
    # column_spec(2:7, "0.65in") %>%
    kable_styling(latex_options="scale_down")
}
