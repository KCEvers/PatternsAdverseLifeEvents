##### SUPPORTING FUNCTIONS FOR ACCUMULATION OF ADVERSE LIFE EVENT COUNTS #####

# Main function
run_accumulation = function(name_dataset, filepath_base, P, rerun = F){
  
  # Prepare data
  datalist = prepare_data(name_dataset, filepath_base, P$event_dicts, rerun = rerun)
  
  # Number of trials for Polya urn = number of possible events
  if (name_dataset == "SHP"){
    n_trials = P$event_dicts[[name_dataset]] %>% dplyr::filter(valence == "negative", 
                                                               !grepl("PL01", var_name), # Remove all personal accidents and illnesses, as only one could occur
                                                               !grepl("PL36", var_name) # Remove "other" events
    ) %>% nrow() + 1 # Personal accident or illness
  } else if (name_dataset == "HILDA"){
    n_trials = P$event_dicts[[name_dataset]] %>% dplyr::filter(valence == "negative") %>% nrow()
  }
  
  # Create dataframe with only complete data for 20 years
  df_accum <<- datalist$df_nr_negevents_pp_py %>%
    dplyr::filter(nr_years_missing_in_timespan == 0, nr_years_obs >= 20, t_id <= 20) %>%
    dplyr::select(p_id, crosswave_h_id, t_id, nr_occur) %>%
    arrange(p_id, t_id) %>%
    dplyr::mutate(nr_no_occur = n_trials - nr_occur)
  
  # Because of an error in the insight package, icc() and r2() sometimes don't work because the function cannot retrieve the data from the environment when creating the null model. Assign data to global environment:
  # assign("df_accum", df_accum, envir = .GlobalEnv)
  
  # Fit models
  mod_poisson <- glmmTMB(
    nr_occur ~ 1,
    data = df_accum, family=poisson(link = "log"))
  
  mod_frailty <- glmmTMB(
    nr_occur ~ 1 + (1 | p_id) + (1 | crosswave_h_id),
    data = df_accum, family=poisson(link = "log"))    
  
  mod_polya <- glmmTMB(
    cbind(nr_occur, nr_no_occur) ~ 1 + (1 | p_id) + (1 | crosswave_h_id),
    data = df_accum, family=betabinomial(link = "logit"))
  
  # Transform estimates to alpha and beta
  mu = plogis(mod_polya$fit$par[["beta"]]) # mu = a / (a + b) = prob
  phi = exp(mod_polya$fit$par[["betadisp"]]) # phi = a + b
  alpha = mu * phi
  beta = phi - alpha
  
  # Number of individuals and households
  nr_people = ranef(mod_frailty)$cond$p_id %>% nrow()
  nr_households = ranef(mod_frailty)$cond$crosswave_h_id %>% nrow()
  
  
  # Model fit
  df_ABIC = cbind(AIC(mod_poisson, mod_frailty, mod_polya),
                  BIC(mod_poisson, mod_frailty, mod_polya))
  
  model_names <- c("poisson", "frailty", "polya")  # Label models
  df_comp_AIC = outer(df_ABIC$AIC %>% stats::setNames(model_names), 
                      df_ABIC$AIC %>% stats::setNames(model_names), FUN = "-")
  df_comp_BIC = outer(df_ABIC$BIC %>% stats::setNames(model_names), 
                      df_ABIC$BIC %>% stats::setNames(model_names), FUN = "-")
  
  
  # Summarise fit and check model assumptions
  mod_summ_poisson = summarise_model(mod_poisson)
  mod_summ_frailty = summarise_model(mod_frailty)
  mod_summ_polya = summarise_model(mod_polya)
  
  return(list(datalist = datalist,
              df_accum = df_accum,
              nr_people = nr_people,
              nr_households = nr_households,
              mod_poisson = mod_poisson,
              mod_frailty = mod_frailty,
              mod_polya = mod_polya,
              mu = mu, phi = phi, alpha = alpha, beta = beta,
              mod_summ_poisson = mod_summ_poisson,
              mod_summ_frailty = mod_summ_frailty,
              mod_summ_polya = mod_summ_polya,
              df_ABIC = df_ABIC,
              df_comp_AIC = df_comp_AIC,
              df_comp_BIC = df_comp_BIC
  ))
}

# Format estimates of accumulation analysis
format_accumulation = function(df_est, effect, group = NULL, transform = ""){
  df = df_est %>% 
    dplyr::filter(effect == !!effect)  %>%
    dplyr::rename(lower = "2.5 %", upper = "97.5 %") %>%
    dplyr::mutate(term = stringr::str_replace_all(term, "\\(Intercept\\)", "intercept"))
  
  if (effect == "ran_pars"){
    df = df %>% dplyr::filter(group == !!group)
  }
  
  if (transform == "exp"){
    df = df %>%
      dplyr::mutate_at(c("estimate", "lower", "upper"), exp)
  } else if (transform == "prob"){
    df = df %>%
      dplyr::mutate_at(c("estimate", "lower", "upper"), ~ exp(.x) / (1 + exp(.)))
  }
  
  df %>% mutate(label = sprintf("%.2f [%.2f, %.2f]", estimate, lower, upper)) %>%
    pull(label) %>% return()
  
}

# Create dataframe with estimates of accumulation analysis
accumulation_create_df = function(name_dataset, dataset){
  
  data.frame(dataset = name_dataset,
             poisson_intercept = format_accumulation(dataset$mod_summ_poisson$df_est, "fixed", transform = "exp"),
             frailty_intercept = format_accumulation(dataset$mod_summ_frailty$df_est, "fixed", transform = "exp"),
             frailty_sd_person = format_accumulation(dataset$mod_summ_frailty$df_est, "ran_pars", group = "p_id"),
             frailty_sd_household = format_accumulation(dataset$mod_summ_frailty$df_est, "ran_pars", group = "crosswave_h_id"),
             delta_AIC_poisson_frailty = round(AIC(dataset$mod_poisson) - AIC(dataset$mod_frailty), 2),
             delta_BIC_poisson_frailty = round(BIC(dataset$mod_poisson) - BIC(dataset$mod_frailty), 2),
             polya_prob = format_accumulation(dataset$mod_summ_polya$df_est, "fixed", transform = "prob"),
             polya_dispersion = round(exp(dataset$mod_polya$fit$par[["betadisp"]]), 2),
             polya_sd_person = format_accumulation(dataset$mod_summ_polya$df_est, "ran_pars", group = "crosswave_h_id"),
             polya_sd_household = format_accumulation(dataset$mod_summ_polya$df_est, "ran_pars", group = "crosswave_h_id"),            
             delta_AIC_frailty_polya = round(AIC(dataset$mod_frailty) - AIC(dataset$mod_polya), 2),
             delta_BIC_frailty_polya = round(BIC(dataset$mod_frailty) - BIC(dataset$mod_polya), 2)) 
  
}

# Create table with estimates of accumulation analysis
accumulation_table = function(SHP, HILDA){
  
  rbind(accumulation_create_df("SHP", SHP), accumulation_create_df("HILDA", HILDA)) %>% kbl(format = "latex",
                                                                                            escape=F,
                                                                                            caption = "Estimates of the Poisson, frailty, and Polya urn models in the accumulation analysis (profile 95\\% profile confidence intervals in brackets). Estimates are on the natural scale. Fit comparisons subtract the fit of the model on the right from that of the model on the left. $\\lambda$ = Yearly rate of adverse life events; $p$ = Per-trial probability of adverse life events; $\\phi$ = Dispersion; $i$ = Individual; $h$ = Household (Source: Swiss Household Panel, SHP and Household, Income and Labour Dynamics in Australia Survey, HILDA).",
                                                                                            label = "accumulation",
                                                                                            booktabs = TRUE,
                                                                                            col.names = NULL) %>%
    add_header_above(c(" ", "$\\\\lambda$" = 1, "$\\\\lambda$" = 1, "$i$" = 1, "$h$" = 1, "$\\\\Delta$AIC" = 1, "$\\\\Delta$BIC" = 1,
                       "$p$" = 1, "$\\\\phi$" = 1, "$i$" = 1, "$h$" = 1, "$\\\\Delta$AIC" = 1, "$\\\\Delta$BIC" = 1),
                     escape = F) %>%
    add_header_above(c(" ", "Fixed" = 1, "Fixed" = 1, "Random $\\\\sigma$" = 2, "Fit" = 2,  "Fixed" = 2, "Random $\\\\sigma$" = 2, "Fit" = 2),
                     escape = F)  %>%
    add_header_above(c(" ", "Poisson" = 1, "Frailty" = 5, "Polya urn" = 6)) %>%
    column_spec(1, "0.4in") %>%
    column_spec(2:13, "0.27in") %>%
    kable_styling()
}

# Formulation simulation results
format_sim = function(sim_df, df_accum, sim_type){
  # To each list, bind participant and time information
  sim_df_merged = purrr::map(1:ncol(sim_df), function(i){
    cbind(df_accum[,c("p_id", "t_id")], nr_occur = sim_df[,i], i = i)
  }) %>% do.call(rbind, .) %>% as.data.frame()
  
  sim_df_merged %>% 
    arrange(i, p_id, t_id) %>%
    group_by(i, p_id) %>%
    mutate(cumsum_nr_occur = cumsum(nr_occur)) %>%
    ungroup() %>%
    group_by(t_id, cumsum_nr_occur) %>%
    dplyr::summarise(frequency = n(), .groups = 'drop') %>%
    group_by(t_id) %>%
    dplyr::mutate(frequency_norm = frequency / sum(frequency)) %>%
    ungroup() %>%
    mutate(sim_type = !!sim_type) %>% return()
}

# Simulate data from poisson, frailty, and Polya urn model
simulate_data = function(models, n_sim, sim_types){
  # Simulate data from all three models
  sim_poisson <- simulate(models$mod_poisson, nsim = n_sim) %>% as.matrix() %>% 
    format_sim(., models$df_accum, "Poisson")
  sim_frailty <- simulate(models$mod_frailty, nsim = n_sim) %>% as.matrix() %>%
    format_sim(., models$df_accum, "Frailty")
  sim_polya <- purrr::map(1:n_sim, function(i){simulate(models$mod_polya)}) %>%
    # Only keep number of successes
    purrr::map(function(x){as.matrix(as.matrix(x)[,1], ncol = 1)}) %>% 
    do.call(cbind, .) %>% format_sim(., models$df_accum, "Polya urn")
  
  df_accum_plot = models$df_accum %>% 
    arrange(p_id, t_id) %>%
    mutate(cumsum_nr_occur = cumsum(nr_occur), .by = p_id) %>%
    group_by(t_id, cumsum_nr_occur) %>%
    dplyr::summarise(frequency = n(), .groups = 'drop') %>%
    group_by(t_id) %>%
    dplyr::mutate(frequency_norm = frequency / sum(frequency),
                  sim_type = "Empirical") %>%
    ungroup() %>%
    rbind(sim_poisson, sim_frailty, sim_polya) %>%
    mutate(sim_type = factor(sim_type, levels = sim_types, ordered = F)) %>%
    mutate(time_label = sprintf("Year %02d", t_id) %>% stringr::str_wrap(width = 9)) 
  
  return(df_accum_plot)
}


# Plot cumulative distribution
plot_cum_distr = function(name_dataset, df, P, col_values, fill_values, size_text = 10, chosen_t_ids = c(1,5,10,15,20),
                      clip_to_emp = T){
  
  df = df %>% filter(t_id %in% chosen_t_ids) %>%
    dplyr::mutate(dataset = name_dataset)
  
  
  # Clip simulated data to range of observed empirical data
  if (clip_to_emp){
    max_emp = df %>% dplyr::filter(sim_type == "Empirical") %>%
      group_by(t_id) %>% dplyr::summarise(max = max(cumsum_nr_occur), .groups = 'drop') %>%
      pull(max, t_id)
    
    df = df %>% group_by(t_id) %>%
      dplyr::filter(cumsum_nr_occur <= max_emp[names(max_emp) == as.character(unique(t_id))] ) %>%
      dplyr::ungroup()
  }
  
  pl_distr = df %>%
    filter(
      # frequency_norm != 0,
           sim_type != "Empirical") %>%
    ggplot() +
    
    geom_line( aes(x=cumsum_nr_occur, y=frequency_norm, col = sim_type) ,alpha = 1, linewidth=1.1) +
    
    # geom_segment( aes(x=cumsum_nr_occur, xend=cumsum_nr_occur, y=0, yend=frequency_norm, col = sim_type) ,alpha = .7, linewidth=.5) +
    geom_segment(data = df %>% filter(sim_type == "Empirical"), 
                 aes(x=cumsum_nr_occur, xend=cumsum_nr_occur, y=0, yend=frequency_norm, col = sim_type) ,alpha = .7, linewidth=.5) +
    # geom_point(aes(x = cumsum_nr_occur, y = frequency_norm, col = sim_type), alpha = .8, size = 1.75) +
    geom_point(data = df %>% filter(sim_type == "Empirical"), aes(x = cumsum_nr_occur, y = frequency_norm, col = sim_type), size = 1.5, alpha = .8) +
    ggh4x::facet_grid2(time_label ~ dataset,  
                       axes = "all", switch = "y",
                       independent = "x",
                       # remove_labels = "x",
                       scales = "free", 
                       # strip.position = "left", ncol = 1
    ) +
    scale_color_manual(name = "", values = col_values, breaks = names(col_values) ) +
    # scale_fill_manual(name = "", values = fill_values , breaks = names(fill_values) ) +
    P$own_theme +
    scale_y_continuous(position = "right", expand = expansion(mult = c(0.02, .08)), n.breaks = 3) +
    scale_x_continuous(n.breaks = 4,
                       expand = expansion(mult = c(.015, .015))) +
    labs(y = "Frequency (normalized per year)", x = "Number of cumulative events")   +
    theme(legend.position = "top") +
    theme(
      panel.grid.major = element_line(color = "gray90", size = 0.3),
      panel.border = element_rect(colour = "grey30", fill=NA, linewidth=1),        panel.spacing.y = unit(0.75, "lines"),
      panel.spacing.x = unit(.01, "lines"),
      axis.text.x = element_text(size=size_text),
      strip.text.y = element_text(size=size_text),
      axis.text.y = element_text(size=size_text, vjust = 0),
      axis.title = element_text(size=size_text+4),
      legend.text = element_text(size=size_text+4)
    ) +
    theme(plot.margin = unit(c(0,0.5,0,0), 'cm'))
  
  pl_distr
  
  return(pl_distr)
}


# Plot cumulative distribution
plot_cum_distr_animate = function(name_dataset, df, P, col_values, fill_values, 
                                  clip_to_emp = T){
  
  size_text = 40
  
  # Clip simulated data to range of observed empirical data
  if (clip_to_emp){
    max_emp = df %>% dplyr::filter(sim_type == "Empirical") %>%
      group_by(t_id) %>% dplyr::summarise(max = max(cumsum_nr_occur), .groups = 'drop') %>%
      dplyr::pull(max, t_id)
    
    df = df %>% group_by(t_id) %>%
      dplyr::filter(cumsum_nr_occur <= max_emp[names(max_emp) == as.character(unique(t_id))] ) %>%
      dplyr::ungroup()
  }
  
  df = df %>% dplyr::mutate(sim_type = factor(sim_type, levels = levels(.data$sim_type), ordered = T ) ) 

  pl_distr = df %>%
    dplyr::filter(sim_type != "Empirical") %>%
    dplyr::arrange(sim_type, t_id, cumsum_nr_occur) %>%
    group_by(t_id, sim_type)%>%
    # dplyr::mutate(frequency_norm_smooth = stats::smooth.spline(cumsum_nr_occur, frequency_norm,spar= 0.4)$y) %>%
    # dplyr::mutate(frequency_norm_smooth = lowess(cumsum_nr_occur, frequency_norm,f= 0.9)$y) %>%
    # dplyr::mutate(frequency_norm_smooth = zoo::rollmean(frequency_norm, k = 2,
    #                                                     c(frequency_norm[1], mean(frequency_norm), dplyr::last(frequency_norm)))) %>%
    # dplyr::mutate(frequency_norm_smooth = forecast::ses(frequency_norm)) %>%
    # mutate(frequency_norm_smooth = c( as.numeric(stats::HoltWinters(frequency_norm, alpha = 0.8, beta = FALSE, gamma = FALSE)$fitted[, "xhat"]), dplyr::last(frequency_norm))) %>%
    
    do({
      # Create a finer grid for interpolation
      x_new <- seq(min(.$cumsum_nr_occur), max(.$cumsum_nr_occur), length.out = 1000)
      # Spline interpolation
      spline_fit <- spline(.$cumsum_nr_occur, .$frequency_norm, xout = x_new, method = "fmm")
      data.frame(
        cumsum_nr_occur = spline_fit$x,
        frequency_norm_smooth = spline_fit$y,
        sim_type = unique(.$sim_type)
      )
    }) %>%
    
    dplyr::ungroup() %>%
    # dplyr::mutate(dataset = name_dataset) %>% group_by(sim_type, t_id) %>%
    # dplyr::mutate(frequency_norm_norm = frequency_norm / max(frequency_norm)) %>%
    # dplyr::ungroup() %>%
    # dplyr::filter(
    #   # frequency_norm != 0, 
    #               sim_type != "Empirical") %>%
    ggplot() +

    geom_line( aes(x=cumsum_nr_occur, y=frequency_norm_smooth, col = sim_type) ,alpha = 1, linewidth=1.5) +
    # geom_smooth( aes(x=cumsum_nr_occur, y=frequency_norm, col = sim_type) ,alpha = 1, linewidth=1.5,
    #              se=FALSE, method="loess") +
    
  geom_bar(data = df %>%
                 dplyr::filter(sim_type == "Empirical"),
               aes(x=cumsum_nr_occur, y=frequency_norm), stat = 'identity',
           col = col_values[names(col_values)== "Empirical"],
           fill = fill_values[names(fill_values)== "Empirical"],
               position = "identity", alpha = .5, linewidth=.5) +
        # 
  # geom_segment(data = df %>% filter(sim_type == "Empirical"), 
  #              aes(x=cumsum_nr_occur, xend=cumsum_nr_occur, y=0, yend=frequency_norm, col = sim_type) ,alpha = .7, linewidth=.5) +
  # geom_point(data = df %>% filter(sim_type == "Empirical"), aes(x = cumsum_nr_occur, y = frequency_norm, col = sim_type), size = 1.5, alpha = .8) +
  # ggh4x::facet_grid2(time_label ~ dataset,  
  #                    axes = "all", switch = "y",
  #                    independent = "x",
  #                    # remove_labels = "x",
  #                    scales = "free", 
  #                    # strip.position = "left", ncol = 1
  # ) +
    # ggh4x::facet_grid2(. ~ time_label) +
  scale_color_manual(name = "", values = col_values[names(col_values) != "Empirical"] ) +
    scale_fill_manual(name = "", values =  fill_values[names(fill_values) != "Empirical"] ) +
    P$own_theme +
    scale_y_continuous( expand = expansion(mult = c(0.0, 0.02)), n.breaks = 3) +
    scale_x_continuous(n.breaks = 4,
                       limits = c(-.5, 75),
                       expand = expansion(mult = c(0, 0))) +
    labs(y = "Frequency (normalized)", x = "Number of cumulative adverse life events")   +
    theme(legend.position = "bottom") +
    theme(
      # panel.grid.major = element_line(color = "gray90", size = 0.3),
      # panel.border = element_rect(colour = "grey30", fill=NA, linewidth=1), 
      panel.border = element_blank(), 
      panel.spacing.y = unit(0.75, "lines"),
      
     panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      panel.spacing.x = unit(.01, "lines"),
      axis.text.x = element_text(size=size_text),
      strip.text.y = element_text(size=size_text),
      # axis.text.y = element_text(size=size_text, vjust = 0),
      axis.text.y = element_blank(),
     axis.ticks.y = element_blank(),
     axis.title = element_text(size=size_text+4),
     title = element_text(size=size_text+10),
     legend.text = element_text(size=size_text+4)
    ) +
    theme(plot.margin = unit(c(0,0.5,0,0), 'cm'))
  
  pl_distr
  
  anim = pl_distr + gganimate::transition_time(t_id) + gganimate::view_follow(fixed_x = T) +
    labs(title ='Year: {frame_time}') + theme(
      title = element_text(size=size_text+10)
    )
  gganimate::animate(anim, nframes = 100, height = 12, width = 18, units = "cm", res = 300)
  
  filepath = file.path(filepath_base, "figs", name_dataset, "cumulative_nr_events.gif")
  gganimate::anim_save(filepath)
  
  
  return(filepath)
}


# Fit heavy-tailed distributions to cumulative counts at twenty years
fit_heavy_tails = function(counts, set_xmin = NULL){

  # Step 2: Fit a power-law distribution
  pl_model <- displ$new(counts)  # Discrete Power-law, 1 parameter
  if (is.null(set_xmin)){
    est <- estimate_xmin(pl_model)  # Estimate xmin and alpha
  } else {
    est = set_xmin
  }
  pl_model$setXmin(est)
  pl_model$setPars(estimate_pars(pl_model))
  cat("Power-law fit: xmin =", pl_model$getXmin(), ", alpha =", pl_model$pars, "\n")

  # Step 3: Fit alternative distributions
  # Lognormal
  ln_model <- dislnorm$new(counts) # Discrete Log-normal, 2 parameters
  if (is.null(set_xmin)){
    ln_est <- estimate_xmin(ln_model)
  } else {
    ln_est = set_xmin
  }
  ln_model$setXmin(ln_est)
  ln_model$setPars(estimate_pars(ln_model))
  cat("Lognormal fit: xmin =", ln_model$getXmin(), ", meanlog =", ln_model$pars[1],
      ", sdlog =", ln_model$pars[2], "\n")

  # Exponential
  exp_model <- disexp$new(counts) # Discrete Exponential, 1 parameter
  if (is.null(set_xmin)){
    exp_est <- estimate_xmin(exp_model)
  } else {
    exp_est = set_xmin
  }
  exp_model$setXmin(exp_est)
  exp_model$setPars(estimate_pars(exp_model))
  cat("Exponential fit: xmin =", exp_model$getXmin(), ", rate =", exp_model$pars, "\n")

  # Poisson
  pois_model <- dispois$new(counts) # Poisson, 1 parameter
  if (is.null(set_xmin)){
    pois_est <- estimate_xmin(pois_model)
  } else {
    pois_est = set_xmin
  }
  pois_model$setXmin(pois_est)
  pois_model$setPars(estimate_pars(pois_model))
  cat("Poisson fit: xmin =", pois_model$getXmin(), ", rate =", pois_model$pars, "\n")

  # Step 5: Plot the empirical counts and fitted distributions
  # Main plot: Empirical counts (CDF)
  plot(pl_model, main = "Fitted Distributions", xlab = "Value", ylab = "CDF", pch = 19)

  # Add fitted distributions
  lines(pl_model, col = "red", lwd = 3, lty = 1)  # Power law
  lines(ln_model, col = "blue", lwd = 3, lty = 2)  # Lognormal
  lines(exp_model, col = "green", lwd = 3, lty = 3)  # Exponential
  lines(pois_model, col = "magenta", lwd = 3, lty = 3)  # Poisson

  # Add a legend
  legend("bottomleft",
         legend = c("Power law", "Lognormal", "Exponential", "Poisson"),
         col = c("red", "blue", "green", "magenta"),
         lty = c(1, 2, 3, 4),
         lwd = 2)

  return(list(pl_model = pl_model,
              ln_model = ln_model,
              exp_model = exp_model,
              pois_model = pois_model))

}

# Compare heavy-tailed distributions
compare_distributions_wrap = function(model1_, model2_, 
                                      model1_name, model2_name,
                                      choose_xmin = c("model1", "model2", "median")[3]){
  
  # Copy models, otherwise they are overwritten (shallow copies do not work)
  model1 = model1_$copy()
  model2 = model2_$copy()
  
  # Need same xmin to compare models: Choose median xmin of model 1 and model 2
  if (choose_xmin == "median"){
    xmin_min = round(median(c(model1$getXmin(), model2$getXmin())))
    
    # Set xmin and reestimate parameters
    model1$setXmin(xmin_min)
    model1$setPars(estimate_pars(model1))
    model2$setXmin(xmin_min)
    model2$setPars(estimate_pars(model2))
    
  } else if (choose_xmin == "model1"){
    
    # Alternatively, set to xmin of model 1
    model2$setXmin(model1$getXmin())
    model2$setPars(estimate_pars(model2))
    
  } else if (choose_xmin == "model2"){
    # Alternatively, set to xmin of model 2
    model1$setXmin(model2$getXmin())
    model1$setPars(estimate_pars(model1))
  }

  # Compare models
  comparison_1_to_2 <- compare_distributions(model1, model2)
  comparison_2_to_1 <- compare_distributions(model2, model1)
  
  res = data.frame(loglik_ratio = comparison_1_to_2$test_statistic,
             p_two_sided = comparison_1_to_2$p_two_sided,
             # Two-sided p-value doesn't indicate which model is better
p_1_better_than_2 = comparison_1_to_2$p_one_sided,
p_2_better_than_1 = comparison_2_to_1$p_one_sided
) %>%
  tidyr::gather() %>% 
  dplyr::mutate(key = sprintf("%s_vs_%s_%s", model1_name, model2_name, key))
return(res)

}

# Run fitting of all heavy-tailed distributions
run_heavy_tails = function(name_dataset, df_accum, datalist){
  
  counts = df_accum %>% 
    arrange(p_id, t_id) %>%
    group_by(p_id) %>%
    dplyr::summarise(nr_occur_cumsum = sum(nr_occur), .groups = 'drop') %>% 
    # Power-law fitting cannot work with zero counts - remove, as we're only interested in the tail behaviour anyway
    dplyr::filter(nr_occur_cumsum != 0) %>%
    dplyr::pull(nr_occur_cumsum)
  
  # Summary statistics
  min_cumsum_count_20y = min(counts)
  max_cumsum_count_20y = max(counts)
  mean_cumsum_count_20y = mean(counts)
  median_cumsum_count_20y = median(counts)
  var_to_mean_cumsum_count_20y = var(counts) / mean(counts)
  
  
  # Fit heavy-tailed distributions
  out = fit_heavy_tails(counts)
  
  # Plot
  
  # Empirical CDF
  ecdf_data <- ecdf(counts)
  x_vals <- unique(sort(counts))
  y_vals <- 1 - ecdf_data(x_vals)  # Complementary CDF (CCDF)
  
  # Empirical data
  df_emp = data.frame(x = unique(sort(counts)), 
                      y = 1 - ecdf_data(x_vals)  # Complementary CDF (CCDF)
  ) %>% dplyr::mutate(type = name_dataset)
  
  # Fitted line data
  pl_model_line = lines(out$pl_model, draw = F)
  ln_model_line = lines(out$ln_model, draw = F)
  exp_model_line = lines(out$exp_model, draw = F)
  pois_model_line = lines(out$pois_model, draw = F)
  df = rbind(
    pl_model_line %>% dplyr::mutate(type = "Power-law"),
    ln_model_line %>% dplyr::mutate(type = "Log-normal"),
    exp_model_line %>% dplyr::mutate(type = "Exponential"),
    pois_model_line %>% dplyr::mutate(type = "Poisson"))
  
  # Plot
  pl = df %>%
    dplyr::mutate(dataset = name_dataset) %>%
    ggplot() +
    geom_line(aes(x = x, y = y, col = type), linewidth = 1,
              linetype = 'dashed') +
    geom_point(data = df_emp, aes(x = x, y = y), size = 0.8, col = 'red') +
    ggh4x::facet_grid2(. ~ dataset) +
    scale_x_continuous(trans = "log10", 
                       breaks = scales::trans_breaks("log10", function(x) 10^x),   # Define breaks at powers of 10
                       labels = scales::trans_format("log10", scales::math_format(10^.x)), # Format ticks as 10^x
                       limits = range(df_emp$x)) +
    scale_y_continuous(trans = "log10", 
                       # expand = expansion(add = c(1, 0), mult = c(1, .05)),
                       breaks = scales::trans_breaks("log10", function(x) 10^x),   # Define breaks at powers of 10
                       labels = scales::trans_format("log10", scales::math_format(10^.x)), # Format ticks as 10^xhttp://127.0.0.1:18585/graphics/0a003603-149c-4235-8a42-777569f495df.png
                       # limits = range(df_emp$y)
    ) +
    labs(x = "Cumulative count of events", y = "Probability") +
    viridis::scale_color_viridis(discrete = T, name = "", end = .95) +
    P$own_theme +
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 10),
          legend.text = element_text(size = 10),
    ) #+ guides(
  # linetype =NULL,
  #                    col = guide_legend(ncol = 1,
  #                              title.position = "top"))
  pl
  
  
  # Save results
  res = data.frame(
    powerlaw_xmin = out$pl_model$xmin,
    powerlaw_alpha = out$pl_model$pars,
    lognormal_xmin = out$ln_model$xmin,
    lognormal_meanlog = out$ln_model$pars[1],
    lognormal_sdlog = out$ln_model$pars[2],
    exp_xmin = out$exp_model$xmin,
    exp_rate = out$exp_model$pars,
    pois_xmin = out$pois_model$xmin,
    pois_rate = out$pois_model$pars) %>%
    # tidyr::gather() %>%
    # Compare distributions; need to have the same xmin
    cbind(
      rbind(compare_distributions_wrap(out$ln_model, out$pl_model, "lognormal", "powerlaw"),
            compare_distributions_wrap(out$ln_model, out$exp_model, "lognormal", "exponential"),
            compare_distributions_wrap(out$ln_model, out$pois_model, "lognormal", "poisson")
      ) %>% tidyr::spread(key, value)) %>%
    dplyr::mutate(dataset = name_dataset) 
  res
  res
  
  df_est = res %>%
    dplyr::select(dataset, pois_xmin, pois_rate, exp_xmin, exp_rate,
                  lognormal_xmin, lognormal_meanlog, lognormal_sdlog,
                  powerlaw_xmin, powerlaw_alpha) %>%
    mutate_at(setdiff(colnames(.), c("dataset", "pois_xmin", "exp_xmin", "lognormal_xmin", "powerlaw_xmin")), ~ sprintf("%.2f", .)) %>%
    mutate(lognormal_vs_poisson = sprintf("%.2f ($p$ = %.2f)",
                                          res[["lognormal_vs_poisson_loglik_ratio"]],
                                          res[["lognormal_vs_poisson_p_two_sided"]]),
           lognormal_vs_exp = sprintf("%.2f ($p$ = %.2f)",
                                      res[["lognormal_vs_exponential_loglik_ratio"]],
                                      res[["lognormal_vs_exponential_p_two_sided"]]),
           lognormal_vs_powerlaw = sprintf("%.2f ($p$ = %.2f)",
                                           res[["lognormal_vs_powerlaw_loglik_ratio"]],
                                           res[["lognormal_vs_powerlaw_p_two_sided"]])
    ) %>%
    mutate_at(c("lognormal_vs_poisson",
                "lognormal_vs_exp",
                "lognormal_vs_powerlaw"), ~ stringr::str_replace_all(., "= 0.00", "$<$ .01"))
  df_est
  
  return(list(df_est = df_est, pl = pl))
}

# Create table of estimates of all heavy-tailed distributions
heavy_tails_table = function(SHP_tails, HILDA_tails){
  rbind(SHP_tails$df_est,
        HILDA_tails$df_est) %>%
    kbl(format = "latex",
        escape=F,
        caption = "Model estimates of distributions fit to twenty-year cumulative adverse life event counts (Source: Swiss Household Panel, SHP and Household, Income and Labour Dynamics in Australia Survey, HILDA).",
        booktabs = TRUE,
        label = "heavytails",
        col.names = NULL) %>%
    add_header_above(c(" ", "$x_{min}$" = 1, "$\\\\lambda$" = 1, "$x_{min}$" = 1, "$\\\\lambda$" = 1,
                       "$x_{min}$" = 1, "$\\\\mu$" = 1, "$\\\\sigma$" = 1,
                       "$x_{min}$", "$\\\\alpha$", "Poisson",  "Exponential",  "Power-law"),
                     escape = F) %>%
    add_header_above(c(" ", "Poisson" = 2, "Exponential" = 2, "Log-normal" = 3, "Power-law" = 2, "Log-Likelihood Ratio cf. Log-normal" = 3),
                     escape = F)  %>%
    kable_styling(latex_options="scale_down")
}


