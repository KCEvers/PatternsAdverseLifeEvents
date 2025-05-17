##### SUPPORTING FUNCTIONS FOR CO-OCCURRENCE ANALYSIS OF ADVERSE LIFE EVENT TYPES #####

# Run co-occurrence analysis for one dataset
run_cooccurrence = function(name_dataset, filepath_base, P, 
                            chosen_leads = c(0,1), rerun = F){
  
  # Prepare data
  datalist = prepare_data(name_dataset, filepath_base, P$event_dicts, rerun = rerun)
  
  # Reformat dataframe of events
  df_events = datalist$df_per_event_pp_py %>%
    filter(valence == "negative") %>%
    dplyr::mutate(occur = ifelse(nr_missing_occur > 0, NA, nr_occur)) %>%
    dplyr::select(p_id, crosswave_h_id, wave_nr, event, occur) %>%
    tidyr::spread(event, occur) %>%
    tidyr::complete(p_id, wave_nr) %>%
    arrange(p_id, wave_nr)
  df_events %>% head
  
  # Get event names
  negevent_names = setdiff(colnames(df_events), c("p_id", "crosswave_h_id", "wave_nr")) %>% sort()
  head(df_events)
  nrow(df_events)
  
  # Loop through leads, fit models with and without random effects, and pick converged models
  conv_model_list = list()
  conv_model_est = list()
  nr_obs = list()
  nr_timepoints_per_p_id = list()
  
  for (chosen_lead in chosen_leads){
    
    # Run linear mixed-effect model for this particular lead
    filepath = file.path(datalist$filepath_deriv, sprintf("cooccurrence_event_types_lead%d.RDS", chosen_lead))
    if (!file.exists(filepath)){
      start_t = Sys.time()
      model_list = run_glmer(df_events, negevent_names, 
                             event_dict = P$event_dicts[[name_dataset]], 
                             chosen_lead = chosen_lead,
                             name_dataset = name_dataset)
      end_t = Sys.time()
      print(end_t - start_t)
      
      saveRDS(model_list, filepath)
    }
    model_list = readRDS(filepath)
    format(object.size(model_list), "Mb")
    
    
    # Pick converged models and summarise 
    filepath = file.path(datalist$filepath_deriv, sprintf("cooccurrence_event_types_summ_lead%d.RDS", chosen_lead))
    if (!file.exists(filepath)){
      conv_model_summ = pick_model(model_list, negevent_names, datalist, name_dataset, chosen_lead)
      saveRDS(conv_model_summ, filepath)
    }
    conv_model_summ = readRDS(filepath)
    conv_model_list[[as.character(chosen_lead)]] = conv_model_summ
    
    # Number of individuals and households
    nr_obs[[as.character(chosen_lead)]] = conv_model_summ$conv_models %>% 
      purrr::map(function(x){ranef(x)$cond %>% purrr::map_vec(nrow)})
    
    # Number of timepoints per person for all models
    nr_timepoints_per_p_id[[as.character(chosen_lead)]] = purrr::map(conv_model_summ$conv_models, 
                                                                     function(x){
                                                                       if ("p_id" %in% colnames(x$frame)){
                                                                         x$frame %>% group_by(p_id) %>%
                                                                           dplyr::summarise(nr_waves = n(), .groups = 'drop') %>%
                                                                           group_by(nr_waves) %>% dplyr::summarise(nr_people = n(), .groups = 'drop')
                                                                       }
                                                                     })
    
    # Get plotting dataframe
    model_df = conv_model_summ$conv_models_summ %>%
      purrr::imap(function(x, i){get_plot_df(x, negevent_names[i], 
                                             conv_model_summ$model_check_df)}) %>% do.call(rbind, .)
    
    conv_model_est[[as.character(chosen_lead)]] = model_df
  }
  
  return(list(
    datalist = datalist,
    df_events = df_events,
    nr_obs = nr_obs,
    nr_timepoints_per_p_id = nr_timepoints_per_p_id,
    negevent_names = negevent_names,
    conv_model_list = conv_model_list,
    conv_model_est = conv_model_est
  ))
  
}


# Run linear mixed-effect model
run_glmer = function(df_events, negevent_names, event_dict, chosen_lead, name_dataset){
  
  model_list = foreach::foreach (i = 1:length(negevent_names),
                                 .packages = c("glmmTMB", "dplyr"),
                                 .export = c("df_events", "chosen_lead", "name_dataset", "event_dict")) %dopar% {
                                   
                                   negevent_name = negevent_names[i]
                                   negevent_names_ = negevent_names
                                   
                                   if (name_dataset == "SHP"){
                                     
                                     # Drop other or unspecified illness or accident as predictor for all models, as otherwise the personal illness or accident categories are perfectly collinear (being mutually exclusive subcategories)
                                     negevent_names_ = setdiff(negevent_names_, "unspecified_illness_accident")
                                     
                                     # If the outcome event is a personal illness or accident, remove other personal illnesses or accidents as predictors, as these cannot co-occur
                                     illness_accident_events = event_dict %>% dplyr::filter(grepl("PL01R", var_name)) %>% dplyr::pull(recode_var)
                                     if (negevent_name %in% illness_accident_events & chosen_lead == 0){
                                       negevent_names_ = setdiff(negevent_names_, illness_accident_events)
                                     } 
                                     
                                   }
                                   
                                   formula_no_random = paste0(negevent_name, " ~ ", paste0(negevent_names_[negevent_names_ != negevent_name], collapse = " + "))
                                   formula_p_id_h_id = paste0(formula_no_random, " + (1 | p_id) + (1 | crosswave_h_id)")
                                   formula_p_id = paste0(formula_no_random, " + (1 | p_id)")
                                   
                                   # How is this event predicted by events last year?
                                   if (chosen_lead > 0){
                                     new_name = paste0("lead_", negevent_name)
                                     
                                     df_events_ = df_events %>% group_by(p_id) %>% 
                                       mutate(!!new_name := lead(!!sym(negevent_name), chosen_lead)) %>% ungroup()
                                     
                                     # Include event of interest as predictor
                                     formula_no_random = paste0("lead_", negevent_name, " ~ ", paste0(negevent_names_, collapse = " + "))
                                     formula_p_id_h_id = paste0(formula_no_random, " + (1 | p_id) + (1 | crosswave_h_id)")
                                     formula_p_id = paste0(formula_no_random, " + (1 | p_id)")
                                     
                                   } else {
                                     df_events_ = df_events
                                   }
                                   
                                   start_t = Sys.time()
                                   
                                   model_p_id_h_id <- glmmTMB(
                                     as.formula(formula_p_id_h_id),
                                     data = df_events_, family=binomial(link = "logit"),
                                   )
                                   
                                   model_p_id <- glmmTMB(
                                     as.formula(formula_p_id),
                                     data = df_events_, family=binomial(link = "logit"),
                                   )        
                                   
                                   model_no_random <- glmmTMB(
                                     as.formula(formula_no_random),
                                     data = df_events_, family=binomial(link = "logit")
                                   )
                                   
                                   end_t = Sys.time()
                                   print(end_t - start_t)
                                   
                                   return(list(
                                     model_p_id_h_id = model_p_id_h_id,
                                     model_p_id = model_p_id,
                                     model_no_random = model_no_random
                                   ))
                                 } 
  return(model_list)
}

# Check singularity and convergence
model_check = function(model){
  return(list(singularity = performance::check_singularity(model),
              no_convergence = !performance::check_convergence(model),
              AIC = AIC(model),
              BIC = BIC(model)))
}


pick_model = function(model_list, negevent_names, datalist, name_dataset, chosen_lead){
  
  # Check singularity and convergence
  model_check_df = model_list %>%
    purrr::map_depth(., 2, model_check) %>%
    purrr::imap(function(x, i){cbind(response = negevent_names[i], 
                                     model_type = names(x), dplyr::bind_rows(x))}) %>%
    do.call(rbind, .) %>% as.data.frame()
  
  # Plot AIC
  pl = model_check_df %>%
    ggplot() + geom_point(aes(x = model_type, y = AIC,col=model_type), size = 2) + ggh4x::facet_wrap2(response ~ ., scales = "free_y", ncol = 3)  + P$own_theme + theme(strip.text = element_text(size = 10), axis.text.y = element_text(size = 10), axis.text.x = element_text(size = 10, angle = 90), legend.position = 'none')
  pl
  
  filepath_image = file.path(datalist$filepath_figs_dataset, sprintf("%s_AIC_all_models_lead%d.pdf", name_dataset, chosen_lead))
  save_plot(pl, filepath_image, height = 200)
  
  # Plot BIC
  pl = model_check_df %>%
    ggplot() + geom_point(aes(x = model_type, y = BIC,col=model_type), size = 2) + ggh4x::facet_wrap2(response ~ ., scales = "free_y", ncol = 3)  + P$own_theme + theme(strip.text = element_text(size = 10), axis.text.y = element_text(size = 10), axis.text.x = element_text(size = 10, angle = 90), legend.position = 'none')
  pl
  
  filepath_image = file.path(datalist$filepath_figs_dataset, sprintf("%s_BIC_all_models_lead%d.pdf", name_dataset, chosen_lead))
  save_plot(pl, filepath_image, height = 200)
  
  # Find outcomes for which the full model did not converge
  unconverged_models = model_check_df %>% 
    filter(model_type == "model_p_id_h_id", 
           singularity == TRUE | no_convergence == TRUE) %>%
    pull(response) %>% match(negevent_names)
  
  # For unconverged models, first choice is dropping the household random effect - check whether these models converged
  converged_p_id = model_check_df %>% filter(model_type == "model_p_id", 
                                             singularity == FALSE & no_convergence == FALSE)
  stopifnot("Not all unconverged models reach convergence by dropping the random intercept for household!" = all(negevent_names[unconverged_models] %in% converged_p_id$response))
  
  # Create list of final converged models    
  conv_models = model_list %>% purrr::map("model_p_id_h_id")
  conv_models[unconverged_models] = model_list[unconverged_models] %>% purrr::map("model_p_id")
  
  # Summarise converged models (parallel doesn't seem to work)
  conv_models_summ = foreach::foreach (i = 1:length(negevent_names)) %do% {
    summarise_model(conv_models[[i]])
  }
  
  # Compare AIC and BIC across models
  comp_AIC_BIC = function(x, AIC_or_BIC){
    df_comp = outer(x[[AIC_or_BIC]] %>% stats::setNames(x$model_type), 
                    x[[AIC_or_BIC]] %>% stats::setNames(x$model_type), FUN = "-")
    return(df_comp)
  }
  
  df_comp_AIC = model_check_df %>% group_by(response) %>%
    group_map(~ comp_AIC_BIC(.x, AIC_or_BIC = "AIC"), .keep = T) %>% 
    setNames(negevent_names)
  
  df_comp_BIC = model_check_df %>% group_by(response) %>%
    group_map(~ comp_AIC_BIC(.x, AIC_or_BIC = "BIC"), .keep = T) %>% 
    setNames(negevent_names)
  
  return(list(
    # model_ci = model_ci,
    model_check_df = model_check_df,
    unconverged_models = unconverged_models,
    conv_models = conv_models,
    conv_models_summ = conv_models_summ,
    df_comp_AIC = df_comp_AIC,
    df_comp_BIC = df_comp_BIC
  ))
  
}

# Create dataframe with model estimates
get_plot_df = function(model_summ, negevent_name, model_check_df){
  
  # Get fixed effects estimates
  model_est = model_summ$df_est %>% 
    filter(effect == "fixed") %>%
    dplyr::rename(lower = "2.5 %", upper = "97.5 %") %>%
    dplyr::select(all_of(c("term", "estimate", "lower", "upper"))) %>%
    dplyr::rename(statistic = term) %>%
    mutate(subplot = "fixed")
  
  # Get SD of random intercepts - note that not all models have a household intercept
  sd_intercept = model_summ$df_est %>% 
    filter(effect == "ran_pars") %>%
    dplyr::rename(lower = "2.5 %", upper = "97.5 %") %>%
    dplyr::mutate(term = paste0(group, "_", term) %>% stringr::str_replace_all(., "\\(Intercept\\)", "intercept") %>%
                    stringr::str_replace_all(., "__", "_")) %>%
    dplyr::select(all_of(c("term", "estimate", "lower", "upper"))) %>%
    dplyr::rename(statistic = term) %>%
    mutate(subplot = "sd_intercept")
  
  # Get variance partitioning statistics
  var_part = rbind(
    model_summ[["ICC"]] %>% as.data.frame() %>% dplyr::select("ICC_adjusted") %>% t(),
    model_summ[["R2"]] %>% as.data.frame() %>% dplyr::select(-optional) %>% t()
  ) %>% magrittr::set_colnames("estimate") %>% as.data.frame() %>% 
    tibble::rownames_to_column("statistic") %>%
    mutate(subplot = "var_part")
  
  # Get AIC and BIC comparisons
  AIC_BIC = model_check_df %>% dplyr::filter(response == !!negevent_name) %>%
    dplyr::select(model_type, AIC, BIC) %>%
    tidyr::pivot_wider(names_from = model_type, values_from = c(AIC, BIC)) %>%
    tidyr::gather(statistic, estimate) %>%
    mutate(subplot = "fit")
  
  model_df = cbind(response = negevent_name, 
                   dplyr::bind_rows(model_est, sd_intercept, var_part, AIC_BIC))
  
  return(model_df)
}


# Plot co-occurrence of event types
plot_cooccur = function(name_dataset, datalist, P, model_df, chosen_lead){
  
  plot_df = model_df
  if (chosen_lead == 0){
    # Remove intercept
    plot_df = plot_df %>%
      dplyr::rowwise() %>%
      mutate(statistic = ifelse(statistic == "(Intercept)", response, statistic)) %>%
      ungroup %>%
      dplyr::filter(response != statistic)
    ylab = "Outcome"
  } else if (chosen_lead == 1){
    plot_df = plot_df %>%
      dplyr::filter(statistic != "(Intercept)")
    ylab = "Outcome (next year)"
  }
  
  
  # Dataframe of fixed effects
  fixed_df = plot_df %>%
    filter(subplot == "fixed") %>%
    # Exponentiate to obtain odds ratio
    dplyr::mutate_at(c("estimate", "lower", "upper"), ~ round(exp(.), 2)) %>%
    dplyr::mutate(excludes_1 = (lower < 1 & upper < 1) | (lower > 1 & upper > 1)) %>%
    dplyr::mutate(
      fill = ifelse(excludes_1, estimate, NA),
      label = sprintf("%.2f\n[%.2f, %.2f]", estimate, lower, upper),
      est = sprintf("%.2f\n", estimate),
      ci = sprintf("\n[%.2f, %.2f]", lower, upper)
    ) %>%
    # Rename events with event description
    recode_events(., c("response", "statistic"), P$event_dicts[[name_dataset]], datalist$df_per_event)
  
  # SHP: make "unspecified_illness_accident" last category because it is the reference category for personal illness or accident
  if (name_dataset == "SHP"){
    fixed_df = fixed_df %>%
      mutate_at(c("response"), 
                ~ forcats::fct_relevel(.x, "Other or unspecified illness or accident"))
    
  }
  
  # Set plotting theme
  size_text = ifelse(name_dataset == "SHP", 9.5, 9.5)
  adapt_theme = P$own_theme +
    theme(legend.position = 'bottom') +
    theme(axis.text.x=element_text(size=size_text, angle=90,hjust = 0,vjust = 1),
          axis.text.y=element_text(size=size_text),
          axis.title.x.top=element_text(size=size_text_theme+4, margin = margin(b = 20)),
          axis.title.y=element_text(size=size_text_theme+4, margin = margin(r = 20)),
          legend.text=element_text(size=size_text),
          legend.title=element_text(size=size_text+4),
          axis.title=element_text(size=size_text+4),
          plot.margin = unit(c(.01, .15, .01, .01), 'cm'),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.spacing = unit(0, "cm"),    # Remove space between legend items
          legend.margin = margin(0, 0, 0, 0))
  size_text = ifelse(name_dataset == "SHP", 1.8, 1.5)
  width_legendbar = 10
  height_legendbar = 1
  vjust = 0.6
  
  # Plot fixed effects
  pl_fixed = fixed_df %>%
    ggplot() + 
    geom_tile(aes(x = statistic, y = response,
                  col = fill,
                  fill = fill, 
    ), alpha = .95, linewidth = .5) +
    
    # Significant estimates, plot estimate and CI separately
    geom_text(data = fixed_df %>% filter(excludes_1 == TRUE), aes(x = statistic, y = response, label= ifelse(is.na(estimate), "", ci)), 
              # fontface = "bold",
              hjust = 0.5, vjust = vjust, size = size_text, family = P$font_family) +
    
    geom_text(data = fixed_df %>% filter(excludes_1 == TRUE), aes(x = statistic, y = response, label= ifelse(is.na(estimate), "", est)), 
              fontface = "bold",
              hjust = 0.5, vjust = vjust, size = size_text+1.1, family = P$font_family) +
    
    # Non-significant estimates, plot estimate and CI separately
    geom_text(data = fixed_df %>% filter(excludes_1 == FALSE), aes(x = statistic, y = response, label= ifelse(is.na(estimate), "", ci)), 
              # fontface = "bold",
              color = "grey70",
              hjust = 0.5, vjust = vjust, size = size_text, family = P$font_family) +
    
    geom_text(data = fixed_df %>% filter(excludes_1 == FALSE), aes(x = statistic, y = response, label= ifelse(is.na(estimate), "", est)), 
              fontface = "bold",
              color = "grey70",
              hjust = 0.5, vjust = vjust, size = size_text+1.1, family = P$font_family) +
    
    scale_color_gradient(name = "Adjusted odds ratio",
                         # low="blue", mid = "gray", high="red", na.value = "white", # scale_color_gradient2()
                         low="yellow", high="red", na.value = "white",
                         # midpoint = 1
    ) +
    scale_fill_gradient(name = "Adjusted odds ratio",
                        # low="blue", mid = "gray", high="red", na.value = "white", # scale_color_gradient2()
                        low="yellow", high="red", na.value = "white",
                        # midpoint = 1
    ) +
    guides(col = NULL, fill = guide_colourbar(title.position = "top",
                                              theme = theme(
                                                legend.key.width  = unit(width_legendbar, "lines"),
                                                legend.key.height = unit(height_legendbar, "lines")
                                              ))) +
    labs(y = ylab, x = "Predictor") +
    adapt_theme +
    scale_x_discrete(labels = function(x) sapply(x, wrap_equally, max_chars = 31), limits=rev, position = 'top',
                     expand = expansion(add = c(0,0))) +
    scale_y_discrete(labels = function(x) sapply(x, wrap_equally, max_chars = 31),
                     expand = expansion(add = c(0,0)))  
  pl_fixed
  
  if (name_dataset == "HILDA"){
    height =  200 
  } else if (name_dataset == "SHP"){
    height = 180 #  200
  }
  
  filepath_image = file.path(datalist$filepath_figs_dataset, sprintf("%s_heatmap_odds_ratio_lead%d.pdf", name_dataset, chosen_lead))
  save_plot(pl_fixed, filepath_image, height = height)
  
}


# Plot co-occurrence of selected event types
plot_selected_cooccur = function(name_dataset, datalist, P, model_df, chosen_lead){
  
  plot_df = model_df
  if (chosen_lead == 0){
    # Remove intercept
    plot_df = plot_df %>%
      dplyr::rowwise() %>%
      mutate(statistic = ifelse(statistic == "(Intercept)", response, statistic)) %>%
      ungroup %>%
      dplyr::filter(response != statistic)
    ylab = "Outcome"
  } else if (chosen_lead == 1){
    plot_df = plot_df %>%
      dplyr::filter(statistic != "(Intercept)")
    ylab = "Outcome (next year)"
  }
  
  selected_events = c("jailed", "separated_from_spouse", "victim_physical_violence", "fired", "injury_illness_self")
  
  # Dataframe of fixed effects
  fixed_df = plot_df %>%
    filter(subplot == "fixed") %>%
    # Exponentiate to obtain odds ratio
    dplyr::mutate_at(c("estimate", "lower", "upper"), ~ round(exp(.), 2)) %>%
    dplyr::mutate(excludes_1 = (lower < 1 & upper < 1) | (lower > 1 & upper > 1)) %>%
    dplyr::mutate(
      fill = ifelse(excludes_1, estimate, NA),
      label = sprintf("%.2f\n[%.2f, %.2f]", estimate, lower, upper),
      est = sprintf("%.2f\n", estimate),
      ci = sprintf("\n[%.2f, %.2f]", lower, upper)
    ) %>%
    dplyr::filter(.data$response %in% selected_events, .data$statistic %in% selected_events) %>%
    # Rename events with event description
    recode_events(., c("response", "statistic"), P$event_dicts[[name_dataset]], datalist$df_per_event) 
  
  
  # Set plotting theme
  size_text_theme = ifelse(name_dataset == "SHP", 9.5 * 1.6, 9.5 * 1.6)
  adapt_theme = P$own_theme +
    theme(legend.position = 'bottom') +
    theme(axis.text.x=element_text(size=size_text_theme, angle=90,hjust = 0,vjust = 1),
          axis.text.y=element_text(size=size_text_theme),
          legend.text=element_text(size=size_text_theme-4),
          legend.title=element_text(size=size_text_theme),
          axis.title.x.top=element_text(size=size_text_theme+4, margin = margin(b = 20)),
          axis.title.y=element_text(size=size_text_theme+4, margin = margin(r = 20)),
          plot.margin = unit(c(.01, .15, .01, .01), 'cm'),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.spacing = unit(0, "cm"),    # Remove space between legend items
          legend.margin = margin(0, 0, 0, 0))
  size_text = ifelse(name_dataset == "SHP", 1.8 * 2.6, 1.5 * 2.6)
  width_legendbar = 25
  height_legendbar = 1
  vjust = 0.6
  
  # Plot fixed effects
  pl_fixed = fixed_df %>%
    ggplot() + 
    geom_tile(aes(x = statistic, y = response,
                  col = fill,
                  fill = fill, 
    ), alpha = .95, linewidth = .5) +
    
    # Significant estimates, plot estimate and CI separately
    geom_text(data = fixed_df %>% filter(excludes_1 == TRUE), aes(x = statistic, y = response, label= ifelse(is.na(estimate), "", ci)), 
              # fontface = "bold",
              hjust = 0.5, vjust = vjust, size = size_text, family = P$font_family) +
    
    geom_text(data = fixed_df %>% filter(excludes_1 == TRUE), aes(x = statistic, y = response, label= ifelse(is.na(estimate), "", est)), 
              fontface = "bold",
              hjust = 0.5, vjust = vjust, size = size_text+1.1, family = P$font_family) +
    
    # Non-significant estimates, plot estimate and CI separately
    geom_text(data = fixed_df %>% filter(excludes_1 == FALSE), aes(x = statistic, y = response, label= ifelse(is.na(estimate), "", ci)), 
              # fontface = "bold",
              color = "grey70",
              hjust = 0.5, vjust = vjust, size = size_text, family = P$font_family) +
    
    geom_text(data = fixed_df %>% filter(excludes_1 == FALSE), aes(x = statistic, y = response, label= ifelse(is.na(estimate), "", est)), 
              fontface = "bold",
              color = "grey70",
              hjust = 0.5, vjust = vjust, size = size_text+1.1, family = P$font_family) +
    
    scale_color_gradient(name = "Adjusted odds ratio",
                         # low="blue", mid = "gray", high="red", na.value = "white", # scale_color_gradient2()
                         low="yellow", high="red", na.value = "white",
                         # midpoint = 1
    ) +
    scale_fill_gradient(name = "Adjusted odds ratio",
                        # low="blue", mid = "gray", high="red", na.value = "white", # scale_color_gradient2()
                        low="yellow", high="red", na.value = "white",
                        # midpoint = 1
    ) +
    guides(col = NULL, fill = guide_colourbar(title.position = "top",
                                              theme = theme(
                                                legend.key.width  = unit(width_legendbar, "lines"),
                                                legend.key.height = unit(height_legendbar, "lines")
                                              ))) +
    labs(y = ylab, x = "Predictor") +
    adapt_theme +
    scale_x_discrete(labels = function(x) sapply(x, wrap_equally, max_chars = 15), limits=rev, position = 'top',
                     expand = expansion(add = c(0,0))) +
    scale_y_discrete(labels = function(x) sapply(x, wrap_equally, max_chars = 15),
                     expand = expansion(add = c(0,0)))  
  pl_fixed
  
  if (name_dataset == "HILDA"){
    height =  200 
  } else if (name_dataset == "SHP"){
    height = 200 #  200
  }
  
  filepath_image = file.path(datalist$filepath_figs_dataset, sprintf("%s_heatmap_odds_ratio_selected_lead%d.pdf", name_dataset, chosen_lead))
  save_plot(pl_fixed, filepath_image, height = height)
  
}


# Create Supplementary Table
df_to_latex = function(name_dataset, chosen_lead, model_df, datalist, P){
  
  latex_df = model_df %>% dplyr::filter(subplot != "fixed") %>%
    dplyr::mutate(estimate = ifelse(is.na(estimate), "", ifelse(subplot == "sd_intercept",
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
    kbl(format = "latex",
        col.names = NULL,
        align = 'lccccc',
        caption = sprintf("Model fit and estimates of %s associations between events (Source: %s). The household intercept had to be dropped for models with no indicated household %s. $i$ = Individual; $h$ = Household; Marg. $R^2$ = Marginal explained variance; Cond. $R^2$ = Conditional explained variance; Adj. ICC = Adjusted Intra-class Correlation Coefficient.", 
                          ifelse(chosen_lead == 0, "contemporaneous", sprintf("lag-%d", chosen_lead)),
                          ifelse(name_dataset == "SHP", "Swiss Household Panel, SHP", "Household, Income and Labour Dynamics in Australia Survey, HILDA"),
                          # name_dataset, 
                          "$\\sigma$"
        ),
        label = sprintf("%s_odds_ratio_lead%d", name_dataset, chosen_lead),
        booktabs = TRUE, escape = FALSE) %>%
    add_header_above(c(" ", "$i$" = 1, "$h$" = 1, "Marg. $R^2$", "Cond. $R^2$", "Adj. ICC"), escape = FALSE) %>%
    add_header_above(c(" ", "Random $\\\\sigma$" = 2, "Variance Partitioning" = 3), escape = FALSE) %>%
    column_spec(1, "1.5in") %>% 
    column_spec(2:3, "1.1in") %>% 
    column_spec(4:6, "0.2in") 
  
}


# Helper function to get joint or conditional probabilities
get_prob = function(df_events, combo, 
                    chosen_lag = 0, 
                    lag_event = c("A", "B")[1], 
                    type = c("joint_prob", "cond_prob_given_eventA", "cond_prob_given_eventB")[1]){
  
  cond = foreach(x = combo,
                 .packages = c("dplyr"), .combine = "c"
  ) %do% {
    
    df_ = df_events %>% dplyr::select_at(c("p_id", unname(unlist(x))))
    
    if (chosen_lag > 0){
      # Lag event A
      if (lag_event == "A"){
        df_[["orig"]] <- df_[[x$eventA]] # Keep original for eventA == eventB
        df_[[x$eventA]] <- unsplit(
          lapply(split(df_[[x$eventA]], df_$p_id), function(x) c(rep(NA, chosen_lag), head(x, -chosen_lag))),
          df_$p_id
        )
      } else if (lag_event == "B"){
        # Lag event B
        df_[["orig"]] <- df_[[x$eventB]] # Keep original for eventA == eventB
        df_[[x$eventB]] <- unsplit(
          lapply(split(df_[[x$eventB]], df_$p_id), function(x) c(rep(NA, chosen_lag), head(x, -chosen_lag))),
          df_$p_id
        )
      }
    }
    
    # Remove NA from dataframe
    df_ = df_ %>%
      dplyr::select(-p_id) %>%
      dplyr::slice(which(complete.cases(.)))
    
    # Get base rate if looking at event in combination with itself
    if (x$eventA == x$eventB){
      if (chosen_lag == 0){
        m = mean(df_[[x$eventA]])
      } else if (chosen_lag > 0){
        # Co-occurrence of event with itself at lag
        m = sum(rowSums(df_) == 2) / nrow(df_)
      }
    } else {
      
      # Remove original unlagged event that was kept for eventA == eventB
      if (chosen_lag > 0){
        df_$orig = NULL
      }
      
      if (type == "joint_prob"){
        # Divide number of co-occurrences by total number of possible co-occurrences
        m = sum(rowSums(df_) == 2) / nrow(df_)
      } else if (type == "cond_prob_given_eventA"){
        # Divide number of co-occurrences by number of times event B occurred
        m = mean(rowSums(df_) == 2) / mean(df_[[x$eventA]])
      } else if (type == "cond_prob_given_eventB"){
        # Divide number of co-occurrences by number of times event B occurred
        m = mean(rowSums(df_) == 2) / mean(df_[[x$eventB]])
      }
    }
    rm(df_)
    return(m)
  }
  return(cond)
}

# Compute joint or conditional probability
run_joint_cond_prob = function(name_dataset, filepath_base, P, 
                               chosen_lags = c(0,1), rerun = F){
  
  # Prepare data
  datalist = prepare_data(name_dataset, filepath_base, P$event_dicts, rerun = rerun)
  
  # Reformat dataframe of events
  df_events = datalist$df_per_event_pp_py %>% 
    filter(valence == "negative") %>%
    dplyr::mutate(occur = ifelse(nr_missing_occur > 0, NA, nr_occur)) %>%
    dplyr::select(p_id, wave_nr, event, occur) %>% 
    tidyr::spread(event, occur) %>%
    tidyr::complete(p_id, wave_nr) %>%
    arrange(p_id, wave_nr) 
  df_events %>% head
  
  # Event names
  negevent_names = setdiff(colnames(df_events), c("p_id", "wave_nr"))
  
  # All event combinations
  combo = expand.grid(eventA = negevent_names,
                      eventB = negevent_names
  ) %>%
    dplyr::mutate_at(c("eventA", "eventB"), ~as.character(.)) %>%
    # Create list entry per row
    split(seq(nrow(.)))
  length(combo)
  
  prob_df = lapply(chosen_lags, function(chosen_lag){
    cbind(lag_eventA = chosen_lag, 
          combo %>% do.call(rbind, .),
          cond_prob_given_eventB = get_prob(df_events, combo, chosen_lag, lag_event = "A", type = "cond_prob_given_eventB"),
          joint_prob = get_prob(df_events, combo, chosen_lag, lag_event = "A", type = "joint_prob"))
  }) %>% do.call(rbind, .) %>% as.data.frame()
  
  return(list(datalist = datalist,
              df_events = df_events,
              prob_df = prob_df))
  
}

# Plot joint or conditional probability heatmaps
plot_joint_cond = function(name_dataset, datalist, prob_df, chosen_lag, event_dict, lag_event = c("A", "B")[1], type = c("joint_prob", "cond_prob_given_eventA", "cond_prob_given_eventB")[1], size_text = 9, size_label = 2.8, size_point = 8){
  
  # Pick correct lag and estimate
  plot_df = prob_df %>%
    dplyr::filter(!!sym(sprintf("lag_event%s", lag_event)) == !!chosen_lag) %>%
    dplyr::mutate(estimate = !!sym(type)) 
  
  # Set estimate to NA for events that cannot co-occur
  if (name_dataset == "SHP" & chosen_lag == 0){
    illness_accident_events = event_dict %>% dplyr::filter(grepl("PL01R", var_name)) %>% dplyr::pull(recode_var)
    plot_df = plot_df %>%
      dplyr::mutate(estimate = ifelse(eventA %in% illness_accident_events & eventB %in% illness_accident_events & eventA != eventB, NA, estimate))
  }
  
  # Get plot labels
  if (chosen_lag == 0){
    xlab = "Event B"
    ylab = "Event A"
    x = "eventB"
    y = "eventA"
  } else if (chosen_lag > 0){
    if (lag_event == "A"){
      ylab = sprintf("Event A (%d year%s later)", chosen_lag, ifelse(chosen_lag > 1, "s", ""))
      xlab = "Event B"
      x = "eventB"
      y = "eventA"
    } else if (lag_event == "B"){
      ylab = sprintf("Event B (%d year%s later)", chosen_lag, ifelse(chosen_lag > 1, "s", ""))
      xlab = "Event A"
      x = "eventA"
      y = "eventB"
    }
  }
  
  if (type == "joint_prob"){
    fill_name = "Joint probability"
    transform = "log10"
    plot_df$label = log10(plot_df$estimate)
  } else if (type == "cond_prob_given_eventA"){
    fill_name = "Conditional probability of event B given event A"
    transform = "identity"
    plot_df$label = plot_df$estimate
  } else if (type == "cond_prob_given_eventB"){
    fill_name = "Conditional probability of event A given event B"
    transform = "identity"
    plot_df$label = plot_df$estimate
  }
  
  # Rename events with event description
  plot_df = plot_df %>%
    mutate_at(c("eventA", "eventB"), 
              ~ dplyr::recode_factor(.x, !!!tibble::deframe(P$event_dicts[[name_dataset]] %>%
                                                              dplyr::select(recode_var, description) %>%
                                                              # Order events in event_dict according to frequency
                                                              slice(match(datalist$df_per_event %>% arrange(nr_occur) %>% pull(event) %>% as.character(), recode_var))), .ordered = T)) 
  
  # Adapt theme
  adapt_theme = P$own_theme +
    theme(legend.position = 'bottom') +
    theme(axis.text.x=element_text(size=size_text+2, angle=90,hjust = 0,vjust = 1),
          axis.text.y=element_text(size=size_text+2),
          legend.text=element_text(size=size_text),
          legend.title=element_text(size=size_text+4),
          axis.title=element_text(size=size_text+4),
          plot.margin = unit(c(.01, .15, .01, .01), 'cm'),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.spacing = unit(0, "cm"),    # Remove space between legend items
          legend.margin = margin(0, 0, 0, 0)
    )
  
  
  width_legendbar = 10
  height_legendbar = 1
  pl_prob = plot_df %>%
    ggplot() + 
    geom_point(aes(x = .data[[x]], y = .data[[y]],
                   col = estimate,
    ), shape = 19, size = size_point) +
    geom_text(aes(x = .data[[x]], y = .data[[y]],
                  label= ifelse(is.na(estimate), "", sprintf("%.3f", as.numeric(estimate)))),
              size = size_label, family = P$font_family) +
    
    
    scale_color_gradient(name = fill_name,
                         low="yellow", high="red", na.value = "white"
    ) +
    guides(color = guide_colourbar(title.position = "top",
                                   theme = theme(
                                     legend.key.width  = unit(width_legendbar, "lines"),
                                     legend.key.height = unit(height_legendbar, "lines")
                                   ))) +
    labs(y = ylab, x = xlab) +
    adapt_theme +
    scale_x_discrete(labels = function(x) sapply(x, wrap_equally, max_chars = 31), limits=rev, position = 'top',
                     expand = expansion(add = c(0.6,0.6))) +
    
    scale_y_discrete(labels = function(x) sapply(x, wrap_equally, max_chars = 31), expand = expansion(add = c(0.6,0.6)))  
  pl_prob
  
  if (name_dataset == "HILDA"){
    height =  190 
  } else if (name_dataset == "SHP"){
    height = 190
  }
  
  
  filepath_image = file.path(datalist$filepath_figs_dataset, sprintf("%s_%s_lag%d.pdf", name_dataset, type, chosen_lag))
  save_plot(pl_prob, filepath_image, height = height)
  
}


