##### COMPLETE ANALYSIS OF ADVERSE LIFE EVENTS #####
# Get base filepath
filepath_base = dirname(dirname(rstudioapi::getActiveDocumentContext()$path))
stopifnot("Empty base filepath!" = nzchar(filepath_base))

# Load functions
source(file.path(filepath_base, "code/setup.R"))
source(file.path(filepath_base, "code/prepare_data_func.R"))

# Set up cluster
nr_cl = ifelse(parallel::detectCores() > 13, 10, parallel::detectCores() - 1)
cl <- parallel::makeCluster(nr_cl)
doParallel::registerDoParallel(cl)

# Get parameters
P = get_parameters()
rerun = FALSE#TRUE
name_dataset = "SHP"
events_dict = P$event_dicts[name_dataset]
chosen_leads = c(0,1)
chosen_lead = chosen_leads[2]

##### DESCRIPTIVES #####
source(file.path(filepath_base, "code/descriptives.R"))

results_SHP <- run_descriptives("SHP", filepath_base, P, rerun = rerun)
results_HILDA <- run_descriptives("HILDA", filepath_base, P, rerun = rerun)
plot_combined_event_frequency(results_SHP, results_HILDA, filepath_base)
plot_combined_event_over_time(results_SHP, results_HILDA, filepath_base)
plot_combined_observation_years(results_SHP, results_HILDA, filepath_base)

##### CO-OCCURRENCE OF EVENT TYPES #####
source(file.path(filepath_base, "code/cooccurrence_event_types.R"))

# Run co-occurrence analysis for different leads for both SHP and HILDA
chosen_leads = c(0, 1)
SHP = run_cooccurrence("SHP", filepath_base, P, chosen_leads = chosen_leads,
                       rerun = rerun)
HILDA = run_cooccurrence("HILDA", filepath_base, P, chosen_leads = chosen_leads,
                         rerun = rerun)
parallel::stopCluster(cl)

# Plot co-occurrence in heatmap
for (chosen_lead in chosen_leads){
  plot_cooccur("SHP", SHP$datalist, P,
               SHP$model_dfs[[as.character(chosen_lead)]], chosen_lead)
  plot_cooccur("HILDA", HILDA$datalist, P,
               HILDA$model_dfs[[as.character(chosen_lead)]], chosen_lead)

  # Create Supplementary Table
  df_to_latex_cooccur("SHP", chosen_lead,
                      SHP$model_dfs[[as.character(chosen_lead)]], SHP$datalist, P)
  df_to_latex_cooccur("HILDA", chosen_lead,
                      HILDA$model_dfs[[as.character(chosen_lead)]], HILDA$datalist, P)
}

### Compute joint and conditional probabilities
chosen_lags = c(0, 1)
SHP = run_joint_cond_prob("SHP", filepath_base, P, chosen_lags = chosen_lags)
HILDA = run_joint_cond_prob("HILDA", filepath_base, P, chosen_lags = chosen_lags)

# Plot joint and conditional probabilities in heatmap
for (chosen_lag in chosen_lags){
  for (type in c("joint_prob", "cond_prob_given_eventB")){
    plot_joint_cond("SHP", SHP$datalist, SHP$prob_df, chosen_lag,
                    P$event_dicts[["SHP"]], P, lag_event = "A", type = type,
                    size_point = 9.5, size_text = 9, size_label = 2.5)

    plot_joint_cond("HILDA", HILDA$datalist, HILDA$prob_df, chosen_lag,
                    P$event_dicts[["HILDA"]], P, lag_event = "A", type = type,
                    size_point = 8.5, size_text = 8, size_label = 2.2)
  }
}



##### LAGGED EVENT COUNTS #####
source(file.path(filepath_base, "code/lagged_event_counts.R"))

# Run lagged analysis
SHP = run_lagged("SHP", filepath_base, P)
HILDA = run_lagged("HILDA", filepath_base, P)

# Model convergence
print(SHP$df_convergence)
print(HILDA$df_convergence)

# Best-fitting model
print(SHP$df_comp_AIC)
print(SHP$df_comp_BIC)

print(HILDA$df_comp_AIC)
print(HILDA$df_comp_BIC)

# Choose model: model_random_int_slope fits best for HILDA but model_random_int_slope doesn't converge for SHP, for consistency use model_random_int for both
model_SHP = SHP$model_random_int
model_HILDA = HILDA$model_random_int

summary(model_SHP)$nobs
summary(model_SHP)$ngrps$cond

summary(model_HILDA)$nobs
summary(model_HILDA)$ngrps$cond

# Summarise and check model
summlist_SHP = summarise_model(model_SHP)
summlist_HILDA = summarise_model(model_HILDA)

# Create table with estimates
create_lagged_table(summlist_SHP, summlist_HILDA)


##### ACCUMULATION EVENT COUNTS #####
source(file.path(filepath_base, "code/accumulation_event_counts.R"))

# Fit models
for (nr_years in c(10, 15, 20)){ # Assess for multiple years for robustness
  SHP = run_accumulation("SHP", filepath_base, P, nr_years = nr_years, rerun = rerun)
  HILDA = run_accumulation("HILDA", filepath_base, P, nr_years = nr_years, rerun = rerun)

  # Polya fits better
  SHP$df_comp_AIC
  SHP$df_comp_BIC

  HILDA$df_comp_AIC
  HILDA$df_comp_BIC

  # Create table with estimates
  create_accumulation_table(SHP, HILDA, nr_years)
}


### Simulate cumulative count data for SHP and HILDA
n_sim = 1000
df_SHP = simulate_accumulation_data(SHP, n_sim, P$sim_types) %>%
  dplyr::mutate(dataset = "SHP")

df_HILDA = simulate_accumulation_data(HILDA, n_sim, P$sim_types) %>%
  dplyr::mutate(dataset = "HILDA")


# Plot cumulative distributions of empirical and simulated data
pl_distr_SHP = plot_cum_distr("SHP", df_SHP, P, P$col_values_accum, P$fill_values_accum)
pl_distr_HILDA = plot_cum_distr("HILDA", df_HILDA, P, P$col_values_accum, P$fill_values_accum)

pl_range_SHP = plot_cum_distr_range("SHP", df_SHP, P, P$col_values_accum, P$fill_values_accum)
pl_range_HILDA = plot_cum_distr_range("HILDA", df_HILDA, P, P$col_values_accum, P$fill_values_accum)

# Create combined plot
plot_cum_distr_combo(pl_distr_SHP, pl_distr_HILDA, filepath_base)
plot_cum_range_combo(pl_range_SHP, pl_range_HILDA, filepath_base)

# Create GIF

### Heavy tails
SHP_tails = run_heavy_tails("SHP", SHP$df_accum, SHP$datalist)
HILDA_tails = run_heavy_tails("HILDA", HILDA$df_accum, HILDA$datalist)

# Create table with estimates
create_heavy_tails_table(SHP_tails$df_est, HILDA_tails$df_est)

# Plot heavy tails
plot_heavy_tails(SHP_tails$pl, HILDA_tails$pl, filepath_base)
