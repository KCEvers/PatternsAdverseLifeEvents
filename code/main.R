##### COMPLETE ANALYSIS OF ADVERSE LIFE EVENTS #####
# Get base filepath
filepath_base = dirname(dirname(rstudioapi::getActiveDocumentContext()$path))
stopifnot("Empty base filepath!" = nzchar(filepath_base))

# Load functions
source(file.path(filepath_base, "code/setup.R"))
source(file.path(filepath_base, "code/prepare_data_func.R"))
source(file.path(filepath_base, "code/cooccurrence_event_types.R"))
source(file.path(filepath_base, "code/lagged_event_counts.R"))
source(file.path(filepath_base, "code/accumulation_event_counts.R"))

# Get parameters
P = get_parameters()

##### CO-OCCURRENCE OF EVENT TYPES #####

# Set up cluster
nr_cl = ifelse(parallel::detectCores() > 13, 13, parallel::detectCores() - 1)
cl <- parallel::makeCluster(nr_cl)
doParallel::registerDoParallel(cl)

# Run co-occurrence analysis for different leads for both SHP and HILDA
chosen_leads = c(0, 1)
SHP = run_cooccurrence("SHP", filepath_base, P, chosen_leads = chosen_leads)
HILDA = run_cooccurrence("HILDA", filepath_base, P, chosen_leads = chosen_leads)
parallel::stopCluster(cl)

# Plot co-occurrence in heatmap
for (chosen_lead in chosen_leads){
  plot_cooccur("SHP", SHP$datalist, P, SHP$conv_model_est[[as.character(chosen_lead)]], chosen_lead)
  plot_cooccur("HILDA", HILDA$datalist, P, HILDA$conv_model_est[[as.character(chosen_lead)]], chosen_lead)
}

# Create Supplementary Table
df_to_latex("SHP", 0, SHP$conv_model_est$`0`, SHP$datalist, P)
df_to_latex("HILDA", 0, HILDA$conv_model_est$`0`, HILDA$datalist, P)
df_to_latex("SHP", 1, SHP$conv_model_est$`1`, SHP$datalist, P)
df_to_latex("HILDA", 1, HILDA$conv_model_est$`1`, HILDA$datalist, P)

### Compute joint and conditional probabilities
chosen_lags = c(0,1)
SHP = run_joint_cond_prob("SHP", filepath_base, P,
                          chosen_lags = chosen_lags)
HILDA = run_joint_cond_prob("HILDA", filepath_base, P,
                            chosen_lags = chosen_lags)

# Plot joint and conditional probabilities in heatmap
for (chosen_lag in chosen_lags){
  plot_joint_cond("SHP", SHP$datalist, SHP$prob_df, chosen_lag, P$event_dicts[["SHP"]], lag_event = "A", type = "joint_prob",
                  size_point = 9.5, size_text = 9, size_label = 2.5)

  plot_joint_cond("SHP", SHP$datalist, SHP$prob_df, chosen_lag, P$event_dicts[["SHP"]], lag_event = "A", type = "cond_prob_given_eventB",
                  size_point = 9.5, size_text = 9, size_label = 2.5)

  plot_joint_cond("HILDA", HILDA$datalist, HILDA$prob_df, chosen_lag, P$event_dicts[["HILDA"]], lag_event = "A", type = "joint_prob",
                  size_point = 8.5, size_text = 8, size_label = 2.2)

  plot_joint_cond("HILDA", HILDA$datalist, HILDA$prob_df, chosen_lag, P$event_dicts[["HILDA"]], lag_event = "A", type = "cond_prob_given_eventB",
                  size_point = 8.5, size_text = 8, size_label = 2.2)

}



##### LAGGED EVENT COUNTS #####

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
lagged_table(summlist_SHP, summlist_HILDA)


##### ACCUMULATION EVENT COUNTS #####
# Fit models
SHP = run_accumulation("SHP", filepath_base, P)
HILDA = run_accumulation("HILDA", filepath_base, P)

# Polya fits better
SHP$df_comp_AIC
SHP$df_comp_BIC

HILDA$df_comp_AIC
HILDA$df_comp_BIC

# Create table with estimates
accumulation_table(SHP, HILDA)


### Simulate cumulative count data for SHP and HILDA
n_sim = 1000
sim_types = c("Empirical", "Poisson", "Frailty", "Polya urn")
df_SHP = simulate_data(SHP, n_sim, sim_types) %>%
  dplyr::mutate(dataset = "SHP")

df_HILDA = simulate_data(HILDA, n_sim, sim_types) %>%
  dplyr::mutate(dataset = "HILDA")


# Prepare colors
col_values = c(P$col_data, P$col_poisson, P$col_frailty, P$col_polya) %>%
  setNames(sim_types)
fill_values =  c(P$fill_data, P$fill_poisson, P$fill_frailty, P$fill_polya) %>%
  setNames(sim_types)

# Plot cumulative distributions of empirical and simulated data
pl_distr_SHP = plot_cum_distr("SHP", df_SHP, P, col_values, fill_values)

pl_distr_HILDA = plot_cum_distr("HILDA", df_HILDA, P, col_values, fill_values)

# Create combined plot
pl_combo = ((pl_distr_SHP + theme(legend.margin = margin(t = 10, r = 20, b = 10, l = 250, unit = "pt")) + (pl_distr_HILDA + theme(strip.text.y.left = element_blank()) + theme(legend.margin = margin(t = 10, r = 20, b = 10, l = 250, unit = "pt"))) +
               plot_layout(axis_titles = "collect")) + guide_area() +
              plot_layout(guides = 'collect') + plot_layout(heights = c(1, .1)) )
pl_combo

filepath = file.path(filepath_base, "figs", "cumulative_nr_events.pdf")
save_plot(pl_combo, filepath, height = 170)

# Create GIF

### Heavy tails
SHP_tails = run_heavy_tails("SHP", SHP$df_accum, SHP$datalist)
HILDA_tails = run_heavy_tails("HILDA", HILDA$df_accum, HILDA$datalist)

# Create table with estimates
heavy_tails_table(SHP_tails$df_est, HILDA_tails$df_est)

# Plot heavy tails
pl1 = SHP_tails$pl
pl2 = HILDA_tails$pl

pl = ((pl1 + theme(legend.margin = margin(t = 10, r = 20, b = 10, l = 225, unit = "pt")) +
         (pl2 + theme(legend.margin = margin(t = 10, r = 20, b = 10, l = 225, unit = "pt"))) +
         plot_layout(axis_titles = "collect")) + guide_area() +
        plot_layout(guides = 'collect') + plot_layout(heights = c(1, .1)) )


filepath_image = file.path(filepath_base, "figs", sprintf("heavy_tails.pdf"))
save_plot(pl, filepath_image, height = 100)
