##### SET-UP AND HELPER FUNCTIONS TO ANALYSE ADVERSE LIFE EVENTS #####

# Load libraries
library(dplyr)
library(lme4)
library(glmmTMB)
library(DHARMa) # Model diagnostics
library(performance) # Model diagnostics
library(foreach)
library(ggplot2)
library(patchwork)
library(knitr) # Tables
library(kableExtra) # Tables

# Power-law
if (!requireNamespace("poweRlaw", quietly = TRUE)) {
  install.packages("poweRlaw")
}
library(poweRlaw)

# Fonts
library(showtext)
font_add_google("Nunito Sans", "nunito")
showtext::showtext_auto()

### Functions
get_parameters = function(){

  P = list()

  # Colors
  P$col_palette = "grDevices::Purple-Yellow"
  P$cols = paletteer::paletteer_c(P$col_palette, 22, direction = -1)
  P$cols_darker = ggsci::pal_npg("nrc", alpha = .5)(10)
  P$cols = ggsci::pal_npg("nrc", alpha = .35)(10)
  P$fills = ggsci::pal_npg("nrc", alpha = .1)(10)
  P$facet_col = P$cols[7]
  P$col_facets = "grey30"
  P$fill_facets = P$cols[7]
  P$col_data = P$cols[8]
  P$col_poisson = P$cols[2]
  P$col_frailty = P$cols[3]
  P$col_polya = P$cols[4]
  P$fill_data = P$fills[8]
  P$fill_poisson = P$fills[2]
  P$fill_frailty = P$fills[3]
  P$fill_polya = P$fills[4]
  P$font_family = "nunito"
  P$dpi = 500
  P$width = 180

  P$own_theme = ggpubr::theme_pubr(base_family = P$font_family, base_size = 16, border = F) +
              theme(
                plot.margin = unit(c(.5,.5,.5,.5), 'cm'),
                    plot.title = element_text(hjust = .5, size = 24),
                    plot.subtitle = element_text(hjust = .5, size = 20),
                strip.text = element_text(size = 14, color = "black"),
                strip.text.y.left = element_text(angle = 0),
                axis.text = element_text(size = 16),
                axis.title = element_text(size = 20),
                legend.title = element_text(size = 20, hjust = .5),
                legend.text = element_text(size = 16),
                text = element_text(),
                title = element_text(size = 22),
                    strip.background=element_rect(
                                                  fill=P$fill_facets))

  # Dictionary of event names and descriptions for SHP and HILDA
  P$event_dicts = list("SHP" =
                       matrix(c(
                         c("PL06", "Illness or accident of closely related person", "illness_accident_close_person", "negative", "independent"),
                         c("PL11", "Death of closely related person","death_close_person", "negative", "independent"),
                         c("PL16", "Termination of close relationship","end_relationship", "negative", "dependent"),
                         c("PL21", "Conflicts with or among related persons","conflicts_close_person", "negative", "dependent"),
                         c("PL26", "Problems with own children","problems_with_children", "negative", "dependent"),
                         c("PL01R_1", "Physical illness", "physical_illness", "negative", "independent"),
                         c("PL01R_2", "Mental illness", "mental_illness", "negative", "dependent"),
                         c("PL01R_3", "Work accident", "work_accident", "negative", "independent"),
                         c("PL01R_4", "Road accident", "road_accident", "negative", "independent"),
                         c("PL01R_5", "Accident in home or garden", "home_accident", "negative", "independent"),
                         c("PL01R_6", "Sport accident", "sport_accident", "negative", "independent"),
                         c("PL01R_7", "Other or unspecified illness or accident", "unspecified_illness_accident", "negative", "ambiguous"),

                         c("PL36_1", "Dismissal or unemployment", "dismissal_unemployment", "negative", "dependent"),
                         c("PL36_2", "Retirement", "retirement", "neutral", "dependent"),
                         c("PL36_3", "Changed job", "changed_job", "neutral", "ambiguous"),
                         c("PL36_4", "Problems and conflicts at work", "problems_at_work", "negative", "dependent"),
                         c("PL36_5", "Problems at school or academic failure", "problems_at_school", "negative", "dependent"),
                         c("PL36_6", "Change of apartment or place of residence", "changed_residence", "neutral", "ambiguous"),
                         c("PL36_7", "Problems with neighbours", "problems_with_neighbours", "negative", "dependent"),
                         c("PL36_8", "Environmental pollution", "environmental_pollution", "negative", "independent"),
                         c("PL36_9", "Financial difficulties", "financial_difficulties", "negative", "dependent"),
                         c("PL36_10", "Problems with insurances", "problems_with_insurances", "negative", "dependent"),
                         c("PL36_11", "Problems with the judiciary", "problems_with_judiciary", "negative", "dependent"),
                         c("PL36_12", "Material damage", "material_damage", "negative", "independent"),
                         c("PL36_13", "Psychological trauma", "psychological_trauma", "negative", "ambiguous"),
                         c("PL36_14", "Hospitalization or operation", "hospitalization_operation", "negative", "independent"),
                         c("PL36_15", "Birth", "birth", "positive", "dependent"),
                         c("PL36_16", "Other unspecified life event", "unspecified_event", "neutral", "ambiguous"),
                         c("PL36_17", "Burglary", "burglary", "negative", "independent"),
                         c("PL36_18", "Insult or threat", "insult_threat", "negative", "independent"),
                         c("PL36_19", "Hit or wounded", "hit_wounded", "negative", "independent"),
                         c("PL36_20", "Obliged to change jobs because of closure of company", "obliged_to_change_jobs_closure_company", "negative", "independent"),
                         c("PL36_21", "Obliged to change jobs because of employer (cuts in manpower, redudancy)", "obliged_to_change_jobs_due_to_employer", "negative", "independent")
                       ), ncol = 5, byrow = T ) %>%
                       magrittr::set_colnames(c("var_name", "description", "recode_var", "valence", "in_dependent")) %>% as.data.frame(),
                     "HILDA" =
                       matrix(c(
                       c("lebth", "Birth/adoption of new child", "child_born", "positive", "dependent", "family", 39), # pregnancy score
                       c("ledfr", "Death of a close friend", "death_friend", "negative", "independent", "relationships", 37),
                       c("ledhm", "A weather related disaster damaged or destroyed your home", "natural_disaster", "negative", "independent", "non-interpersonal", NA),
                       c("ledrl", "Death of close relative/family member", "death_family", "negative", "independent", "family", 63),
                       c("ledsc", "Death of spouse or child", "death_spouse_or_child", "negative", "independent", "family", 100),
                       c("lefni", "Major improvement in finances", "improvement_finances", "positive", "dependent", "finances", 38),
                       c("lefnw", "Major worsening in finances", "worsening_finances", "negative", "dependent", "finances", 38),
                       c("lefrd", "Fired or made redundant", "fired", "negative", "dependent", "work", 47),
                       c("leinf", "Serious injury/illness to family member", "illness_injury_family", "negative", "independent", "family", 44),
                       c("leins", "Serious personal injury/illness", "injury_illness_self", "negative", "independent", "health", 53),
                       c("lejlf", "Close family member detained in jail", "family_jailed", "negative", "independent", "family", NA),
                       c("lejls", "Detained in jail", "jailed", "negative", "dependent", "legal", 63),
                       c("lejob", "Changed jobs", "changed_job", "neutral", "ambiguous", "work", 36),
                       c("lemar", "Got married", "married", "positive", "dependent", "family", 50),
                       c("lemvd", "Changed residence", "moved", "neutral", "ambiguous", "non-interpersonal", 20),
                       c("lepcm", "Victim of a property crime", "victim_property_theft", "negative", "independent", "non-interpersonal", NA),
                       c("leprg", "Pregnancy", "pregnant", "positive", "dependent", "family", 39),
                       c("leprm", "Promoted at work", "promoted", "positive", "dependent", "work", 29),
                       c("lercl", "Got back together with spouse", "reconciled_with_spouse", "positive", "dependent", "family", 45),
                       c("lertr", "Retired from the workforce", "retired", "neutral", "ambiguous", "work", 45),
                       c("lesep", "Separated from spouse", "separated_from_spouse", "negative", "dependent", "family", 65),
                       c("levio", "Victim of physical violence", "victim_physical_violence", "negative", "independent", "non-interpersonal", NA)
                     ), ncol = 7, byrow = T ) %>%
                       magrittr::set_colnames(c("var_name", "description", "recode_var", "valence", "in_dependent", "domain", "HolmesRahe")) %>%
                       as.data.frame())
  return(P)
}

# Recode event codes with descriptions, optionally order by frequency
recode_events = function(df, recode_col, event_dict, df_per_event = NULL){

  # Don't order events by frequency
  if (is.null(df_per_event)){
    df = df %>%
    dplyr::mutate_at(recode_col,
            ~ dplyr::recode_factor(.x, !!!tibble::deframe(event_dict %>%
                                                            dplyr::select(recode_var, description)),
                                   .ordered = T))
  } else {
    # Order events by frequency
    df = df %>%
      dplyr::mutate_at(recode_col,
                ~ dplyr::recode_factor(.x, !!!tibble::deframe(event_dict %>%
                                                                dplyr::select(recode_var, description) %>%
                                                                # Order events in event_dict according to frequency
                                                                dplyr::slice(match(df_per_event %>% dplyr::arrange(nr_occur) %>% dplyr::pull(event) %>% as.character(), recode_var))), .ordered = T))
  }
  return(df)
}

# Summarize linear mixed-effects model
summarise_model = function(model, run_dharma = FALSE, run_diagnose = FALSE){

  out = list()

  # Fixed effect coefficients
  out[["fixef_coef"]] <- broom.mixed::tidy(model, effects = "fixed")
  out[["check.std.error"]] = all(!is.na(out[["fixef_coef"]]$std.error))

  # Number of observations
  out[["nobs"]] = summary(model)$nobs
  out[["ngrps"]] = summary(model)$ngrps
  out[["vcov"]] = summary(model)$vcov
  out[["varcor"]] = summary(model)$varcor

  # # Missing values
  # keep_p_id = unique(row.names(ranef(model_list[[1]])$p_id))
  # print(length(keep_p_id))

  # Wald CI
  out[["Wald_CI"]] = confint(model, method = "Wald")

  # Multicollinearity
  if (nrow(out[["fixef_coef"]]) > 1){
    out[["multicollinearity"]] = performance::check_collinearity(model, component = "all")
  }

  # If the model contains random effects
  if (insight::is_mixed_model(model)){

    # Intraclass correlation and marginal and conditional explained variance
    out[["ICC"]] = performance::icc(model)  # Measures how much variance is explained by individual effects
    # Unadjusted (Marginal) ICC: This tells you how much variance is explained by the random effects without considering fixed effects.
    # Adjusted (Conditional) ICC: This tells you how much variance is still explained by the random effects after including fixed effects.
    out[["R2"]] = performance::r2(model)
  }

  if (insight::is_mixed_model(model)){
    # If random effect variances are very small (close to zero), the model may struggle to estimate them properly
    out[["VarCorr"]] = VarCorr(model)

    if (run_diagnose){
      out[["diagnose"]] = diagnose(model, check_hessian = F) # returns a logical value based on whether anything questionable was found
    }
  }

  # Check singularity (should be FALSE) and convergence (should be TRUE)
  out[["singularity"]] = performance::check_singularity(model)
  out[["check_convergence"]] = performance::check_convergence(model) #model$fit$convergence
  out[["convergence"]] = model$fit$convergence
  out[["message"]] = model$fit$message
  out[["sdr"]] = model$sdr
  out[["pdHess"]] = model$sdr$pdHess

  # Profile confidence intervals
  out[["profile_CI"]] = confint(model)

  # # Check with performance and dharma package
  # # https://easystats.github.io/performance/articles/simulate_residuals.html
  if (run_dharma){
    sim_res <- performance::simulate_residuals(model)
    out[["sim_res"]] = sim_res
    out[["check_residuals"]] = performance::check_residuals(sim_res)
    out[["check_overdispersion"]] = performance::check_overdispersion(sim_res)
    out[["check_zeroinflation"]] = performance::check_zeroinflation(sim_res)
    out[["check_outliers"]] = performance::check_outliers(sim_res, type = "bootstrap")
    # out[["sim_res"]] = performance::check_model(sim_res, size_dot = 1.5) # Takes very long
  }

  out[["df_est"]] = broom.mixed::tidy(model) %>%
    # Term column in Wald CI has different names, keep to check match
    cbind(out[["profile_CI"]] %>% as.data.frame() %>%
            tibble::rownames_to_column("term") %>% dplyr::select(-any_of(c("term", "Estimate")) ))

  return(out)
}

# Save plot
save_plot = function(pl, filepath_image, height, width = 180, dpi = 500, units = "mm"){
  ggsave(
    filename = filepath_image,  # File format: PNG (or TIFF for high-quality print)
    plot = pl,                 # The ggplot object
    dpi = dpi,                # Resolution in dots per inch
    width = width,       # Convert mm to inches (180 mm ÷ 25.4)
    height = height,      # Adjust height accordingly (e.g., 120 mm)
    units = units              # Use inches for precise control
  )
}

# Wrap labels across multiple lines, ensuring the number of characters is approximately equal across lines
wrap_equally <- function(text, max_chars = 20) {
  if (nchar(text) <= max_chars) {
    return(text)  # No wrapping needed
  }

  words <- unlist(strsplit(text, " "))  # Split text into words
  n <- length(words)

  if (n <= 2) {
    return(text)  # Avoid unnecessary wrapping for short text
  }

  # Compute cumulative character count per word
  char_counts <- cumsum(nchar(words) + 1)  # +1 accounts for spaces

  # Find the best split point where total characters are nearly equal
  total_chars <- sum(nchar(words)) + (n - 1)  # Total characters including spaces
  split_point <- which.min(abs(char_counts - total_chars / 2))  # Closest to half

  # Create two balanced lines
  line1 <- paste(words[1:split_point], collapse = " ")
  line2 <- paste(words[(split_point + 1):n], collapse = " ")

  return(paste(line1, "\n", line2, sep = ""))
}
