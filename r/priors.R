# nolint start: object_usage_linter

source("r/multi_level_hierarchy.R")

get_tumor_priors <- function(stan_data, coef_elicited_priors) {
  # Directly specified priors (simplified)
  # Choose log-total rate prior similar to historical center; adjust if needed.
  # Updated: shifted mean from -1.0 to -2.0 to reduce prior-posterior conflict
  tr_loc_pop_mean <- -2.0 # centered closer to typical posterior (~-2.6)
  tr_loc_pop_sd <- 0.8 # tightened from 1.2 to reduce total variance with hierarchies
  # Fraction (logit) prior - updated to reduce severe prior-posterior conflict
  # Old prior at -1.0 (~27% fraction) conflicted with posterior at +2.4 (~92%)
  frac_logit_loc_pop_mean <- 1.5 # shifted from -1.0 to reduce 4+ SD conflict
  frac_logit_loc_pop_sd <- 1.0 # slightly wider to allow data to inform
  # Initial proportion logit
  init_logit_loc_pop_mean <- 0.0 # formerly pop_decrease_prop_logis_mean
  init_logit_loc_pop_sd <- 0.8 # tightened from 1.5 (was extremely wide!)

  # Get n_levels from stan_data (default 2 for backward compat)
  n_levels <- stan_data$n_levels %||% 2L
  n_covar <- stan_data$n_covar

  lst(
    # GP hyperparameters
    pop_tumor_gp_rho_meanlog = 3,
    pop_tumor_gp_rho_sdlog = 0.6,
    log_patient_tumor_gp_rho_sd_sd = 1.75,
    pop_decrease_process_sd_sd = 0.1,
    pop_growth_process_sd_sd = 0.1,
    process_corr_param = 2.0,
    # inv_gamma prior for measure_sd keeps mass away from zero
    # mode = beta/(alpha+1) = 0.75/6 = 0.125 (at typical posterior)
    measure_sd_alpha = 5,
    measure_sd_beta = 0.75,

    # Process parameters
    decrease_process_alpha = 9.7,
    decrease_process_beta = 38.4,
    growth_process_alpha = 9.7,
    growth_process_beta = 38.4,

    # Total rate module hyperparams (multi-level): c(trial, patient)
    tr_loc_pop_mean = tr_loc_pop_mean,
    tr_loc_pop_sd = tr_loc_pop_sd,
    tr_sd_level_intercept_sd = c(0.35, 0.35),
    tr_coef_qr_pop_mean = as.array(rep(0, n_covar)),
    tr_coef_qr_pop_sd = as.array(rep(1, n_covar)),
    tr_sd_level_slope_sd = list(
      trial = rep(0.15, n_covar),
      patient = rep(0.10, n_covar)
    ),

    # Patient-level process noise (AR(1) time-varying rates per patient)
    # Note: These are deviations in log-rates (decrease/growth), which integrate over time
    # Even small rate deviations accumulate into substantial tumor trajectory effects
    tr_log_sd_pop_process_noise_mean      = -3,    # log(0.05) ≈ -3, median σ ≈ 0.05 (5% deviations)
    tr_log_sd_pop_process_noise_sd        = 0.5,   # allows σ ~[0.02, 0.13] (95% CI)
    tr_logit_phi_pop_process_noise_mean   = 1.4,   # logit(0.8) ≈ 1.39, favors high correlation
    tr_logit_phi_pop_process_noise_sd     = 0.5,   # allows φ ~[0.6, 0.9] (95% CI)
    tr_log_sd_patient_process_noise_sd    = 0.3,   # patient-level variation in log(σ)
    tr_phi_patient_process_noise_sd       = 0.3,   # patient-level variation on logit(φ) scale

    # Population-level process noise (shared AR(1) temporal trend across all patients)
    tr_log_sd_pop_process_noise_pop_mean  = -3,    # log(0.05), conservative
    tr_log_sd_pop_process_noise_pop_sd    = 0.5,   # allows σ ~[0.02, 0.13]
    tr_logit_phi_pop_process_noise_pop_mean = 1.4, # favors high correlation
    tr_logit_phi_pop_process_noise_pop_sd = 0.5,   # allows φ ~[0.6, 0.9]

    # Fraction module hyperparams (multi-level): c(trial, patient)
    frac_logit_loc_pop_mean = frac_logit_loc_pop_mean,
    frac_logit_loc_pop_sd = frac_logit_loc_pop_sd,
    frac_sd_level_intercept_sd = c(0.25, 0.25),
    frac_coef_qr_pop_mean = as.array(coef_elicited_priors$coef_mean),
    frac_coef_qr_pop_sd = as.array(coef_elicited_priors$coef_sd),
    frac_sd_level_slope_sd = list(
      trial = rep(0.05, n_covar),
      patient = rep(0.03, n_covar)
    ),

    # Initial state proportion module hyperparams (multi-level): c(trial, patient)
    init_logit_loc_pop_mean = init_logit_loc_pop_mean,
    init_logit_loc_pop_sd = init_logit_loc_pop_sd,
    init_sd_level_intercept_sd = c(0.6, 0.5),
    init_coef_qr_pop_mean = as.array(coef_elicited_priors$coef_mean),
    init_coef_qr_pop_sd = as.array(coef_elicited_priors$coef_sd),
    init_sd_level_slope_sd = list(
      trial = rep(0.10, n_covar),
      patient = rep(0.08, n_covar)
    ),

    # Growth lag
    growth_lag_mean = 2.7,
    growth_lag_sd = 0.5,
    patient_log_growth_lag_sd_sd = 1,
    log_growth_transition_rate_sd = 1,

    rate_corr_param = 2.0,

    # Other events baseline hazard GP hyperparameters
    oe_log_lambda_gp_pop_intercept_mean = array(
      -4.5,
      dim = c(stan_data$n_causes)
    ),
    oe_log_lambda_gp_pop_intercept_sd = array(0.5, dim = c(stan_data$n_causes)),
    # inv_gamma(5, 1.3) matches half-normal(0, 0.4) but avoids values near 0
    oe_log_lambda_gp_pop_alpha_alpha = array(5.0, dim = c(stan_data$n_causes)),
    oe_log_lambda_gp_pop_alpha_beta = array(1.3, dim = c(stan_data$n_causes)),
    oe_log_lambda_gp_pop_rho_alpha = array(8.0, dim = c(stan_data$n_causes)),
    oe_log_lambda_gp_pop_rho_beta = array(12.0, dim = c(stan_data$n_causes)),
    # inv_gamma(5, 0.8) matches half-normal(0, 0.25) but avoids values near 0
    oe_log_lambda_gp_trial_alpha_alpha = array(
      5.0,
      dim = c(stan_data$n_causes)
    ),
    oe_log_lambda_gp_trial_alpha_beta = array(0.8, dim = c(stan_data$n_causes)),
    oe_log_lambda_gp_trial_rho_alpha = array(5.0, dim = c(stan_data$n_causes)),
    oe_log_lambda_gp_trial_rho_beta = array(7.0, dim = c(stan_data$n_causes)),
    oe_log_lambda_gp_trial_intercept_sd_sd = array(
      0.3,
      dim = c(stan_data$n_causes)
    ),

    # Other events covariate effect hyperparameters
    oe_tumor_coef_pop_mean = array(
      rep(0, stan_data$n_tumor_covar),
      dim = c(stan_data$n_causes, stan_data$n_tumor_covar)
    ),
    oe_tumor_coef_pop_sd = array(
      rep(0.5, stan_data$n_tumor_covar),
      dim = c(stan_data$n_causes, stan_data$n_tumor_covar)
    ),
    oe_covar_coef_qr_pop_mean = array(
      -coef_elicited_priors$coef_mean,
      dim = c(stan_data$n_causes, stan_data$n_covar)
    ),
    oe_covar_coef_qr_pop_sd = array(
      coef_elicited_priors$coef_sd,
      dim = c(stan_data$n_causes, stan_data$n_covar)
    ),
    # Multi-level random slope SD hyperpriors for non-tumor covariates
    # Dimensions: [n_causes, n_levels, n_covar]
    oe_sd_level_slope_sd = array(
      rep(0.15, stan_data$n_causes * n_levels * stan_data$n_covar),
      dim = c(stan_data$n_causes, n_levels, stan_data$n_covar)
    ),

    log_lod_sd = 0.2
  )
}

get_pfs_priors <- function() {
  lst(
    log_lambda_gp_intercept_mean = -4.5,
    log_lambda_gp_intercept_sd = 0.25,
    log_lambda_gp_alpha_sd = 0.5,
    log_lambda_gp_rho_alpha = 7.3,
    log_lambda_gp_rho_beta = 7.5,
    log_lambda_gp_trial_alpha_sd = 0.25,
    log_lambda_gp_trial_intercept_sd_sd = 0.25,

    tumor_stim_pop_intercept_sd = 0.35,
    tumor_stim_pop_coef_sd = c(0.15, 0.15, 0.15, 0.05, 0.05),
    tumor_stim_trial_coef_sd_sd = c(0.3, 0.125, 0.125, 0.125, 0.125, 0.125),
    tumor_stim_location_coef_sd_sd = tumor_stim_trial_coef_sd_sd,

    orr_pop_coef_sd = 0.125
  )
}

get_pfs_conf_resp_priors <- function(stan_data) {
  get_pfs_priors() |>
    list_assign(
      log_lambda_gp_intercept_mean = -2.5,
      log_lambda_gp_intercept_sd = 0.5,
      log_lambda_gp_alpha_sd = 1,
      log_lambda_gp_rho_alpha = 3,
      log_lambda_gp_rho_beta = 10,

      covar_effect_mean = rep(0, stan_data$n_covar),
      covar_effect_sd = rep(0.15, stan_data$n_covar),
      tumor_stim_pop_coef_sd = c(0.2, 0.2, 0.2, 0.15, 0.15),
      conf_resp_effect_mean = 0,
      conf_resp_effect_sd = 0.3,

      covar_trial_sd_sd = 0.1,
      covar_trial_corr_eta = 2
    )
}

get_confirmed_resp_priors <- function() {
  lst(
    crcr_covar_effect_mean = 0,
    crcr_covar_effect_sd = 0.2,
    crcr_tumor_stim_pop_coef_sd = c(0.5, 0.5, 0.2, 0.15, 0.15),

    log_crcr_lambda_gp_intercept_mean = rep(-2.5, 2),
    log_crcr_lambda_gp_intercept_sd = 0.75,

    log_crcr_lambda_gp_alpha_sd = 0.6,
    log_crcr_lambda_gp_rho_alpha = 3,
    log_crcr_lambda_gp_rho_beta = 10,
    log_crcr_lambda_gp_trial_alpha_sd = 0.25,
    log_crcr_lambda_gp_trial_intercept_sd_sd = 0.25,

    crcr_covar_trial_sd_sd = 0.1,
    crcr_covar_trial_corr_eta = 2
  )
}

# nolint end: object_usage_linter
