# nolint start: object_usage_linter

get_tumor_priors <- function(stan_data, coef_elicited_priors) {
  # Directly specified priors (simplified)
  # Choose log-total rate prior similar to historical center; adjust if needed.
  tr_loc_pop_mean <- -1.0   # simplified fixed value (formerly pop_log_total_rate_mean)
  tr_loc_pop_sd   <-  0.8   # tightened from 1.2 to reduce total variance with hierarchies
  # Fraction (logit) prior keeps strong shrinkage bias (>0.5 fraction)
  frac_logit_loc_pop_mean <- -1.0  # formerly pop_decrease_frac_logit_mean
  frac_logit_loc_pop_sd   <- 0.8   # tightened from 1.2 (critical for logit scale)
  # Initial proportion logit
  init_logit_loc_pop_mean <- 0.0   # formerly pop_decrease_prop_logis_mean
  init_logit_loc_pop_sd   <- 0.8   # tightened from 1.5 (was extremely wide!)

  lst(
    # GP hyperparameters 
    pop_tumor_gp_rho_meanlog = 3,
    pop_tumor_gp_rho_sdlog = 0.6,
    log_patient_tumor_gp_rho_sd_sd = 1.75,
    pop_decrease_process_sd_sd = 0.1,
    pop_growth_process_sd_sd = 0.1,
    process_corr_param = 2.0,
    measure_sd_sd = 0.2,

    # Process parameters 
    decrease_process_alpha = 9.7,
    decrease_process_beta = 38.4,
    growth_process_alpha = 9.7,
    growth_process_beta = 38.4,

    # Total rate module hyperparams
    tr_loc_pop_mean = tr_loc_pop_mean,
    tr_loc_pop_sd   = tr_loc_pop_sd,
    tr_sd_trial_intercept_sd   = 0.35,  # tightened from 0.6
    tr_sd_patient_intercept_sd = 0.35,  # tightened from 0.5
    tr_coef_qr_pop_mean = as.array(rep(0, stan_data$n_covar)),
    tr_coef_qr_pop_sd   = as.array(rep(1, stan_data$n_covar)),
    tr_sd_trial_slope_sd   = as.array(rep(0.15, stan_data$n_covar)),
    tr_sd_patient_slope_sd = as.array(rep(0.10, stan_data$n_covar)),

    # Fraction module hyperparams
    frac_logit_loc_pop_mean = frac_logit_loc_pop_mean,
    frac_logit_loc_pop_sd   = frac_logit_loc_pop_sd,
    frac_sd_trial_intercept_sd   = 0.25,  # tightened from 0.3 (critical for logit)
    frac_sd_patient_intercept_sd = 0.25,  # tightened from 0.3 (critical for logit)
    frac_coef_qr_pop_mean = as.array(coef_elicited_priors$coef_mean),
    frac_coef_qr_pop_sd   = as.array(coef_elicited_priors$coef_sd),
    frac_sd_trial_slope_sd   = as.array(rep(0.05, stan_data$n_covar)),
    frac_sd_patient_slope_sd = as.array(rep(0.03, stan_data$n_covar)),

    # Initial state proportion module hyperparams
    init_logit_loc_pop_mean = init_logit_loc_pop_mean,
    init_logit_loc_pop_sd   = init_logit_loc_pop_sd,
    init_sd_trial_intercept_sd   = 0.6, # tightened from 1.5 (was way too wide!)
    init_sd_patient_intercept_sd = 0.5, # tightened from 1.0
    init_coef_qr_pop_mean = as.array(coef_elicited_priors$coef_mean),
    init_coef_qr_pop_sd   = as.array(coef_elicited_priors$coef_sd),
    init_sd_trial_slope_sd   = as.array(rep(0.10, stan_data$n_covar)),
    init_sd_patient_slope_sd = as.array(rep(0.08, stan_data$n_covar)),

    # Growth lag
    growth_lag_mean = 2.7,
    growth_lag_sd   = 0.5,
    patient_log_growth_lag_sd_sd = 1,
    log_growth_transition_rate_sd = 1,

    rate_corr_param = 2.0,

    # Other events baseline hazard GP hyperparameters
    oe_log_lambda_gp_pop_intercept_mean = array(-4.5, dim = c(stan_data$n_causes)),
    oe_log_lambda_gp_pop_intercept_sd   = array(0.5, dim = c(stan_data$n_causes)),
    oe_log_lambda_gp_pop_alpha_sd       = array(0.4, dim = c(stan_data$n_causes)),
    oe_log_lambda_gp_pop_rho_alpha      = array(8.0, dim = c(stan_data$n_causes)),
    oe_log_lambda_gp_pop_rho_beta       = array(12.0, dim = c(stan_data$n_causes)),
    oe_log_lambda_gp_trial_alpha_sd     = array(0.25, dim = c(stan_data$n_causes)),
    oe_log_lambda_gp_trial_rho_alpha    = array(5.0, dim = c(stan_data$n_causes)),
    oe_log_lambda_gp_trial_rho_beta     = array(7.0, dim = c(stan_data$n_causes)),
    oe_log_lambda_gp_trial_intercept_sd_sd = array(0.3, dim = c(stan_data$n_causes)),
    
    # Other events covariate effect hyperparameters
    oe_tumor_coef_qr_pop_mean = array(rep(0, stan_data$n_tumor_covar), dim = c(stan_data$n_causes, stan_data$n_tumor_covar)),
    oe_tumor_coef_qr_pop_sd   = array(rep(1, stan_data$n_tumor_covar), dim = c(stan_data$n_causes, stan_data$n_tumor_covar)),
    oe_covar_coef_qr_pop_mean = array(rep(0, stan_data$n_covar), dim = c(stan_data$n_causes, stan_data$n_covar)),
    oe_covar_coef_qr_pop_sd   = array(rep(1, stan_data$n_covar), dim = c(stan_data$n_causes, stan_data$n_covar)),
    oe_sd_trial_tumor_slope_sd = array(rep(0.15, stan_data$n_tumor_covar), dim = c(stan_data$n_causes, stan_data$n_tumor_covar)),
    oe_sd_trial_slope_sd = array(rep(0.15, stan_data$n_covar), dim = c(stan_data$n_causes, stan_data$n_covar)),

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