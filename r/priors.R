get_tumor_priors <- function(stan_data) {
  lst(
    # GP hyperparameters
    pop_tumor_gp_rho_meanlog = 3, 
    pop_tumor_gp_rho_sdlog = 0.6,
    log_patient_tumor_gp_rho_sd_sd = 1.75,
    
    # Process noise parameters
    pop_decrease_process_sd_sd = 0.1,
    pop_growth_process_sd_sd = 0.1,
    process_corr_param = 2.0,
    measure_sd_sd = 0.2,
    
    # Process parameters
    decrease_process_alpha = 9.7,
    decrease_process_beta = 38.4,
    growth_process_alpha = 9.7,  
    growth_process_beta = 38.4,
    
    # Rate parameters - ALIGNED WITH POSTERIOR EVIDENCE
    pop_log_net_rate_mean = -3.5,        # Much lower net rate to shift growth left
    pop_log_net_rate_sd = 1.8,           # Keep wide range
    pop_log_rate_ratio_mean = 2.2,       # Keep high ratio: exp(2.2) ≈ 9 (d = 9g)
    pop_log_rate_ratio_sd = 1.2,         # Very wide to allow flexibility
    patient_log_net_rate_sd_sd = 0.5,    # Allow patient variation
    patient_log_rate_ratio_sd_sd = 0.3,  # Allow patient variation  
    trial_log_net_rate_sd_sd = 0.6,      # Allow trial variation
    
    # Growth lag parameters
    growth_lag_mean = 2.7,
    growth_lag_sd = 0.5,
    patient_log_growth_lag_sd_sd = 1,
    log_growth_transition_rate_sd = 1,
    
    # Correlation parameters
    rate_corr_param = 2.0,
    
    # Proportion parameters
    pop_decrease_prop_logis_mean = -1.0,    # logit^{-1}(-1.0) ≈ 0.27
    pop_decrease_prop_logis_sd = 1.5,       # Wider but not extreme
    patient_decrease_prop_logis_sd_sd = 1.0,
    trial_decrease_prop_logis_sd_sd = 1.5,
    
    # Baseline hazard GP parameters
    log_lambda_gp_pop_intercept_mean = array(-4.5),
    log_lambda_gp_pop_intercept_sd = array(0.5),
    log_lambda_gp_pop_alpha_sd = array(0.4),
    log_lambda_gp_pop_rho_alpha = array(8.0),
    log_lambda_gp_pop_rho_beta = array(12.0),
    log_lambda_gp_trial_alpha_sd = array(0.25),
    log_lambda_gp_trial_rho_alpha = array(5.0),
    log_lambda_gp_trial_rho_beta = array(7.0),
    log_lambda_gp_trial_intercept_sd_sd = array(0.3),
    
    # Covariate effect priors for net rates (for scaled but not centered covariates)
    pop_log_net_rate_coef_mean = rep(0, stan_data$n_covar),
    pop_log_net_rate_coef_sd = rep(0.3, stan_data$n_covar),
    trial_log_net_rate_coef_sd_sd = 0.2,
    
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