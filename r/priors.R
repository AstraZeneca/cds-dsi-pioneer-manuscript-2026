get_tumor_priors <- function() {
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
    growth_process_alpha = 9.7,  # Same as decrease_process_alpha
    growth_process_beta = 38.4,  # Same as decrease_process_beta
    
    # Rate parameters
    pop_log_net_rate_mean = -2.5,        # Changed from -4 to -2.5 (more reasonable net rate)
    pop_log_net_rate_sd = 0.8,           # Changed from 2 to 0.8 (much tighter)
    pop_log_rate_ratio_mean = 1.2,       # Changed from 2 to 1.2 (smaller ratio)
    pop_log_rate_ratio_sd = 0.4,         # Changed from 1 to 0.4 (much tighter)
    patient_log_net_rate_sd_sd = 0.3,    # Changed from 0.5 to 0.3 (less patient variation)
    patient_log_rate_ratio_sd_sd = 0.2,  # Changed from 0.5 to 0.2 (less patient variation)
    trial_log_net_rate_sd_sd = 0.4,      # Changed from 1.0 to 0.4 (less trial variation)
    
    # Growth lag parameters
    growth_lag_mean = 2.7,
    growth_lag_sd = 0.5,
    patient_log_growth_lag_sd_sd = 1,
    log_growth_transition_rate_sd = 1,
    
    # Correlation parameters
    rate_corr_param = 2.0,
    
    # Proportion parameters
    pop_decrease_prop_logis_mean = -1.5,    # logit^{-1}(-1.5) ≈ 0.18
    pop_decrease_prop_logis_sd = 0.5,       # Much tighter
    patient_decrease_prop_logis_sd_sd = 0.75,  # Also tighten patient-level
    trial_decrease_prop_logis_sd_sd = 1.5,    # Between patient (2.5) and population (2) 
    
    # Baseline hazard GP parameters - Population level (wrapped in array() for n_causes)
    log_lambda_gp_pop_intercept_mean = array(-4.5),    # Baseline log-hazard
    log_lambda_gp_pop_intercept_sd = array(0.5),       # Moderate uncertainty
    
    # GP variance (alpha) - controls overall variability of hazard over time
    log_lambda_gp_pop_alpha_sd = array(0.4),           # Moderate temporal variation
    
    # GP length-scale (rho) - controls smoothness of hazard over time
    log_lambda_gp_pop_rho_alpha = array(8.0),          # Shape parameter
    log_lambda_gp_pop_rho_beta = array(12.0),          # Rate parameter (mean rho ≈ 1.5 weeks)
    
    # Trial-level hierarchical effects (when add_trial_level_baseline_hazard = 1)
    log_lambda_gp_trial_alpha_sd = array(0.25),        # Smaller trial-level variation
    log_lambda_gp_trial_rho_alpha = array(5.0),        # Trial-level smoothness
    log_lambda_gp_trial_rho_beta = array(7.0),         # Structure
    log_lambda_gp_trial_intercept_sd_sd = array(0.3),   # Trial baseline variation
    
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