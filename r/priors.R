get_tumor_priors <- function(stan_data, coef_elicited_priors) {
  # Directly specified priors (simplified)
  # Choose log-total rate prior similar to historical center; adjust if needed.
  pop_log_total_rate_mean <- -3.2   # simplified fixed value
  pop_log_total_rate_sd   <-  2.0   # broader to absorb previous mapping uncertainty
  # Fraction (logit) prior keeps strong shrinkage bias (>0.5 fraction)
  pop_decrease_frac_logit_mean <- 2.2
  pop_decrease_frac_logit_sd   <- 1.2

  # Covariate effects: interpret 'worsens' as increasing fraction-driven shrinkage.
  dir_sign <- with(coef_elicited_priors, case_when(
    fct_match(effect_direction, "worsens")  ~  1,
    fct_match(effect_direction, "improves") ~ -1,
    TRUE                                     ~  0
  ))
  pop_decrease_frac_logit_coef_mean <- coef_elicited_priors$mean * dir_sign
  pop_decrease_frac_logit_coef_sd   <- coef_elicited_priors$sd

  lst(
    # GP hyperparameters
    pop_tumor_gp_rho_meanlog = 3,
    pop_tumor_gp_rho_sdlog = 0.6,
    log_patient_tumor_gp_rho_sd_sd = 1.75,
    pop_decrease_process_sd_sd = 0.1,
    pop_growth_process_sd_sd = 0.1,
    process_corr_param = 2.0,
    measure_sd_sd = 0.05,

    # Process parameters
    decrease_process_alpha = 9.7,
    decrease_process_beta = 38.4,
    growth_process_alpha = 9.7,
    growth_process_beta = 38.4,

    # Rate parameters (derived mapping)
    pop_log_total_rate_mean = pop_log_total_rate_mean,
    pop_log_total_rate_sd   = pop_log_total_rate_sd,
    pop_decrease_frac_logit_mean = pop_decrease_frac_logit_mean,
    pop_decrease_frac_logit_sd   = pop_decrease_frac_logit_sd,
    patient_log_total_rate_sd_sd      = 0.5,
    patient_decrease_frac_logit_sd_sd = 0.3,
    trial_log_total_rate_sd_sd        = 0.6,

    # Growth lag
    growth_lag_mean = 2.7,
    growth_lag_sd   = 0.5,
    patient_log_growth_lag_sd_sd = 1,
    log_growth_transition_rate_sd = 1,

    rate_corr_param = 2.0,

    # Initial state proportions
    pop_decrease_prop_logis_mean = -1.0,
    pop_decrease_prop_logis_sd   = 1.5,
    patient_decrease_prop_logis_sd_sd = 1.0,
    trial_decrease_prop_logis_sd_sd   = 1.5,

    # Baseline hazard GP
    log_lambda_gp_pop_intercept_mean = array(-4.5),
    log_lambda_gp_pop_intercept_sd   = array(0.5),
    log_lambda_gp_pop_alpha_sd       = array(0.4),
    log_lambda_gp_pop_rho_alpha      = array(8.0),
    log_lambda_gp_pop_rho_beta       = array(12.0),
    log_lambda_gp_trial_alpha_sd     = array(0.25),
    log_lambda_gp_trial_rho_alpha    = array(5.0),
    log_lambda_gp_trial_rho_beta     = array(7.0),
    log_lambda_gp_trial_intercept_sd_sd = array(0.3),

    # Covariate effects (fraction & initial state proportion)
    pop_decrease_frac_logit_coef_mean = pop_decrease_frac_logit_coef_mean,
    pop_decrease_frac_logit_coef_sd   = pop_decrease_frac_logit_coef_sd,
    pop_decrease_prop_logis_coef_mean = with(coef_elicited_priors,
      mean * if_else(fct_match(effect_direction, "worsens"), -1, 1)),
    pop_decrease_prop_logis_coef_sd = coef_elicited_priors$sd,
    trial_decrease_frac_logit_coef_sd_sd = rep(0.05, stan_data$n_covar),

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
