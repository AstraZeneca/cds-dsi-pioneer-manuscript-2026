get_tumor_priors <- function() {
  lst(
    tumor_mean_mean = 3,
    tumor_mean_sd = 5,
    tumor_measure_error_sd_sd = 0.5,
    
    pop_tumor_gp_alpha_sd = 2.5,
    pop_tumor_gp_rho_alpha = 4,
    pop_tumor_gp_rho_beta = 8, 
    
    trial_tumor_gp_alpha_sd = 2.5,
    trial_tumor_gp_rho_alpha = pop_tumor_gp_rho_alpha,
    trial_tumor_gp_rho_beta = pop_tumor_gp_rho_beta,
    trial_tumor_gp_intercept_sd_sd = 2,
    
    patient_tumor_gp_alpha_sd = 3,
    patient_tumor_gp_rho_alpha = 1.5,
    patient_tumor_gp_rho_beta = trial_tumor_gp_rho_beta,
    patient_tumor_gp_rho_meanlog = 3, 
    patient_tumor_gp_rho_sdlog = 0.6,
    patient_tumor_gp_intercept_sd_sd = 3,
    
    log_trial_rho_sd_sd = 1,
    log_patient_rho_sd_sd = 1.75,
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