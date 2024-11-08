create_crcr_initializer <- function(stan_data, n_causes = 2) {
  sep_trial <- if (stan_data$separate_baseline_hazard) stan_data$n_trials else 1 
  
  function(chain_id) {
    init_vals <- lst(
      log_crcr_lambda_gp_intercept = matrix(
        with(stan_data, rnorm(n_causes * sep_trial, t(log_crcr_lambda_gp_intercept_mean), t(log_crcr_lambda_gp_intercept_sd))),
        nrow = sep_trial, byrow = TRUE
      )
    )
    
    if (with(stan_data, add_trial_level_baseline_hazard && !separate_baseline_hazard)) {
      init_vals <- init_vals |> 
        list_assign(
          log_crcr_lambda_gp_trial_intercept_sd = with(stan_data, abs(rnorm(n_causes, sd = log_crcr_lambda_gp_trial_intercept_sd_sd)))
        )
    }
    
    return(init_vals) 
  }
}

create_crcr_pfs_initializer <- function(stan_data, n_causes = 2) {
  crcr_init_fun <- create_crcr_initializer(stan_data, n_causes)
  sep_trial <- if (stan_data$separate_baseline_hazard) stan_data$n_trials else 1 
  
  function(chain_id) {
    init_vals <- crcr_init_fun(chain_id) |> 
      list_assign(
        log_lambda_gp_intercept = with(stan_data, rnorm(sep_trial, log_lambda_gp_intercept_mean, log_lambda_gp_intercept_sd))
      )
    
    if (with(stan_data, add_trial_level_baseline_hazard && !separate_baseline_hazard)) {
      init_vals <- init_vals |>  
        list_assign(
          log_lambda_gp_trial_intercept_sd = with(stan_data, abs(rnorm(1, sd = log_lambda_gp_trial_intercept_sd_sd)))
        )  
    }
    
    if (with(stan_data, add_trial_level_prop_hazard && !separate_prop_hazard)) {
      init_vals <- init_vals |>  
        list_assign(
          covar_trial_sd = with(stan_data, abs(rnorm(n_tumor_covar + n_covar + 1, sd = covar_trial_sd_sd)))
        )  
    }
    
    return(init_vals)
  }
}
