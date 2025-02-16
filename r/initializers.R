create_crcr_initializer <- function(stan_data, n_causes = 2) {
  base_sep_trial <- if (stan_data$separate_baseline_hazard) stan_data$n_trials else 1 
  prop_sep_trial <- if (stan_data$separate_prop_hazard) stan_data$n_trials else 1 
  
  function(chain_id) {
    init_vals <- lst(
      # array[n_causes] vector[n_base_separate_trials] log_crcr_lambda_gp_intercept;
      log_crcr_lambda_gp_intercept = matrix(
        with(stan_data, rnorm(n_causes * base_sep_trial, t(log_crcr_lambda_gp_intercept_mean), t(log_crcr_lambda_gp_intercept_sd))),
        ncol = base_sep_trial, byrow = FALSE 
      ),
      # array[n_causes] vector<lower = 0>[n_base_separate_trials] log_crcr_lambda_gp_alpha;
      # log_crcr_lambda_gp_alpha = matrix(0.001 + abs(rnorm(n_causes * sep_trial, sd = stan_data$log_crcr_lambda_gp_alpha_sd)), nrow = n_causes)
    )
    
    # if (!stan_data$no_prop_hazard) {
    #   init_vals <- init_vals |> 
    #     list_assign(
    #       # array[n_causes, n_prop_separate_trials] vector[no_prop_hazard ? 0 : n_tumor_covar] crcr_tumor_stim_pop_coef;
    #       crcr_tumor_stim_pop_coef = with(
    #         stan_data, 
    #         rnorm(n_causes * n_tumor_covar, 0, crcr_tumor_stim_pop_coef_sd) |>
    #           array(c(n_causes, prop_sep_trial, n_tumor_covar))  
    #       ),
    #       
    #       # array[n_causes, n_prop_separate_trials] vector[no_prop_hazard ? 0 : n_covar] crcr_covar_effect;  
    #       crcr_covar_coef = with(
    #         stan_data, 
    #         rnorm(n_causes * prop_sep_trial * n_covar, crcr_covar_effect_mean, crcr_covar_effect_sd) |>
    #           array(c(n_causes, prop_sep_trial, n_covar))  
    #       ),
    #     )
    #   
    #   if (stan_data$add_trial_level_prop_hazard) {
    #     init_vals <- init_vals |> 
    #       list_assign(
    #         # vector<lower = 0>[add_trial_level_prop_hazard && !no_prop_hazard ? n_tumor_covar + n_covar : 0] crcr_covar_trial_sd;
    #         crcr_covar_trial_sd = with(stan_data, rnorm(n_tumor_covar + n_covar, 0, crcr_covar_trial_sd_sd)),
    #         
    #         # array[n_causes, add_trial_level_prop_hazard && !no_prop_hazard ? n_trials : 0] vector[n_tumor_covar + n_covar] raw_crcr_covar_trial_coef;
    #         raw_crcr_covar_trial_coef = with(stan_data, rnorm(n_causes * prop_sep_trial * (n_tumor_covar + n_covar)) |> array(c(n_causes, prop_sep_trial, n_tumor_covar + n_covar)))
    #       )
    #   }
    #     
    # }
    
    if (with(stan_data, add_trial_level_baseline_hazard && !separate_baseline_hazard)) {
      init_vals <- init_vals |> 
        list_assign(
          # vector[add_trial_level_baseline_hazard ? n_causes : 0] log_crcr_lambda_gp_trial_intercept_sd;
          log_crcr_lambda_gp_trial_intercept_sd = with(stan_data, abs(rnorm(n_causes, sd = log_crcr_lambda_gp_trial_intercept_sd_sd)))
        )
    }
    
    return(init_vals) 
  }
}

create_crcr_pfs_initializer <- function(stan_data, n_causes = 2) {
  prop_hazard <- !stan_data$no_prop_hazard 
  crcr_init_fun <- create_crcr_initializer(stan_data, n_causes)
  baseline_sep_trial <- if (stan_data$separate_baseline_hazard) stan_data$n_trials else 1 
  prop_sep_trial <- if (stan_data$separate_prop_hazard) stan_data$n_trials else 1 
  
  function(chain_id) {
    init_vals <- crcr_init_fun(chain_id) |> 
      list_assign(
        # vector[n_base_separate_trials] log_lambda_gp_intercept;
        log_lambda_gp_intercept = with(stan_data, rnorm(baseline_sep_trial, log_lambda_gp_intercept_mean, log_lambda_gp_intercept_sd)),
        
        tumor_stim_pop_coef = if (prop_hazard) { 
          with(stan_data, matrix(rnorm(prop_sep_trial * n_tumor_covar, sd = c(t(tumor_stim_pop_coef_sd))), byrow = TRUE, nrow = prop_sep_trial))
        } else {
          array(dim = c(0, stan_data$n_tumor_covar))
        },
        covar_effect = if (prop_hazard) { 
          with(stan_data, matrix(rnorm(prop_sep_trial * n_covar, c(t(covar_effect_mean)), c(t(covar_effect_sd))), byrow = TRUE, nrow = prop_sep_trial))
        },
        conf_resp_effect = if (prop_hazard) { 
          with(stan_data, rnorm(prop_sep_trial, conf_resp_effect_mean, conf_resp_effect_sd))
        }
      ) |> 
      compact()
    
    if (with(stan_data, add_trial_level_baseline_hazard && !separate_baseline_hazard)) {
      init_vals <- init_vals |>  
        list_assign(
          log_lambda_gp_trial_intercept_sd = with(stan_data, abs(rnorm(1, sd = log_lambda_gp_trial_intercept_sd_sd))),
          log_lambda_gp_trial_alpha = with(stan_data, abs(rnorm(1, sd = log_lambda_gp_trial_alpha_sd)))
        )  
    }
    
    if (with(stan_data, add_trial_level_prop_hazard && !separate_prop_hazard) && prop_hazard) {
      init_vals <- init_vals |>  
        list_assign(covar_trial_sd = with(stan_data, abs(rnorm(n_tumor_covar + n_covar + 1, sd = covar_trial_sd_sd))))  
    }
    
    return(init_vals)
  }
}
