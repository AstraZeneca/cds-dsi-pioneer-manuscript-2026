functions {
  #include "../extern_util.stan"
  #include "../util.stan"
  #include "../extern_pfs_functions.stan"
  #include "../pfs_functions.stan"
  #include "crcr_functions.stan"
}

data {
  // Model settings
  int<lower = 0, upper = 1> fit_data; // If 0, just do prior prediction
  int<lower = 0, upper = 1> crcr_ignore_interval_censoring; // Treat observed intervals as true intervals 
  int<lower = 0, upper = 1> gen_log_lik;
  int<lower = 0, upper = 1> prior_sense;
  
  // Hierarchical settings 
  int<lower = 0, upper = 1> add_trial_level;

  // This is the data that is shared with the tumor model 
  #include "../base_data.stan"
  #include "crcr_data.stan"
 
  // Hyperparam
  #include "crcr_hyperparam.stan"
}

transformed data {
  #include "../base_transformed_data.stan" 
  #include "crcr_transformed_data.stan"
  
  int grain_size = 83;
}

parameters {
  #include "crcr_parameters.stan"
}

transformed parameters {
  #include "crcr_transformed_parameters.stan"
}

model {
  profile("priors") {
    #include "crcr_priors.stan" 
  }
  
  if (fit_data) {
    profile("loglik") {
      // last_unclassified_response_week ~ comp_risk_pch(confirmed_response_cause, early_confirmed_response_censored, log_crcr_cond_prob_surv, max_confresp_week);
      target += reduce_sum(
        partial_sum_crcr_lupmf, last_unclassified_response_week, grain_size,
        confirmed_response_cause, 
        early_confirmed_response_censored, crcr_ignore_interval_censoring ? zeros_int_array(n_patients) : confirmed_response_interval_censored, 
        log_crcr_cond_prob_surv, max_confresp_week
      );
    }
  }
}

generated quantities {
  #include "crcr_gen_quants.stan"
  
  vector[gen_log_lik || prior_sense ? n_patients : 0] log_lik = rep_vector(0, gen_log_lik || prior_sense ? n_patients : 0);
  real lprior = 0;
  
  if (gen_log_lik || prior_sense) {
    log_lik += calc_comp_risk_pch_loglik(
      last_unclassified_response_week,  
      confirmed_response_cause, 
      early_confirmed_response_censored, crcr_ignore_interval_censoring ? zeros_int_array(n_patients) : confirmed_response_interval_censored, 
      log_crcr_cond_prob_surv, max_confresp_week
    );
  }
  
  if (prior_sense) {
    lprior += normal_lpdf(log_crcr_lambda_gp_alpha | 0, log_crcr_lambda_gp_alpha_sd) +
      inv_gamma_lpdf(log_crcr_lambda_gp_rho | log_crcr_lambda_gp_rho_alpha, log_crcr_lambda_gp_rho_beta) +
      std_normal_lpdf(to_vector(log_crcr_lambda_gp_eta)) +
      normal_lpdf(log_crcr_lambda_gp_intercept | log_crcr_lambda_gp_intercept_mean, log_crcr_lambda_gp_intercept_sd); 

    for (k in 1:n_causes) {
      lprior += normal_lpdf(crcr_tumor_stim_pop_coef[, k] | 0, crcr_tumor_stim_pop_coef_sd) + normal_lpdf(crcr_covar_effect[, k] | 0, crcr_covar_effect_sd);
    }

    // if (add_trial_level) { 
    //   lprior += normal_lpdf(log_crcr_lambda_gp_trial_alpha | 0, log_crcr_lambda_gp_trial_alpha_sd) +
    //     inv_gamma_lpdf(log_crcr_lambda_gp_trial_rho | log_crcr_lambda_gp_rho_alpha, log_crcr_lambda_gp_rho_beta) +
    //     normal_lpdf(log_crcr_lambda_gp_trial_intercept_sd | 0, log_crcr_lambda_gp_trial_intercept_sd_sd) +
    //     std_normal_lpdf(to_vector(raw_log_crcr_lambda_gp_trial_intercept)) +
    //     normal_lpdf(crcr_covar_trial_sd | 0, crcr_covar_trial_sd_sd) +
    //     lkj_corr_cholesky_lpdf(L_crcr_covar_trial_corr | crcr_covar_trial_corr_eta);
    //   
    //   for (s in 1:n_trials) {
    //     lprior += std_normal_lpdf(to_vector(raw_crcr_covar_trial_coef[s])) + std_normal_lpdf(to_vector(log_crcr_lambda_gp_trial_eta[s]));
    //   }
    // }
  }
}
