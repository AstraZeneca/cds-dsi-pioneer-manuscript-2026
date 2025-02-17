if (gen_log_lik || prior_sense) {
  if (leave_out_trial > 0) {
    log_lik += calc_pch_loglik(
      last_unclassified_response_week[training_patients],
      confirmed_response_cause[training_patients], 
      early_confirmed_response_censored[training_patients], 
      confirmed_response_interval_censored[training_patients], 
      crcr_ignore_interval_censoring,
      log_crcr_cond_prob_surv[training_patients]
    );
  } else {
    log_lik += calc_pch_loglik(
      last_unclassified_response_week,
      confirmed_response_cause, 
      early_confirmed_response_censored, confirmed_response_interval_censored, crcr_ignore_interval_censoring,
      log_crcr_cond_prob_surv
    );
  }
}

if (prior_sense) {
  for (k in 1:n_causes) {
    lprior += normal_lpdf(log_crcr_lambda_gp_alpha[k] | 0, log_crcr_lambda_gp_alpha_sd) +
      inv_gamma_lpdf(log_crcr_lambda_gp_rho[k] | log_crcr_lambda_gp_rho_alpha, log_crcr_lambda_gp_rho_beta); 
      
    for (s in 1:n_base_separate_trials) {
      lprior += normal_lpdf(log_crcr_lambda_gp_intercept[k] | log_crcr_lambda_gp_intercept_mean[s, k], log_crcr_lambda_gp_intercept_sd[s]);  
    }
      
    for (s in 1:n_prop_separate_trials) {
      lprior += normal_lpdf(crcr_tumor_stim_pop_coef[k, s] | 0, crcr_tumor_stim_pop_coef_sd[s]) + 
        normal_lpdf(crcr_covar_effect[k, s] | crcr_covar_effect_mean[s], crcr_covar_effect_sd[s]);
    }
  }
  
  if (add_trial_level_baseline_hazard) {
    lprior += normal_lpdf(log_crcr_lambda_gp_trial_alpha | 0, log_crcr_lambda_gp_trial_alpha_sd) +
      inv_gamma_lpdf(log_crcr_lambda_gp_trial_rho | log_crcr_lambda_gp_rho_alpha[1], log_crcr_lambda_gp_rho_beta[1]) +
      normal_lpdf(log_crcr_lambda_gp_trial_intercept_sd | 0, log_crcr_lambda_gp_trial_intercept_sd_sd);
  }
  
  if (add_trial_level_prop_hazard && !no_prop_hazard) {
    lprior += normal_lpdf(crcr_covar_trial_sd | 0, crcr_covar_trial_sd_sd); // + lkj_corr_cholesky_lpdf(L_crcr_covar_trial_corr | crcr_covar_trial_corr_eta);
  }
}