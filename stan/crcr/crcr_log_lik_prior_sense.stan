if (gen_log_lik || prior_sense) {
  if (leave_out_trial > 0) {
    log_lik += calc_comp_risk_pch_loglik(
      last_unclassified_response_week[training_patients],
      confirmed_response_cause[training_patients], 
      early_confirmed_response_censored[training_patients], 
      crcr_ignore_interval_censoring ? zeros_int_array(n_training_patients) : confirmed_response_interval_censored[training_patients], 
      log_crcr_cond_prob_surv[training_patients], max_confresp_week
    );
  } else {
    log_lik += calc_comp_risk_pch_loglik(
      last_unclassified_response_week,
      confirmed_response_cause, 
      early_confirmed_response_censored, crcr_ignore_interval_censoring ? zeros_int_array(n_patients) : confirmed_response_interval_censored, 
      log_crcr_cond_prob_surv, max_confresp_week
    );
  }
}

if (prior_sense) {
  for (s in 1:n_base_separate_trials) {
    lprior += normal_lpdf(to_vector(log_crcr_lambda_gp_alpha[s]) | 0, log_crcr_lambda_gp_alpha_sd[s]) +
      inv_gamma_lpdf(to_vector(log_crcr_lambda_gp_rho[s]) | log_crcr_lambda_gp_rho_alpha[s], log_crcr_lambda_gp_rho_beta[s]) +
      normal_lpdf(to_vector(log_crcr_lambda_gp_intercept[s]) | log_crcr_lambda_gp_intercept_mean[s], log_crcr_lambda_gp_intercept_sd[s]); 
  }
  
  if (add_trial_level_baseline_hazard) {
    lprior += normal_lpdf(log_crcr_lambda_gp_trial_alpha | 0, log_crcr_lambda_gp_trial_alpha_sd) +
      inv_gamma_lpdf(log_crcr_lambda_gp_trial_rho | log_crcr_lambda_gp_rho_alpha[1], log_crcr_lambda_gp_rho_beta[1]) +
      normal_lpdf(log_crcr_lambda_gp_trial_intercept_sd | 0, log_crcr_lambda_gp_trial_intercept_sd_sd);
  }
  
  for (s in 1:n_prop_separate_trials) {
    for (k in 1:n_causes) {
      lprior += normal_lpdf(crcr_tumor_stim_pop_coef[s, , k] | 0, crcr_tumor_stim_pop_coef_sd[s]) + 
        normal_lpdf(crcr_covar_effect[s, , k] | crcr_covar_effect_mean[s], crcr_covar_effect_sd[s]);
    }
  }
  
  if (add_trial_level_prop_hazard && !no_prop_hazard) {
    lprior += normal_lpdf(crcr_covar_trial_sd | 0, crcr_covar_trial_sd_sd) + lkj_corr_cholesky_lpdf(L_crcr_covar_trial_corr | crcr_covar_trial_corr_eta);
  }
}