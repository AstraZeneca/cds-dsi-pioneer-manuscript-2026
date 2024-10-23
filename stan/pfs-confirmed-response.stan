functions {
  #include "extern_util.stan"
  #include "util.stan"
  #include "pfs_functions.stan"
  #include "extern_pfs_functions.stan"
  #include "crcr/crcr_functions.stan"
  #include "bootstrap/leave_out_trial_bootstrap_functions.stan"
}

data {
  // Model settings
  int<lower = 0, upper = 1> fit_data; // If 0, just do prior prediction
  int<lower = 0, upper = 1> gen_interval_censored; // Should the generated PFS be interval censored?
  int<lower = 0, upper = 1> crcr_ignore_interval_censoring; // Treat observed confirmed response week as true and ignore t_measure.
  int<lower = 0, upper = 1> pfs_ignore_interval_censoring; // Treat observed PFS as true pfs and ignore t_measure.
  int<lower = 0, upper = 1> gen_log_lik;
  int<lower = 0, upper = 1> prior_sense;
  
  // Hierarchical settings 
  int<lower = 0, upper = 1> add_trial_level;
  
  // This is the data that is shared with the tumor model 
  #include "base_data.stan"
  
  #include "bootstrap/leave_out_trial_bootstrap_data.stan"

  array[n_patients] int<lower = 0> pfs; // How many periods after baseline did patient survive.
  array[n_patients] int<lower = 0> death_week; 
  array[n_patients] int<lower = 0, upper = 1> right_censored;
  array[n_patients] int<lower = 0> interval_censored;
 
  #include "crcr/crcr_data.stan" 

  // Hyperparam
  #include "baseline_hazard/baseline_hazard_hyperparam.stan"
  #include "crcr/crcr_hyperparam.stan"
  
  vector<lower = 0>[2] tumor_stim_pop_coef_sd;
  real<lower = 0> conf_resp_effect_sd;
  vector<lower = 0>[n_covar] covar_effect_sd;
 
  real<lower = 0> covar_trial_sd_sd;
  real<lower = 0> covar_trial_corr_eta; 
}

transformed data {
  int gen_pfs = 1;
  
  #include "base_transformed_data.stan" 
  #include "pfs_transformed_data.stan"
  #include "crcr/crcr_transformed_data.stan"
  #include "bootstrap/leave_out_trial_bootstrap_transformed_data.stan"
  #include "fixed_bootstrap_transformed_data.stan"
  
  int n_missing_confirmed_response = sum(confirmed_response_censored); 
  int n_obs_confirmed_response = n_patients - n_missing_confirmed_response; 
  array[n_obs_confirmed_response] int<lower = 1, upper = n_patients> obs_confirmed_response;
  array[n_missing_confirmed_response] int<lower = 1, upper = n_patients> missing_confirmed_response;
  
  (obs_confirmed_response, missing_confirmed_response) = get_mask_idx(confirmed_response_censored);
  
  int crcr_grain_size = 83;
  
  array[n_patients + 1] int<lower = 1> patient_pfs_interval_pos = linspaced_int_array(n_patients + 1, 1, n_patients * max_all_t + 1);
  
}

parameters {
  #include "baseline_hazard/baseline_hazard_parameters.stan"
  #include "crcr/crcr_parameters.stan"
  
  vector[n_tumor_covar] tumor_stim_pop_coef;
  vector[n_covar] covar_effect;  
  real conf_resp_effect;
  
  vector<lower = 0>[add_trial_level ? n_tumor_covar + n_covar + 1 : 0] covar_trial_sd;
  cholesky_factor_corr[add_trial_level ? n_tumor_covar + n_covar + 1 : 0] L_covar_trial_corr;
  matrix[n_tumor_covar + n_covar + 1, add_trial_level ? n_trials : 0] raw_covar_trial_coef;
}


transformed parameters {
  #include "baseline_hazard/baseline_hazard_transformed_parameters.stan"
  #include "crcr/crcr_transformed_parameters.stan"
  
  matrix<upper = 0>[n_time_periods, n_causes] log_cond_prob_surv;
  matrix[n_patients, n_causes] time_invariant_log_hazard_ratio = rep_matrix(tumor_sum_covar * tumor_stim_pop_coef + covar_design_matrix * covar_effect, n_causes);
  
  // The last one is the confirmed response effect 
  time_invariant_log_hazard_ratio[, n_causes] += conf_resp_effect; 
  
  matrix[n_tumor_covar + n_covar + 1, add_trial_level ? n_trials : 0] covar_trial_coef;
 
  if (add_trial_level) {
    covar_trial_coef = diag_pre_multiply(covar_trial_sd, L_covar_trial_corr) * raw_covar_trial_coef;
  }
  
  { // Calculate patient-interval conditional probability of disease progression.
    int patient_pos = 1;
    
    profile("log surv loop") {
      for (s in 1:n_trials) {
        int patient_end = patient_pos + n_trial_patients[s] - 1; 
        
        if (add_trial_level) {
          time_invariant_log_hazard_ratio[patient_pos:patient_end] += rep_matrix(
            tumor_sum_covar[patient_pos:patient_end] * covar_trial_coef[:n_tumor_covar, s] +
            covar_design_matrix[patient_pos:patient_end] * covar_trial_coef[(n_tumor_covar + 1):(n_tumor_covar + n_covar), s],
            n_causes
          );

          // The last one is the confirmed response effect
          time_invariant_log_hazard_ratio[patient_pos:patient_end, n_causes] += covar_trial_coef[n_tumor_covar + n_covar + 1, s];
        }
    
        for (i in patient_pos:patient_end) {
          int pfs_interval_pos = patient_pfs_interval_pos[i];
          int pfs_interval_end = patient_pfs_interval_pos[i + 1] - 1;
          
          log_cond_prob_surv[pfs_interval_pos:pfs_interval_end] = 
            - exp(rep_matrix(log_trial_lambda[s, 1:max_all_t], n_causes) + rep_matrix(time_invariant_log_hazard_ratio[i], max_all_t));
        }
        
        patient_pos = patient_end + 1; 
      }
    }
  }
  
  matrix<upper = 0>[n_patients, n_causes] patient_response_lp = append_col( 
    calc_pch_loglik2(pfs, right_censored, interval_censored, pfs_ignore_interval_censoring, log_cond_prob_surv[, 1], max_all_t, rep_array(1, n_patients)),
    calc_pch_loglik2(pfs, right_censored, interval_censored, pfs_ignore_interval_censoring, log_cond_prob_surv[, 2], max_all_t, rep_array(1, n_patients))
  );
}

model {
  // Priors
  
  profile("crcr priors") {
    #include "crcr/crcr_priors.stan"
  }
  
  #include "baseline_hazard/baseline_hazard_priors.stan"
  
  profile("pfs priors") {
    tumor_stim_pop_coef ~ normal(0, tumor_stim_pop_coef_sd);
    conf_resp_effect ~ normal(0, conf_resp_effect_sd); 
    covar_effect ~ normal(0, covar_effect_sd);
    
    if (add_trial_level) {
      covar_trial_sd ~ normal(0, covar_trial_sd_sd);
      L_covar_trial_corr ~ lkj_corr_cholesky(covar_trial_corr_eta);
      to_vector(raw_covar_trial_coef) ~ std_normal();
    }
  }
  
  // Likelihood
  
  if (fit_data) {
    // Confirmed response model
    
    profile("crcr loglik") {
      if (leave_out_trial > 0) {
        target += reduce_sum(
          partial_sum_crcr_lupmf, last_unclassified_response_week[training_patients], crcr_grain_size,
          confirmed_response_cause[training_patients], 
          early_confirmed_response_censored[training_patients], 
          crcr_ignore_interval_censoring ? zeros_int_array(n_training_patients) : confirmed_response_interval_censored[training_patients], 
          log_crcr_cond_prob_surv[training_crcr_intervals], max_confresp_week
        );
      } else {
        target += reduce_sum(
          partial_sum_crcr_lupmf, last_unclassified_response_week, crcr_grain_size,
          confirmed_response_cause, 
          early_confirmed_response_censored, crcr_ignore_interval_censoring ? zeros_int_array(n_patients) : confirmed_response_interval_censored, 
          log_crcr_cond_prob_surv, max_confresp_week
        );
      }
    }
    
    // PFS model
      
    for (s in 1:n_trials) {
      int patient_pos = trial_patient_pos[s];
      int patient_end = trial_patient_pos[s + 1] - 1; 
      
      if (s != leave_out_trial) {
        for (i in patient_pos:patient_end) {
          if (confirmed_response_censored[i]) { // Unclassified
            int conf_resp_interval_pos = patient_conf_resp_interval_pos[i];
            int conf_resp_interval_end = patient_conf_resp_interval_pos[i + 1] - 1;
          
            real prob_non_response = calc_cif(1, log_crcr_cond_prob_surv[conf_resp_interval_pos:conf_resp_interval_end], max_confresp_week).2[1, 1];
          
            target += log_mix(prob_non_response, patient_response_lp[i, 1], patient_response_lp[i, 2]);
          } else {
            target += patient_response_lp[i, confirmed_response_cause[i]];
          }
        }
      }
    }
  }
}

generated quantities {
  #include "crcr/crcr_gen_quants.stan"
  #include "bootstrap/leave_out_trial_bootstrap_gen_quants.stan"
  
  array[n_patients] int<lower = 0, upper = 1> sim_confirmed_response = confirmed_response;
  
  sim_confirmed_response[missing_confirmed_response] = bernoulli_rng(prob_cause[missing_confirmed_response, 2]); 
  
  array[n_patients] int<lower = 0> sim_pfs; 
  array[n_patients] int<lower = 0, upper = 1> sim_censored; 
  
  real<lower = 0> sim_median_pfs;
  vector<lower = 0>[n_trials] sim_trial_median_pfs;
  
  vector<lower = 0, upper = 1>[max_all_t + 1] km_est; 
  array[n_trials] vector<lower = 0, upper = 1>[max_all_t + 1] trial_km_est; 

  { 
    array[n_causes] matrix[max_all_t, n_patients] mat_log_cond_prob_surv;
    
    profile("gen_pfs") {
      for (i in 1:n_patients) {
        int pfs_interval_pos = patient_pfs_interval_pos[i];
        int pfs_interval_end = patient_pfs_interval_pos[i + 1] - 1;
        
        for (k in 1:n_causes) {
          mat_log_cond_prob_surv[k, , i] = log_cond_prob_surv[pfs_interval_pos:pfs_interval_end, k];
        }
        
        (sim_pfs[i], sim_censored[i]) = survival_time_rng(mat_log_cond_prob_surv[sim_confirmed_response[i] + 1, , i]); 
      }
      
      sim_median_pfs = survival_median(sim_pfs, max_all_t).1; 
      km_est = estimate_kaplan_meier(sim_pfs, sim_censored, max_all_t).1; 
      
      for (s in 1:n_trials) {
        int patient_pos = trial_patient_pos[s];
        int patient_end = trial_patient_pos[s + 1] - 1;
        
        sim_trial_median_pfs[s] = survival_median(sim_pfs[patient_pos:patient_end], max_all_t).1; 
        trial_km_est[s] = estimate_kaplan_meier(sim_pfs[patient_pos:patient_end], sim_censored[patient_pos:patient_end], max_all_t).1; 
      }
    }
  }
  
  vector[gen_log_lik || prior_sense ? n_training_patients : 0] log_lik = rep_vector(0, gen_log_lik || prior_sense ? n_training_patients : 0);
  real lprior = 0;

  #include "crcr/crcr_log_lik_prior_sense.stan"
  
  if (gen_log_lik || prior_sense) {
    int log_lik_pos = 1;
    
    for (s in 1:n_trials) {
      int patient_pos = trial_patient_pos[s];
      int patient_end = trial_patient_pos[s + 1] - 1; 
      
      if (s != leave_out_trial) {
        for (i in patient_pos:patient_end) {
          if (confirmed_response_censored[i]) { // Unclassified
            int conf_resp_interval_pos = patient_conf_resp_interval_pos[i];
            int conf_resp_interval_end = patient_conf_resp_interval_pos[i + 1] - 1;
          
            real prob_non_response = calc_cif(1, log_crcr_cond_prob_surv[conf_resp_interval_pos:conf_resp_interval_end], max_confresp_week).2[1, 1];
          
            log_lik[log_lik_pos] += log_mix(prob_non_response, patient_response_lp[i, 1], patient_response_lp[i, 2]);
          } else {
            log_lik[log_lik_pos] += patient_response_lp[i, confirmed_response_cause[i]];
          }
          
          log_lik_pos += 1;
        }
      }
    }
  }
  
  if (prior_sense) {
    lprior += normal_lpdf(log_lambda_gp_alpha | 0, log_lambda_gp_alpha_sd) + inv_gamma_lpdf(log_lambda_gp_rho | log_lambda_gp_rho_alpha, log_lambda_gp_rho_beta) +
      normal_lpdf(log_lambda_gp_intercept | log_lambda_gp_intercept_mean, log_lambda_gp_intercept_sd);
    
    if (add_trial_level) {
      lprior += normal_lpdf(log_lambda_gp_trial_intercept_sd | 0, log_lambda_gp_trial_intercept_sd_sd) + normal_lpdf(log_lambda_gp_trial_alpha | 0, log_lambda_gp_trial_alpha_sd) +
      inv_gamma_lpdf(log_lambda_gp_trial_rho | log_lambda_gp_rho_alpha, log_lambda_gp_rho_beta);
    }
  }
}
