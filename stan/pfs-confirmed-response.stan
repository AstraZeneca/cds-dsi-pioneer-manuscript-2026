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
  int<lower = 0, upper = 1> gen_log_lik; // Calculate log likelihood for LOO-CV
  int<lower = 0, upper = 1> prior_sense; // For prior sensitivity using {priorsense}
  int<lower = 0, upper = 1> pfs_only; // Ignore the confirmed response model; don't include as predictor in proportional hazard.
  int<lower = 0, upper = 1> no_tumor_effects; // Don't include tumor size as a predictor in proportional hazard.
  
  // Hierarchical settings 
  int<lower = 0, upper = 1> add_trial_level_baseline_hazard;
  int<lower = 0, upper = 1> add_trial_level_prop_hazard;
  int<lower = 0, upper = 1 - add_trial_level_baseline_hazard> separate_baseline_hazard;
  int<lower = 0, upper = 1 - add_trial_level_prop_hazard> separate_prop_hazard;
  
  // This is the data that is shared with the tumor model 
  #include "base_data.stan"
  
  #include "bootstrap/leave_out_trial_bootstrap_data.stan"
  
  // Calculating log likelihood for a single trial. Useful if you want to compare the preformance of a model using a single trial with one that is multilevel. 
  int<lower = 0, upper = n_trials> log_lik_trial;

  array[n_patients] int<lower = 0> pfs; // How many periods after baseline did patient survive.
  array[n_patients] int<lower = 0, upper = 1> right_censored;
  array[n_patients] int<lower = 0> interval_censored;
 
  #include "crcr/crcr_data.stan" 

  // Hyperparam
  #include "baseline_hazard/baseline_hazard_hyperparam.stan"
  #include "crcr/crcr_hyperparam.stan"
  
  array[separate_prop_hazard ? n_trials : 1] vector<lower = 0>[2] tumor_stim_pop_coef_sd;
  array[separate_prop_hazard ? n_trials : 1] real conf_resp_effect_mean;
  array[separate_prop_hazard ? n_trials : 1] real<lower = 0> conf_resp_effect_sd;
  array[separate_prop_hazard ? n_trials : 1] vector[n_covar] covar_effect_mean;
  array[separate_prop_hazard ? n_trials : 1] vector<lower = 0>[n_covar] covar_effect_sd;

  // Multilevel hyperparameters for the proportional hazard parameters 
  real<lower = 0> covar_trial_sd_sd;
  real<lower = 0> covar_trial_corr_eta; // Correlation between parameters
}

transformed data {
  int gen_pfs = 1;
  
  #include "base_transformed_data.stan" 
  #include "pfs_transformed_data.stan"
  #include "crcr/crcr_transformed_data.stan"
  #include "bootstrap/leave_out_trial_bootstrap_transformed_data.stan"
  #include "fixed_bootstrap_transformed_data.stan"
 
  // Indices of the observed and missing confirmed response values 
  int n_missing_confirmed_response = sum(confirmed_response_censored); 
  int n_obs_confirmed_response = n_patients - n_missing_confirmed_response; 
  array[n_obs_confirmed_response] int<lower = 1, upper = n_patients> obs_confirmed_response;
  array[n_missing_confirmed_response] int<lower = 1, upper = n_patients> missing_confirmed_response;
  
  (obs_confirmed_response, missing_confirmed_response) = get_mask_idx(confirmed_response_censored);
  
  int crcr_grain_size = 83; // For reduce_sum()
  
  array[n_patients + 1] int<lower = 1> patient_pfs_interval_pos = linspaced_int_array(n_patients + 1, 1, n_patients * max_all_t + 1);
}

parameters {
  #include "baseline_hazard/baseline_hazard_parameters.stan"
  #include "crcr/crcr_parameters.stan"
  
  // Proportional hazard parameters
  
  array[separate_prop_hazard ? n_trials : 1] vector[n_tumor_covar] tumor_stim_pop_coef;
  array[separate_prop_hazard ? n_trials : 1] vector[n_covar] covar_effect;  
  vector[separate_prop_hazard ? n_trials : 1] conf_resp_effect;
  
  vector<lower = 0>[add_trial_level_prop_hazard ? n_tumor_covar + n_covar + 1 : 0] covar_trial_sd;
  cholesky_factor_corr[add_trial_level_prop_hazard ? n_tumor_covar + n_covar + 1 : 0] L_covar_trial_corr;
  matrix[n_tumor_covar + n_covar + 1, add_trial_level_prop_hazard ? n_trials : 0] raw_covar_trial_coef;
}


transformed parameters {
  #include "baseline_hazard/baseline_hazard_transformed_parameters.stan"
  #include "crcr/crcr_transformed_parameters.stan"
 
 // Log conditional probability of survival at each interval. Separate columns for confirmed responders and non-responders, so we can calculate the mixture log likelihood. 
  matrix<upper = 0>[n_time_periods, n_causes] log_cond_prob_surv; 
  matrix[n_patients, n_causes] time_invariant_log_hazard_ratio; // Log proportional hazard 
  
  matrix[n_tumor_covar + n_covar + 1, add_trial_level_prop_hazard ? n_trials : 0] covar_trial_coef_residual; // Multilevel variations
  array[n_trials] vector[n_tumor_covar + n_covar + 1] covar_trial_coef;
  
  if (add_trial_level_prop_hazard) {
    covar_trial_coef_residual = diag_pre_multiply(covar_trial_sd, L_covar_trial_corr) * raw_covar_trial_coef;
  }
  
  { // Calculate patient-interval conditional probability of disease progression.
    int patient_pos = 1;
    
    profile("log surv loop") {
      for (s in 1:n_trials) {
        int patient_end = patient_pos + n_trial_patients[s] - 1; 
       
        // Tumor size effect
        covar_trial_coef[s, :n_tumor_covar] = no_tumor_effects ? rep_vector(0, n_tumor_covar) : tumor_stim_pop_coef[separate_prop_hazard ? s : 1];
        // Othe patient level covariates
        covar_trial_coef[s, (n_tumor_covar + 1):(n_tumor_covar + n_covar)] = covar_effect[separate_prop_hazard ? s : 1];
        // Confirmed response status, to impute if not observed
        covar_trial_coef[s, n_tumor_covar + n_covar + 1] = pfs_only ? 0 : conf_resp_effect[separate_prop_hazard ? s : 1]; 
          
        if (add_trial_level_prop_hazard) {
          covar_trial_coef[s, (n_tumor_covar + 1):(n_tumor_covar + n_covar)] += covar_trial_coef_residual[(n_tumor_covar + 1):(n_tumor_covar + n_covar), s];
          
          if (!no_tumor_effects) {
            covar_trial_coef[s, :n_tumor_covar] += covar_trial_coef_residual[:n_tumor_covar, s];
          }
          
          if (!pfs_only) {
            covar_trial_coef[s, n_tumor_covar + n_covar + 1] += covar_trial_coef_residual[n_tumor_covar + n_covar + 1, s];
          }
        }
        
        time_invariant_log_hazard_ratio[patient_pos:patient_end] = rep_matrix(
          tumor_sum_covar[patient_pos:patient_end] * covar_trial_coef[s, :n_tumor_covar] + 
          covar_design_matrix[patient_pos:patient_end] * covar_trial_coef[s, (n_tumor_covar + 1):(n_tumor_covar + n_covar)], 
          n_causes);
         
        if (!pfs_only) { 
          time_invariant_log_hazard_ratio[patient_pos:patient_end, n_causes] += covar_trial_coef[s, n_tumor_covar + n_covar + 1];
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
  
  matrix[n_patients, n_causes] patient_response_lp = append_col( 
    // If non-responder
    calc_pch_loglik(pfs, right_censored, interval_censored, pfs_ignore_interval_censoring, log_cond_prob_surv[, 1], max_all_t, rep_array(1, n_patients)),
    // If responder
    calc_pch_loglik(pfs, right_censored, interval_censored, pfs_ignore_interval_censoring, log_cond_prob_surv[, 2], max_all_t, rep_array(1, n_patients))
  );
  
  vector<lower = 0, upper = 1>[n_patients] prob_non_response;
  
  for (i in 1:n_patients) {
    int conf_resp_interval_pos = patient_conf_resp_interval_pos[i];
    int conf_resp_interval_end = patient_conf_resp_interval_pos[i + 1] - 1;
    
    prob_non_response[i] = calc_cif(1, log_crcr_cond_prob_surv[conf_resp_interval_pos:conf_resp_interval_end], max_confresp_week).2[1, 1];
  }
}

model {
  // Priors
  
  profile("crcr priors") {
    #include "crcr/crcr_priors.stan"
  }
  
  #include "baseline_hazard/baseline_hazard_priors.stan"
  
  profile("pfs priors") {
    for (s in 1:(separate_prop_hazard ? n_trials : 1)) {
      tumor_stim_pop_coef[s] ~ normal(0, tumor_stim_pop_coef_sd[s]);
      conf_resp_effect[s] ~ normal(conf_resp_effect_mean[s], conf_resp_effect_sd[s]); 
      covar_effect[s] ~ normal(covar_effect_mean[s], covar_effect_sd[s]);
    }
    
    if (add_trial_level_prop_hazard) {
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
          
            target += log_mix(prob_non_response[i], patient_response_lp[i, 1], patient_response_lp[i, 2]);
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
 
  // Impute confirmed response status if needed 
  array[n_patients] int<lower = 0, upper = 1> sim_confirmed_response = confirmed_response;
  sim_confirmed_response[missing_confirmed_response] = bernoulli_rng(prob_cause[missing_confirmed_response, 2]); 
 
  // Posterior predicted PFS and censoring status 
  array[n_patients] int<lower = 0> sim_pfs; 
  array[n_patients] int<lower = 0, upper = 1> sim_censored; 
  real<lower = 0> sim_median_pfs;
  vector<lower = 0>[n_trials] sim_trial_median_pfs;
  vector<lower = 0, upper = 1>[max_all_t + 1] km_est; 
  array[n_trials] vector<lower = 0, upper = 1>[max_all_t + 1] trial_km_est; 
 
  // Forecast PFS and censoring status 
  array[n_patients] int<lower = 0> forecast_pfs; 
  array[n_patients] int<lower = 0, upper = 1> forecast_censored; 
  real<lower = 0> forecast_median_pfs;
  vector<lower = 0>[n_trials] forecast_trial_median_pfs;
  vector<lower = 0, upper = 1>[max_all_t + 1] forecast_km_est; 
  array[n_trials] vector<lower = 0, upper = 1>[max_all_t + 1] forecast_trial_km_est; 
  
  real<lower = 0, upper = 1> sim_pfs6, sim_pfs9, forecast_pfs6, forecast_pfs9;
  vector<lower = 0, upper = 1>[n_trials] sim_trial_pfs6, sim_trial_pfs9, forecast_trial_pfs6, forecast_trial_pfs9;

  { 
    array[n_causes] matrix[max_all_t, n_patients] mat_log_cond_prob_surv;
    
    profile("gen_pfs") {
      for (i in 1:n_patients) {
        int pfs_interval_pos = patient_pfs_interval_pos[i];
        int pfs_interval_end = patient_pfs_interval_pos[i + 1] - 1;
        
        for (k in 1:n_causes) {
          mat_log_cond_prob_surv[k, , i] = log_cond_prob_surv[pfs_interval_pos:pfs_interval_end, k];
        }
        
        if (right_censored[i] || interval_censored[i]) {
          (forecast_pfs[i], forecast_censored[i]) = survival_time_rng(mat_log_cond_prob_surv[sim_confirmed_response[i] + 1, , i], pfs[i], right_censored[i], interval_censored[i]);
        } else {
          forecast_censored[i] = 0;
          forecast_pfs[i] = pfs[i];
        }
        
        (sim_pfs[i], sim_censored[i]) = survival_time_rng(mat_log_cond_prob_surv[sim_confirmed_response[i] + 1, , i]); 
      }
      
      sim_median_pfs = survival_median(sim_pfs, max_all_t).1; 
      km_est = estimate_kaplan_meier(sim_pfs, sim_censored, max_all_t).1; 
      sim_pfs6 = calc_pfs_n(sim_pfs, months_to_weeks(6));
      sim_pfs9 = calc_pfs_n(sim_pfs, months_to_weeks(9));
      
      forecast_median_pfs = survival_median(forecast_pfs, max_all_t).1; 
      forecast_km_est = estimate_kaplan_meier(forecast_pfs, forecast_censored, max_all_t).1; 
      forecast_pfs6 = calc_pfs_n(forecast_pfs, months_to_weeks(6));
      forecast_pfs9 = calc_pfs_n(forecast_pfs, months_to_weeks(9));
      
      for (s in 1:n_trials) {
        int patient_pos = trial_patient_pos[s];
        int patient_end = trial_patient_pos[s + 1] - 1;
        
        sim_trial_median_pfs[s] = survival_median(sim_pfs[patient_pos:patient_end], max_all_t).1; 
        trial_km_est[s] = estimate_kaplan_meier(sim_pfs[patient_pos:patient_end], sim_censored[patient_pos:patient_end], max_all_t).1; 
        sim_trial_pfs6[s] = calc_pfs_n(sim_pfs[patient_pos:patient_end], months_to_weeks(6));
        sim_trial_pfs9[s] = calc_pfs_n(sim_pfs[patient_pos:patient_end], months_to_weeks(9));
        
        forecast_trial_median_pfs[s] = survival_median(forecast_pfs[patient_pos:patient_end], max_all_t).1; 
        forecast_trial_km_est[s] = estimate_kaplan_meier(forecast_pfs[patient_pos:patient_end], forecast_censored[patient_pos:patient_end], max_all_t).1; 
        forecast_trial_pfs6[s] = calc_pfs_n(forecast_pfs[patient_pos:patient_end], months_to_weeks(6));
        forecast_trial_pfs9[s] = calc_pfs_n(forecast_pfs[patient_pos:patient_end], months_to_weeks(9));
      }
    }
  }
  
  vector<lower = 0, upper = 1>[n_trials] trial_c_index;
  
  for (s in 1:n_trials) {
    int patient_pos = trial_patient_pos[s];
    int patient_end = trial_patient_pos[s + 1] - 1;
    
    trial_c_index[s] = calc_c_index(
      pfs[patient_pos:patient_end], right_censored[patient_pos:patient_end], 
      confirmed_response[patient_pos:patient_end], confirmed_response_censored[patient_pos:patient_end],
      time_invariant_log_hazard_ratio,
      prob_non_response[patient_pos:patient_end] 
    );
  }
  
  vector[gen_log_lik || prior_sense ? n_training_patients : 0] log_lik = rep_vector(0, gen_log_lik || prior_sense ? n_training_patients : 0);
  real lprior = 0;

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
          
            log_lik[log_lik_pos] += log_mix(prob_non_response[i], patient_response_lp[i, 1], patient_response_lp[i, 2]);
          } else {
            log_lik[log_lik_pos] += patient_response_lp[i, confirmed_response_cause[i]];
          }
          
          log_lik_pos += 1;
        }
      }
    }
  }
  
  if (prior_sense) {
    for (s in 1:(separate_prop_hazard ? n_trials : 1)) {
      lprior += normal_lpdf(tumor_stim_pop_coef[s] | 0, tumor_stim_pop_coef_sd[s]) + normal_lpdf(covar_effect[s] | covar_effect_mean[s], covar_effect_sd[s]) +
        normal_lpdf(conf_resp_effect[s] | conf_resp_effect_sd[s], conf_resp_effect_sd[s]); 
    }
    
    if (add_trial_level_prop_hazard) {
      lprior += normal_lpdf(covar_trial_sd | 0, covar_trial_sd_sd) + lkj_corr_cholesky_lpdf(L_covar_trial_corr | covar_trial_corr_eta);
    }
  }
  
  #include "crcr/crcr_log_lik_prior_sense.stan"
  #include "baseline_hazard/baseline_hazard_log_lik_prior_sense.stan"
  
  vector[(gen_log_lik || prior_sense) && log_lik_trial > 0 && leave_out_trial == 0 ? n_trial_patients[log_lik_trial] : 0] trial_log_lik;
  
  if (rows(trial_log_lik) > 0) {
    int pos = trial_patient_pos[log_lik_trial];
    int end = trial_patient_pos[log_lik_trial + 1] - 1;
    
    trial_log_lik = log_lik[pos:end];
  }
}
