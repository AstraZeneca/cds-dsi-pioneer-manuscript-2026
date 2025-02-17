functions {
  #include "functions.stan"
}

data {
  #include "data.stan"
}

transformed data {
  #include "transformed_data.stan"
 
  // The calendar_week used to be an offset within trials and now it is global. Be careful with old bootstrap code that might have relied on the 
  // within trial aspect of this.
  array[n_patients] int<lower = 1> sorted_calendar_week;
  
  for (s in 1:n_trials) {
    int patient_pos = trial_patient_pos[s];
    int patient_end = trial_patient_pos[s + 1] - 1;
    
    sorted_calendar_week[patient_pos:patient_end] = sort_asc(calendar_week[patient_pos:patient_end]);
  }
  
  #include "../bootstrap/leave_out_trial_bootstrap_transformed_data.stan"
  // #include "../fixed_bootstrap_transformed_data.stan"
}

parameters {
  #include "parameters.stan"
}

transformed parameters {
  #include "transformed_parameters.stan"
  
  matrix[n_patients, no_prop_hazard || pfs_only ? 1 : n_causes] patient_response_lp; 
  
  patient_response_lp[, 1] = calc_pch_loglik(pfs, right_censored, interval_censored, pfs_ignore_interval_censoring, log_cond_prob_surv[1]);
  
  if (!(no_prop_hazard || pfs_only)) {
    patient_response_lp[, 2] = calc_pch_loglik(pfs, right_censored, interval_censored, pfs_ignore_interval_censoring, log_cond_prob_surv[2]);
  }
}

model {
  #include "priors.stan"
  
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
          log_crcr_cond_prob_surv[training_crcr_intervals]
        );
      } else {
        target += reduce_sum(
          partial_sum_crcr_lupmf, last_unclassified_response_week, crcr_grain_size,
          confirmed_response_cause, 
          early_confirmed_response_censored, crcr_ignore_interval_censoring ? zeros_int_array(n_patients) : confirmed_response_interval_censored, 
          log_crcr_cond_prob_surv
        );
      }
    }
    
    // PFS model
      
    for (s in 1:n_trials) {
      int patient_pos = trial_patient_pos[s];
      int patient_end = trial_patient_pos[s + 1] - 1; 
      
      if (s != leave_out_trial) {
        if (no_prop_hazard || pfs_only) {
          target += sum(patient_response_lp[patient_pos:patient_end, 1]);
        } else {
          for (i in patient_pos:patient_end) {
            if (confirmed_response_censored[i]) { // Unclassified
              target += log_sum_exp(log_cif[1, i, max_confresp_week] + patient_response_lp[i, 1], log_cif[2, i, max_confresp_week] + patient_response_lp[i, 2]) -
                log_sum_exp(log_cif[1, i, max_confresp_week], log_cif[2, i, max_confresp_week]);
            } else {
              target += patient_response_lp[i, confirmed_response_cause[i]];
            }
          }
        }
      }
    }
  }
}

generated quantities {
  #include "../crcr/crcr_gen_quants.stan"
  // #include "../bootstrap/leave_out_trial_bootstrap_gen_quants.stan"
 
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
  
  matrix<lower = 0, upper = 1>[fit_data ? n_patients : 0, max_all_t] admin_brier_score;

  { 
    profile("gen_pfs") {
      for (i in 1:n_patients) {
        int pfs_interval_pos = patient_pfs_interval_pos[i];
        int pfs_interval_end = patient_pfs_interval_pos[i + 1] - 1;
        
        if (right_censored[i] || interval_censored[i]) {
          (forecast_pfs[i], forecast_censored[i]) = 
            survival_time_rng(log_cond_prob_surv[no_prop_hazard ? 1 : sim_confirmed_response[i] + 1, i], pfs[i], right_censored[i], interval_censored[i]);
        } else {
          forecast_censored[i] = 0;
          forecast_pfs[i] = pfs[i];
        }
        
        (sim_pfs[i], sim_censored[i]) = survival_time_rng(log_cond_prob_surv[no_prop_hazard ? 1 : sim_confirmed_response[i] + 1, i]); 
      }
     
      if (fit_data) { 
        admin_brier_score = calc_admin_brier_score(
          pfs, admin_right_censored_week, interval_censored, 
          confirmed_response_cause, confirmed_response_censored,
          pfs_ignore_interval_censoring, log_cif[, , max_confresp_week], log_cond_prob_surv); 
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
  
  vector<lower = 0, upper = 1>[no_prop_hazard ? 0 : n_trials] trial_c_index;
 
  if (!no_prop_hazard) { 
    for (s in 1:n_trials) {
      int patient_pos = trial_patient_pos[s];
      int patient_end = trial_patient_pos[s + 1] - 1;
      
      trial_c_index[s] = calc_c_index(
        pfs[patient_pos:patient_end], right_censored[patient_pos:patient_end], 
        confirmed_response[patient_pos:patient_end], confirmed_response_censored[patient_pos:patient_end],
        time_invariant_log_hazard_ratio,
        // prob_non_response[patient_pos:patient_end] 
        log_cif[,, max_confresp_week] 
      );
    }
  }
  
  vector[gen_log_lik || prior_sense ? n_training_patients : 0] log_lik = rep_vector(0, gen_log_lik || prior_sense ? n_training_patients : 0);
  real lprior = 0;
  
  #include "../crcr/crcr_log_lik_prior_sense.stan"
  
  vector[gen_log_lik || prior_sense ? n_training_patients : 0] crcr_log_lik = log_lik; 

  if (gen_log_lik || prior_sense) {
    int log_lik_pos = 1;
    
    for (s in 1:n_trials) {
      int patient_pos = trial_patient_pos[s];
      int patient_end = trial_patient_pos[s + 1] - 1; 
      
      if (s != leave_out_trial) {
        for (i in patient_pos:patient_end) {
          if (no_prop_hazard || pfs_only) {
            log_lik[log_lik_pos] += patient_response_lp[i, 1];
          } else if (confirmed_response_censored[i]) { // Unclassified
            int conf_resp_interval_pos = patient_conf_resp_interval_pos[i];
            int conf_resp_interval_end = patient_conf_resp_interval_pos[i + 1] - 1;
          
            log_lik[log_lik_pos] += log_sum_exp(log_cif[1, i, max_confresp_week] + patient_response_lp[i, 1], log_cif[2, i, max_confresp_week] + patient_response_lp[i, 2]) -
              log_sum_exp(log_cif[1, i, max_confresp_week], log_cif[2, i, max_confresp_week]);
            
          } else {
            log_lik[log_lik_pos] += patient_response_lp[i, confirmed_response_cause[i]];
          }
          
          log_lik_pos += 1;
        }
      }
    }
  }
  
  if (prior_sense && !no_prop_hazard) {
    for (s in 1:n_prop_separate_trials) {
      lprior += normal_lpdf(tumor_stim_pop_coef[s] | 0, tumor_stim_pop_coef_sd[s]) + normal_lpdf(covar_effect[s] | covar_effect_mean[s], covar_effect_sd[s]) +
        normal_lpdf(conf_resp_effect[s] | conf_resp_effect_sd[s], conf_resp_effect_sd[s]); 
    }
    
    if (add_trial_level_prop_hazard) {
      lprior += normal_lpdf(covar_trial_sd | 0, covar_trial_sd_sd); // + lkj_corr_cholesky_lpdf(L_covar_trial_corr | covar_trial_corr_eta);
    }
  }
  
  #include "../baseline_hazard/baseline_hazard_log_lik_prior_sense.stan"
  
  vector[(gen_log_lik || prior_sense) && log_lik_trial > 0 && leave_out_trial == 0 ? n_trial_patients[log_lik_trial] : 0] trial_log_lik;
  vector[(gen_log_lik || prior_sense) && log_lik_trial > 0 && leave_out_trial == 0 ? n_trial_patients[log_lik_trial] : 0] trial_crcr_log_lik;
  matrix<lower = 0, upper = 1>[(gen_log_lik || prior_sense) && log_lik_trial > 0 && leave_out_trial == 0 ? n_trial_patients[log_lik_trial] : 0, max_all_t] trial_admin_brier_score;
  
  if (rows(trial_log_lik) > 0) {
    int pos = trial_patient_pos[log_lik_trial];
    int end = trial_patient_pos[log_lik_trial + 1] - 1;
    
    trial_log_lik = log_lik[pos:end];
    trial_crcr_log_lik = crcr_log_lik[pos:end];
    trial_admin_brier_score = admin_brier_score[pos:end];
  }
}
