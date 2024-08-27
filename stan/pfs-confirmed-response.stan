functions {
  #include "extern_util.stan"
  #include "util.stan"
  #include "pfs_functions.stan"
  #include "extern_pfs_functions.stan"
  #include "crcr/crcr_functions.stan"
  
  array[] int forecast_pfs_rng(
    int pred_calendar_week, array[] int start_calendar_week,
    array[] int patient_ids,
    array[] int pfs_interval_pos, array[] int pfs, array[] int pfs_right_censored, matrix log_cond_prob_surv,
    array[] int conf_resp_interval_pos, array[] int conf_resp_week, array[] int conf_resp, array[] int conf_resp_censored, matrix log_crcr_cond_prob_surv
  ) {
    int n_sample_patients = size(patient_ids);
    array[n_sample_patients] int pred_pfs = pfs[patient_ids];
    
    for (i in 1:n_sample_patients) {
      if (pfs_right_censored[patient_ids[i]] || (start_calendar_week[i] + pred_pfs[i] - 1) > pred_calendar_week) {
        int conf_resp_cause = conf_resp[patient_ids[i]] + 1;
        
        int max_weeks_observed = pred_calendar_week - start_calendar_week[i] + 1;
        
        if (conf_resp_censored[patient_ids[i]] || (start_calendar_week[i] + conf_resp_week[patient_ids[i]] - 1) > pred_calendar_week) {
          int n_weeks_obs_unclassified = min(conf_resp_week[patient_ids[i]], max_weeks_observed); 
          int curr_conf_resp_interval_pos = conf_resp_interval_pos[patient_ids[i]] + n_weeks_obs_unclassified;
          int curr_conf_resp_interval_end = conf_resp_interval_pos[patient_ids[i] + 1] - 1;
         
          conf_resp_cause = competing_risks_survival_time_rng(log_crcr_cond_prob_surv[curr_conf_resp_interval_pos:curr_conf_resp_interval_end]).3; 
        }
        
        int n_weeks_obs_surv = min(pred_pfs[i], max_weeks_observed); 
        int curr_pfs_interval_pos = pfs_interval_pos[patient_ids[i]] + n_weeks_obs_surv;
        int curr_pfs_interval_end = pfs_interval_pos[patient_ids[i] + 1] - 1;
        
        pred_pfs[i] = n_weeks_obs_surv + survival_time_rng(log_cond_prob_surv[curr_pfs_interval_pos:curr_pfs_interval_end, conf_resp_cause]).1; 
      }
    }
    
    return pred_pfs;
  }
}

data {
  // Model settings
  int<lower = 0, upper = 1> fit_data; // If 0, just do prior prediction
  int<lower = 0, upper = 1> gen_interval_censored; // Should the generated PFS be interval censored?
  int<lower = 0, upper = 1> ignore_interval_censoring; // Treat observed PFS as true pfs and ignore t_measure.
  
  // Hierarchical settings 
  int<lower = 0, upper = 1> add_trial_level;
  int<lower = 0, upper = 1> add_trial_level_glm;
  
  // This is the data that is shared with the tumor model 
  #include "base_data.stan"
  
  int<lower = 0, upper = n_trials> leave_out_trial;
  int<lower = 0> n_bootstrap_sample;
  int<lower = 0> n_bootstrap_cr_maturity_rates;
  vector<lower = 0, upper = 1>[n_bootstrap_cr_maturity_rates] bootstrap_cr_maturity_rates;
  
  array[n_patients] int<lower = 0> pfs; // How many periods after baseline did patient survive.
  array[n_patients] int<lower = 0> death_week; 
  array[n_patients] int<lower = 0, upper = 1> right_censored;
 
  #include "crcr/crcr_data.stan" 

  // Hyperparam
  #include "baseline_hazard_hyperparam.stan"
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
  
  int n_missing_confirmed_response = sum(confirmed_response_censored); 
  int n_obs_confirmed_response = n_patients - n_missing_confirmed_response; 
  array[n_obs_confirmed_response] int<lower = 1, upper = n_patients> obs_confirmed_response;
  array[n_missing_confirmed_response] int<lower = 1, upper = n_patients> missing_confirmed_response;
  
  (obs_confirmed_response, missing_confirmed_response) = get_mask_idx(confirmed_response_censored);
  
  int crcr_grain_size = 83;
  
  array[n_patients + 1] int<lower = 1> patient_pfs_interval_pos = linspaced_int_array(n_patients + 1, 1, n_patients * max_all_t + 1);
}

parameters {
  #include "baseline_hazard_parameters.stan"
  #include "crcr/crcr_parameters.stan"
  #include "recruit/recruit_parameters.stan"
  
  vector[n_tumor_covar] tumor_stim_pop_coef;
  vector[n_covar] covar_effect;  
  real conf_resp_effect;
  
  vector<lower = 0>[add_trial_level ? n_tumor_covar + n_covar + 1 : 0] covar_trial_sd;
  // cholesky_factor_corr[add_trial_level ? n_tumor_covar + n_covar + 1 : 0] L_covar_trial_corr;
  matrix[n_tumor_covar + n_covar + 1, add_trial_level ? n_trials : 0] raw_covar_trial_coef;
}


transformed parameters {
  #include "baseline_hazard_transformed_parameters.stan"
  #include "crcr/crcr_transformed_parameters.stan"
  
  matrix<upper = 0>[n_time_periods, n_causes] log_cond_prob_surv;
  matrix[n_patients, n_causes] time_invariant_log_hazard_ratio = rep_matrix(tumor_sum_covar * tumor_stim_pop_coef + covar_design_matrix * covar_effect, n_causes);
  
  time_invariant_log_hazard_ratio[, n_causes] += conf_resp_effect; 
  
  matrix[n_tumor_covar + n_covar + 1, add_trial_level ? n_trials : 0] covar_trial_coef;
 
  if (add_trial_level) {
    covar_trial_coef = diag_pre_multiply(covar_trial_sd, raw_covar_trial_coef); 
    // covar_trial_coef = diag_pre_multiply(covar_trial_sd, L_covar_trial_corr) * raw_covar_trial_coef; 
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
}

model {
  // Priors
  
  profile("crcr priors") {
    #include "crcr/crcr_priors.stan"
  }
  
  #include "baseline_hazard_priors.stan"
  #include "recruit/recruit_priors.stan"
  
  profile("pfs priors") {
    tumor_stim_pop_coef ~ normal(0, tumor_stim_pop_coef_sd);
    conf_resp_effect ~ normal(0, conf_resp_effect_sd); 
    covar_effect ~ normal(0, covar_effect_sd);
    
    if (add_trial_level) {
      covar_trial_sd ~ normal(0, covar_trial_sd_sd);
      // L_covar_trial_corr ~ lkj_corr_cholesky(covar_trial_corr_eta);
      to_vector(raw_covar_trial_coef) ~ std_normal(); 
    }
  }
  
  // Likelihood
  
  if (fit_data) {
    profile("crcr loglik") {
      target += reduce_sum(
        partial_sum_crcr_lupmf, last_unclassified_response_week, crcr_grain_size,
        confirmed_response_cause, early_confirmed_response_censored, log_crcr_cond_prob_surv, max_confresp_week
      );
    }
    
    experiment_start_week ~ neg_binomial_2(recruit_lambda[patient_trial], recruit_phi[patient_trial]);  
  
    profile("pfs loglik") {  
      matrix[n_patients, n_causes] response_lp = append_col( 
        calc_pch_loglik2(
          pfs, right_censored, interval_censored, ignore_interval_censoring, log_cond_prob_surv[, 1], max_all_t, rep_array(1, n_patients)
        ),
        calc_pch_loglik2(
          pfs, right_censored, interval_censored, ignore_interval_censoring, log_cond_prob_surv[, 2], max_all_t, rep_array(1, n_patients)
        )
      );
      
      for (s in 1:n_trials) {
        int patient_pos = trial_patient_pos[s];
        int patient_end = trial_patient_pos[s + 1] - 1; 
        
        if (s != leave_out_trial) {
          for (i in patient_pos:patient_end) {
            if (confirmed_response_censored[i]) { // Unclassified
              int conf_resp_interval_pos = patient_conf_resp_interval_pos[i];
              int conf_resp_interval_end = patient_conf_resp_interval_pos[i + 1] - 1;
            
              real prob_non_response = calc_cif(1, log_crcr_cond_prob_surv[conf_resp_interval_pos:conf_resp_interval_end], max_confresp_week).2[1, 1];
            
              target += log_mix(prob_non_response, response_lp[i, 1], response_lp[i, 2]);
            } else {
              target += response_lp[i, confirmed_response_cause[i]];
            }
          }
        }
      }
    }
  }
}

generated quantities {
  #include "crcr/crcr_gen_quants.stan"
  
  array[n_patients] int<lower = 0, upper = 1> sim_confirmed_response = confirmed_response;
  
  sim_confirmed_response[missing_confirmed_response] = bernoulli_rng(prob_cause[missing_confirmed_response, 2]); 
  
  array[n_patients] int<lower = 0> sim_pfs; 
  array[n_patients] int<lower = 0, upper = 1> sim_censored; 
  
  real<lower = 0> sim_median_pfs;
  vector<lower = 0>[n_trials] sim_trial_median_pfs;
  
  vector<lower = 0, upper = 1>[max_all_t + 1] km_est; 
  array[n_trials] vector<lower = 0, upper = 1>[max_all_t + 1] trial_km_est; 
  
  array[n_bootstrap_cr_maturity_rates] int<lower = 1> bs_prediction_calendar_week = rep_array(1000000, n_bootstrap_cr_maturity_rates);
  array[n_bootstrap_cr_maturity_rates] int n_bs_sample_classified = rep_array(0, n_bootstrap_cr_maturity_rates); 
  array[n_bootstrap_cr_maturity_rates] int n_bs_sample_unclassified = rep_array(0, n_bootstrap_cr_maturity_rates); 
  vector<lower = 0>[n_bootstrap_cr_maturity_rates] bs_median_pfs = rep_vector(max_all_t, n_bootstrap_cr_maturity_rates);
 
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
    
    array[n_bootstrap_sample] int bs_sample_idx;
    array[n_bootstrap_sample] int bs_start_calendar_week; 
    array[n_bootstrap_sample] int bs_cr_mature_calendar_week; 
    array[n_bootstrap_sample] int bs_cr_censored; 
    vector[n_bootstrap_sample] bs_cr_maturity_rate;
    
    profile("bootstrap") {
      if (leave_out_trial > 0 && n_bootstrap_sample > 0 && n_bootstrap_cr_maturity_rates > 0) {
        int patient_pos = trial_patient_pos[leave_out_trial];
        int patient_end = trial_patient_pos[leave_out_trial + 1] - 1;
        
        bs_sample_idx = discrete_range_rng(rep_array(patient_pos, n_bootstrap_sample), rep_array(patient_end, n_bootstrap_sample)); 
        
        for (bsi in 1:n_bootstrap_sample) {
          bs_cr_mature_calendar_week[bsi] = sorted_experiment_start_week[patient_pos + bsi - 1] + confirmed_response_week[bs_sample_idx[bsi]] - 1;
        }
        
        array[n_bootstrap_sample] int cr_mature_sorted_idx = sort_indices_asc(bs_cr_mature_calendar_week);
        
        bs_sample_idx = bs_sample_idx[cr_mature_sorted_idx]; 
        bs_cr_mature_calendar_week = bs_cr_mature_calendar_week[cr_mature_sorted_idx];
        bs_cr_censored = confirmed_response_censored[bs_sample_idx];
        bs_start_calendar_week = sorted_experiment_start_week[patient_pos:(patient_pos + n_bootstrap_sample - 1)][cr_mature_sorted_idx];
        
        int mature_rate_pos = 1;
        
        for (bsi in 1:n_bootstrap_sample) {
          bs_cr_maturity_rate[bsi] = 1.0 * (1 - bs_cr_censored[bsi]) / n_bootstrap_sample; 
          
          if (bsi > 1) {
            bs_cr_maturity_rate[bsi] += bs_cr_maturity_rate[bsi - 1];
          }
          
          if (mature_rate_pos <= n_bootstrap_cr_maturity_rates && bs_cr_maturity_rate[bsi] >= bootstrap_cr_maturity_rates[mature_rate_pos]) {
            n_bs_sample_classified[mature_rate_pos] = bsi;
            bs_prediction_calendar_week[mature_rate_pos] = bs_cr_mature_calendar_week[bsi]; 
            
            int n_remaining = n_bootstrap_sample - bsi;
            array[n_remaining] int remaining_patients, remaining_sort_idx;
            
            if (n_remaining > 0) {
              remaining_sort_idx = sort_indices_asc(bs_start_calendar_week[(bsi + 1):]);
              remaining_patients = bs_sample_idx[(bsi + 1):][remaining_sort_idx];
            }
            
            for (ri in 1:n_remaining) {
              if (bs_start_calendar_week[(bsi + 1):][remaining_sort_idx[ri]] <= bs_prediction_calendar_week[mature_rate_pos]) {
                n_bs_sample_unclassified[mature_rate_pos] += 1; 
              } else {
                break;
              }
            }
            
            int n_bs_sample = n_bs_sample_classified[mature_rate_pos] + n_bs_sample_unclassified[mature_rate_pos];
            array[n_bs_sample] int bs_pfs;
            
            bs_pfs[:n_bs_sample_classified[mature_rate_pos]] = forecast_pfs_rng(
              bs_prediction_calendar_week[mature_rate_pos], 
              bs_start_calendar_week[:n_bs_sample_classified[mature_rate_pos]],
              bs_sample_idx[:n_bs_sample_classified[mature_rate_pos]],
              patient_pfs_interval_pos, pfs, right_censored, log_cond_prob_surv,
              patient_conf_resp_interval_pos, confirmed_response_week, confirmed_response, confirmed_response_censored, log_crcr_cond_prob_surv
            );
            
            bs_pfs[(n_bs_sample_classified[mature_rate_pos] + 1):] = forecast_pfs_rng(
              bs_prediction_calendar_week[mature_rate_pos], 
              bs_start_calendar_week[(bsi + 1):][remaining_sort_idx[:n_bs_sample_unclassified[mature_rate_pos]]],
              remaining_patients[:n_bs_sample_unclassified[mature_rate_pos]],
              patient_pfs_interval_pos, pfs, right_censored, log_cond_prob_surv,
              patient_conf_resp_interval_pos, confirmed_response_week, confirmed_response, confirmed_response_censored, log_crcr_cond_prob_surv
            );
            
            bs_median_pfs[mature_rate_pos] = survival_median(bs_pfs, max_all_t).1;
            
            mature_rate_pos += 1;
          }
        }
      }
    }
  }
}
