functions {
  #include "../util.stan"
  #include "../pos.stan"
  #include "../gp.stan"
  #include "../pfs_functions.stan"
  #include "../lfo.stan"
  #include "sf-ssls_functions.stan"
  #include "recist.stanfunctions"
}  

data {
  #include "../base_data.stan"
  #include "base_data.stan"
  
  #include "sf-ssls-outcomes_info.stan"
  
  #include "sf-ssls-hyperparam.stan" 
} 

transformed data {
  #include "../base_transformed_data.stan"
  #include "tumor_transformed_data.stan"
  #include "sf-transformed_data.stan"
  #include "other_events_transformed_data.stan"
  #include "mature_cutoffs_transformed_data.stan"
  #include "sf-ssls-outcomes_info_transformed_data.stan"
}

parameters {
  #include "other_events_parameters.stan"
  #include "sf-ssls-parameters.stan"
}

transformed parameters {
  #include "other_events_transformed_parameters.stan"
  #include "sf-ssls-transformed_parameters.stan"
  
  matrix[n_train_patients, n_causes] patient_response_lp = rep_matrix(0, n_train_patients, n_causes); 

  // Non-target progression
  patient_response_lp[, 1] = calc_pch_loglik(
    non_target_pfs[train_patients_pos:train_patients_end], 
    non_target_right_censored[train_patients_pos:train_patients_end], 
    zeros_int_array(n_train_patients),
    0, 
    log_cond_prob_surv[1]
  );
}

model {
  #include "other_events_priors.stan"
  #include "sf-ssls-priors.stan"

  profile("loglik") { 
    if (fit_tumor_data) {
      for (i in train_patients_pos:train_patients_end) {
        int train_idx = i - train_patients_pos + 1;
        int train_visit_start, train_visit_end;
        (train_visit_start, train_visit_end) = get_pos(train_patient_visit_pos, train_idx);
        int visit_pos, visit_end;
        (visit_pos, visit_end) = get_pos(patient_visit_pos, i);
        normalized_sld[visit_pos:visit_end] ~ sf_log_space_obs(states[train_visit_start:train_visit_end], measure_sd, log_lod - log(sum_tumor_size[visit_pos]));
      }

      for (s in 1:n_trials) {
        for (k in 1:n_causes) {
          target += sum(get_sub_vector(patient_response_lp[, k], train_trial_patient_pos, s));
        }
      }
    }
  }
}

generated quantities {
  // Parameters //////////////////////////////////////////////
  
  real pop_log_growth_rate = pop_log_net_rate - log_diff_exp(pop_log_rate_ratio, 0);
  vector[n_trials] trial_log_growth_rate = pop_log_growth_rate + trial_log_net_rate_effect;
  vector[n_trials] trial_log_growth_rate_residual = trial_log_growth_rate - pop_log_growth_rate;
  vector[n_train_patients] patient_log_growth_rate_residual = patient_log_growth_rate - trial_log_growth_rate[patient_trial[train_patients_pos:train_patients_end]];
  
  real pop_log_decrease_rate = pop_log_growth_rate + pop_log_rate_ratio;
  vector[n_trials] trial_log_decrease_rate = trial_log_growth_rate + pop_log_rate_ratio; 
  vector[n_trials] trial_log_decrease_rate_residual = trial_log_decrease_rate - pop_log_decrease_rate;
  vector[n_train_patients] patient_log_decrease_rate_residual = patient_log_decrease_rate - trial_log_decrease_rate[patient_trial[train_patients_pos:train_patients_end]];
  
  vector[n_train_patients] patient_decrease_prop_residual = inv_logit(patient_decrease_prop_logis) - inv_logit(trial_decrease_prop_logis[patient_trial[train_patients_pos:train_patients_end]]);
  vector[n_trials] trial_decrease_prop_residual = inv_logit(trial_decrease_prop_logis) - inv_logit(pop_decrease_prop_logis);
  real pop_log_decrease_prop = -log1p_exp(-pop_decrease_prop_logis);  
  real pop_log_growth_prop = pop_log_decrease_prop - pop_decrease_prop_logis;
  vector[n_trials] trial_log_decrease_prop = -log1p_exp(-trial_decrease_prop_logis);
  vector[n_trials] trial_log_growth_prop = trial_log_decrease_prop - trial_decrease_prop_logis;
  
  corr_matrix[independ_cross_process_noise ? 0 : 2] process_corr;
  
  if (!independ_cross_process_noise) {
    process_corr = L_process_corr * L_process_corr';
  }
  
  vector<lower = 0, upper = 1>[max_all_t] all_growth_factor = get_growth_lag_factor(all_tumor_measure_t, exp(pop_log_growth_lag), exp(pop_log_growth_transition_rate));
  
  matrix[max_all_t, 2] all_scaled_process_sd = scale_process_sd(all_tumor_measure_t, pop_process_sd);

  // Latent states //////////////////////////////////////////
  
  matrix[n_total_train_visits_m1, 2] obs_patient_process_noise;
  matrix[n_total_train_forecast_visits, 2] forecast_patient_process_noise;
  
  matrix[n_total_train_forecast_visits, 2] forecast_patient_states;
  
  vector[n_total_train_visits] rep_patient_log_sld;
  vector[n_total_train_forecast_visits] forecast_patient_log_sld;
  
  array[n_total_train_visits] int<lower = CR, upper = PD + 1> rep_recist = rep_array(PD + 1, n_total_train_visits);
  array[n_total_train_forecast_visits] int<lower = CR, upper = PD> forecast_recist;
   
  // Endpoints (PFS, ORR, Median PFS, PFSn, ...) ////////////
 
  array[n_train_patients] int<lower = 0> sample_target_pfs, spop_target_pfs, spop_non_target_pfs, spop_pfs, spop_target_obs_cens_pfs; 
  array[n_train_patients] int<lower = 0, upper = 1> 
    sample_target_right_censored, spop_target_right_censored, spop_non_target_right_censored, spop_right_censored, spop_target_obs_cens_right_censored; 

  // Forecasting for right censored patients 
  array[n_train_right_censored_patients] int<lower = 0> forecast_target_pfs, forecast_non_target_pfs, forecast_pfs; 
  array[n_train_right_censored_patients] int<lower = 0, upper = 1> 
    forecast_target_right_censored, forecast_non_target_right_censored, forecast_right_censored; 
  
  array[n_trials] vector<lower = 0, upper = 1>[max_all_t + 1] sample_km_est, // sample_target_km_est, sample_non_target_km_est,
                                                              spop_target_km_est, spop_non_target_km_est, spop_km_est, spop_target_obs_cens_km_est;

  array[n_cond_group] vector<lower = 0, upper = 1>[max_all_t + 1] cond_sample_km_est, 
    cond_spop_target_km_est, cond_spop_non_target_km_est, cond_spop_km_est, cond_spop_target_obs_cens_km_est;

  array[n_trials] vector<lower = 0, upper = 1>[n_pfs_timepoints] forecast_target_pfs_n; //, sample_pfs_n; 
  array[n_cond_group] vector<lower = 0, upper = 1>[n_pfs_timepoints] cond_forecast_target_pfs_n; //, cond_sample_pfs_n;

  array[n_train_patients] int<lower = 0, upper = 1> forecast_confirmed_response;
  vector<lower = 0, upper = 1>[n_trials] forecast_target_orr;
  vector<lower = 0, upper = 1>[n_cond_group] cond_forecast_target_orr;

  vector<lower = 0>[n_trials] trial_median_pfs = zeros_vector(n_trials);
  vector<lower = 0>[n_cond_group] cond_median_pfs = zeros_vector(n_cond_group);

  // Cutoffs at different maturity times //////////////////////////////////////////////
  
  array[n_trials, n_mature_cutoffs_calendar_days] vector<lower = 0, upper = 1>[truncated_max_all_t + 1] 
    cutoff_sample_km_est, cutoff_spop_target_km_est; //, cutoff_spop_non_target_km_est, cutoff_spop_km_est, cutoff_spop_target_obs_cens_km_est;

  array[n_cond_group, n_mature_cutoffs_calendar_days] vector<lower = 0, upper = 1>[truncated_max_all_t + 1] 
    cutoff_cond_spop_target_km_est;
   
  array[n_trials, n_mature_cutoffs_calendar_days] vector<lower = 0, upper = 1>[n_pfs_timepoints] cutoff_forecast_target_pfs_n;
  array[n_cond_group, n_mature_cutoffs_calendar_days] vector<lower = 0, upper = 1>[n_pfs_timepoints] cutoff_cond_forecast_target_pfs_n;

  array[n_train_patients, n_mature_cutoffs_calendar_days] int<lower = 0, upper = 1> cutoff_forecast_confirmed_response;
  matrix<lower = 0, upper = 1>[n_trials, n_mature_cutoffs_calendar_days] cutoff_forecast_target_orr;
  matrix<lower = 0, upper = 1>[n_cond_group, n_mature_cutoffs_calendar_days] cutoff_cond_forecast_target_orr;

  matrix<lower = 0>[n_trials, n_mature_cutoffs_calendar_days] cutoff_trial_median_pfs = rep_matrix(0, n_trials, n_mature_cutoffs_calendar_days);
  matrix<lower = 0>[n_cond_group, n_mature_cutoffs_calendar_days] cutoff_cond_median_pfs = rep_matrix(0, n_cond_group, n_mature_cutoffs_calendar_days);
  
  {
    int right_censored_idx = 1;
    
    for (i in train_patients_pos:train_patients_end) {
      int train_idx = i - train_patients_pos + 1;
  
      int train_visit_start, train_screening_visit_end, train_treat_visit_start, train_visit_end;
      (train_visit_start, train_screening_visit_end, train_treat_visit_start, train_visit_end) = get_visit_pos(
        train_patient_visit_pos, train_idx, n_patient_screening_visits[i]);
      int train_visit_size = get_pos_size(train_patient_visit_pos, train_idx);
      int n_obs_treat_visits = train_visit_size - n_patient_screening_visits[i];
  
      int train_visit_m1_start,  train_visit_m1_end;
      (train_visit_m1_start, train_visit_m1_end) = get_pos(train_patient_visit_m1_pos, train_idx);
      int train_visit_m1_size = get_pos_size(train_patient_visit_m1_pos, train_idx);
  
      int train_forecast_visit_start, train_forecast_visit_end;
      (train_forecast_visit_start, train_forecast_visit_end) = get_pos(train_forecast_visits_pos, train_idx);
  
      int visit_pos, visit_screening_end, visit_treat_pos, visit_end;
      (visit_pos, visit_screening_end, visit_treat_pos, visit_end) = get_visit_pos(patient_visit_pos, i, n_patient_screening_visits[i]);

      int visit_m1_start, visit_m1_end;
      (visit_m1_start, visit_m1_end) = get_pos(patient_visit_m1_pos, i);
  
      int forecast_visit_start, forecast_visit_end;
      (forecast_visit_start, forecast_visit_end) = get_pos(forecast_visits_pos, i);
      int forecast_size = get_pos_size(forecast_visits_pos, i);
      
      array[n_patient_visits[i]] int curr_visits = get_int_sub_array(train_patient_visits, train_patient_visit_pos, train_idx);
      
      array[n_patient_forecast_visits[i] + 1] int forecast_time = linspaced_int_array(
        n_patient_forecast_visits[i] + 1, patient_last_obs_visit[i], last_predict_visit);
      
      // Generate patient states using direct tuple assignment
      // Generate zeros for observed process noise (as it was hardcoded before)
      obs_patient_process_noise[train_visit_m1_start:train_visit_m1_end] = rep_matrix(0.0, n_patient_visits[i] - 1, 2);
      
      // Generate states, SLD replications, and forecasts
      (forecast_patient_states[train_forecast_visit_start:train_forecast_visit_end], 
       rep_patient_log_sld[train_visit_start:train_visit_end], 
       forecast_patient_log_sld[train_forecast_visit_start:train_forecast_visit_end]) = 
        generate_patient_states_rng(
          states[train_visit_start:train_visit_end],
          forecast_time,
          patient_log_decrease_rate[train_idx], patient_log_growth_rate[train_idx],
          sum_tumor_size[visit_pos], 
          0.0001, 0.0001, // exp(patient_log_growth_lag[train_idx]), exp(pop_log_growth_transition_rate),
          rep_matrix(0.0, forecast_size, 2), // Hardcode zeros for forecast process noise
          measure_sd
        );
      
      // Generate forecast process noise conditionally (as it was in the original code)
      if (independ_long_process_noise) {
        forecast_patient_process_noise[train_forecast_visit_start:train_forecast_visit_end] = multi_normal_rng(
          forecast_size, pop_process_sd, independ_cross_process_noise ? diag_matrix(ones_vector(2)) : L_process_corr
        );
      } else {
        forecast_patient_process_noise[train_forecast_visit_start:train_forecast_visit_end] = multi_normal_rng(
          obs_patient_process_noise[train_visit_m1_start:train_visit_m1_end],
          get_int_sub_array(t_patient_visits, patient_visit_pos, i)[2:],
          forecast_time[2:],
          exp(log_pop_tumor_gp_rho),
          pop_process_sd,
          independ_cross_process_noise ? diag_matrix(ones_vector(2)) : L_process_corr,
          delta
        );
      }      
      
      (spop_non_target_pfs[train_idx], spop_non_target_right_censored[train_idx]) = survival_time_rng(log_cond_prob_surv[1, train_idx]);
      
      array[n_obs_treat_visits + forecast_size] int full_predict_recist = calculate_target_recist(
        exp(append_row(rep_patient_log_sld[train_visit_start:train_visit_end], 
                      forecast_patient_log_sld[train_forecast_visit_start:train_forecast_visit_end])) * 10,
        n_patient_screening_visits[i]
      );
      
      rep_recist[train_treat_visit_start:train_visit_end] = full_predict_recist[:n_obs_treat_visits]; 
      forecast_recist[train_forecast_visit_start:train_forecast_visit_end] = full_predict_recist[(n_obs_treat_visits + 1):];
     
      int forecast_confirmed_response_week = 0, confirmed_response_censored = 0; 
      array[n_mature_cutoffs_calendar_days] int cutoff_forecast_confirmed_response_week = zeros_int_array(n_mature_cutoffs_calendar_days),
                                                cutoff_confirmed_response_censored = zeros_int_array(n_mature_cutoffs_calendar_days);

      array[train_visit_size + forecast_size] int curr_recist;
      curr_recist[:train_visit_size] = recist[visit_pos:visit_end];
          
      if (n_patient_forecast_visits[i] > 0) {
        curr_recist[(train_visit_size + 1):] = forecast_recist[train_forecast_visit_start:train_forecast_visit_end];

        (sample_target_pfs[train_idx], sample_target_right_censored[train_idx]) = 
          find_first_week(curr_recist, PD, n_patient_screening_visits[i], curr_visits, forecast_time, max_all_t);
        (forecast_confirmed_response_week, confirmed_response_censored) = 
          find_first_week(curr_recist, { PR, CR }, 2, n_patient_screening_visits[i], curr_visits, forecast_time, max_all_t);

        for (m in 1:n_mature_cutoffs_calendar_days) {
          int last_after_cutoff_visit = min(patient_relative_day_at_cutoff[i, m], train_visit_size + forecast_size);
          (cutoff_forecast_confirmed_response_week[m], cutoff_confirmed_response_censored[m]) = find_first_week(
            curr_recist[:last_after_cutoff_visit], { PR, CR }, 2, n_patient_screening_visits[i], curr_visits, forecast_time, max_all_t);
        }
         
        if (right_censored[i]) {
          (forecast_non_target_pfs[right_censored_idx], forecast_non_target_right_censored[right_censored_idx]) = survival_time_rng(
            log_cond_prob_surv[1, train_idx], pfs[i] + interval_censored[i], right_censored[i], 0
          );

          (forecast_target_pfs[right_censored_idx], forecast_target_right_censored[right_censored_idx]) = 
            find_first_week(forecast_recist[train_forecast_visit_start:train_forecast_visit_end], PD, 
                            n_patient_screening_visits[i], curr_visits, forecast_time, max_all_t);
          
          forecast_pfs[right_censored_idx] = min(forecast_target_pfs[right_censored_idx], forecast_non_target_pfs[right_censored_idx]);
          forecast_right_censored[right_censored_idx] = forecast_target_right_censored[right_censored_idx] && forecast_non_target_right_censored[right_censored_idx];
          
          right_censored_idx += 1;
        }
        
        (spop_target_pfs[train_idx], spop_target_right_censored[train_idx]) = find_first_week(
          append_array(rep_recist[train_treat_visit_start:train_visit_end], 
                       forecast_recist[train_forecast_visit_start:train_forecast_visit_end]), 
          PD, n_patient_screening_visits[i], curr_visits, forecast_time, max_all_t); 
      } else {
        (spop_target_pfs[train_idx], spop_target_right_censored[train_idx]) = find_first_week(
          rep_recist[train_treat_visit_start:train_visit_end], PD, n_patient_screening_visits[i], curr_visits, forecast_time, max_all_t); 
        (sample_target_pfs[train_idx], sample_target_right_censored[train_idx]) = find_first_week(
          recist[visit_pos:visit_end], PD, n_patient_screening_visits[i], curr_visits, forecast_time, max_all_t);
        (forecast_confirmed_response_week, confirmed_response_censored) = find_first_week(
          recist[visit_pos:visit_end], { PR, CR }, 2, n_patient_screening_visits[i], curr_visits, forecast_time, max_all_t);

        for (m in 1:n_mature_cutoffs_calendar_days) {
          int last_after_cutoff_visit = min(patient_relative_day_at_cutoff[i, m], train_visit_size);
          (cutoff_forecast_confirmed_response_week[m], cutoff_confirmed_response_censored[m]) = find_first_week(
            recist[:last_after_cutoff_visit], { PR, CR }, 2, n_patient_screening_visits[i], curr_visits, forecast_time, max_all_t);
        }
      }
     
      spop_target_obs_cens_right_censored[train_idx] = right_censored[i]; 
      spop_target_obs_cens_pfs[train_idx] = right_censored[i] ? min(spop_target_pfs[train_idx], pfs[i]) : spop_target_pfs[train_idx];
      
      spop_pfs[train_idx] = min(spop_non_target_pfs[train_idx] + 1, 
                                max(0, spop_target_pfs[train_idx])); // BUG a couple of patients end up with negative weeks. We need to figure out why.
      spop_right_censored[train_idx] = spop_target_right_censored[train_idx] && spop_non_target_right_censored[train_idx]; 

      forecast_confirmed_response[train_idx] = 
        !confirmed_response_censored && forecast_confirmed_response_week < sample_target_pfs[train_idx];

      for (m in 1:n_mature_cutoffs_calendar_days) {
        cutoff_forecast_confirmed_response[train_idx, m] = 
          !cutoff_confirmed_response_censored[m] && cutoff_forecast_confirmed_response_week[m] < sample_target_pfs[train_idx];
      }
    }
      
    for (s in 1:n_trials) {
      if (get_pos_size(train_trial_patient_pos, s) > 0) {
        forecast_target_orr[s] = mean(get_int_sub_array(forecast_confirmed_response, train_trial_patient_pos, s));

        int n_curr_uncensored_obs = get_pos_size(train_trial_right_uncensored_pos, s);
        array[n_curr_uncensored_obs] int curr_uncensored_obs = get_int_sub_array(train_right_uncensored_patients, train_trial_right_uncensored_pos, s);
        int n_curr_right_censored = get_pos_size(train_trial_right_censored_pos, s);

        array[n_curr_uncensored_obs + n_curr_right_censored] int curr_sample_pfs = 
          get_int_sub_array(sample_target_pfs, train_trial_patient_pos, s),
                                                                 curr_sample_right_censored = 
          get_int_sub_array(sample_target_right_censored, train_trial_patient_pos, s);

        sample_km_est[s] = estimate_kaplan_meier(curr_sample_pfs, curr_sample_right_censored, max_all_t).1; 

        array[n_curr_uncensored_obs + n_curr_right_censored] int curr_spop_target_pfs = 
          get_int_sub_array(spop_target_pfs, train_trial_patient_pos, s),
                                                                 curr_spop_target_right_censored = 
          get_int_sub_array(spop_target_right_censored, train_trial_patient_pos, s); 
        
        spop_target_km_est[s] = estimate_kaplan_meier(curr_spop_target_pfs, curr_spop_target_right_censored, max_all_t, 0).1; 
                                               
        spop_target_obs_cens_km_est[s] = estimate_kaplan_meier(get_int_sub_array(spop_target_obs_cens_pfs, train_trial_patient_pos, s), 
                                               get_int_sub_array(spop_target_obs_cens_right_censored, train_trial_patient_pos, s), 
                                               max_all_t, 0).1; 
                                               
        spop_non_target_km_est[s] = estimate_kaplan_meier(get_int_sub_array(spop_non_target_pfs, train_trial_patient_pos, s), 
                                               get_int_sub_array(spop_non_target_right_censored, train_trial_patient_pos, s), 
                                               max_all_t, 0).1; 
        
        spop_km_est[s] = estimate_kaplan_meier(get_int_sub_array(spop_pfs, train_trial_patient_pos, s), 
                                               get_int_sub_array(spop_right_censored, train_trial_patient_pos, s), 
                                               max_all_t, 0).1;
                                               
        trial_median_pfs[s] = km_median(sample_km_est[s]).1;

        for (n in 1:n_pfs_timepoints) {
          forecast_target_pfs_n[s, n] = calc_km_pfs_n(sample_km_est[s], months_to_weeks(pfs_timepoints[n])); 
        }

        for (m in 1:n_mature_cutoffs_calendar_days) {
          cutoff_forecast_target_orr[s, m] = mean(get_int_sub_array(cutoff_forecast_confirmed_response[, m], train_trial_patient_pos, s));

          array[n_curr_uncensored_obs + n_curr_right_censored] int curr_truncated_sample_pfs, curr_truncated_sample_right_censored;
          (curr_truncated_sample_pfs, curr_truncated_sample_right_censored) = truncate_at_max_time(
            curr_sample_pfs, curr_sample_right_censored, patient_relative_week_at_cutoff[, m] 
          ); 

          cutoff_sample_km_est[s, m] = estimate_kaplan_meier(curr_truncated_sample_pfs, curr_truncated_sample_right_censored, truncated_max_all_t).1; 

          array[n_curr_uncensored_obs + n_curr_right_censored] int curr_truncated_spop_target_pfs, curr_truncated_spop_target_right_censored;
          (curr_truncated_spop_target_pfs, curr_truncated_spop_target_right_censored) = truncate_at_max_time(
            curr_spop_target_pfs, curr_spop_target_right_censored, patient_relative_week_at_cutoff[, m] 
          ); 

          cutoff_spop_target_km_est[s, m] = estimate_kaplan_meier(
            curr_truncated_spop_target_pfs, curr_truncated_spop_target_right_censored, truncated_max_all_t).1; 

          cutoff_trial_median_pfs[s, m] = km_median(cutoff_sample_km_est[s, m]).1;

          for (n in 1:n_pfs_timepoints) {
            cutoff_forecast_target_pfs_n[s, m, n] = calc_km_pfs_n(cutoff_sample_km_est[s, m], months_to_weeks(pfs_timepoints[n])); 
          }
        }
      } else {
        spop_target_km_est[s] = zeros_vector(max_all_t + 1);
        spop_target_obs_cens_km_est[s] = zeros_vector(max_all_t + 1);
        spop_non_target_km_est[s] = zeros_vector(max_all_t + 1);
        spop_km_est[s] = zeros_vector(max_all_t + 1);
        sample_km_est[s] = zeros_vector(max_all_t + 1);
      }
    } 
    
    for (c in 1:n_cond_group) {
      array[cond_group_size[c]] int curr_group_patients = get_int_sub_array(cond_group, cond_group_pos, c);
      
      cond_forecast_target_orr[c] = mean(forecast_confirmed_response[curr_group_patients]);

      cond_sample_km_est[c] = estimate_kaplan_meier(sample_target_pfs[curr_group_patients], 
                                             sample_target_right_censored[curr_group_patients],
                                             max_all_t, 0).1;
    
      cond_spop_target_km_est[c] = estimate_kaplan_meier(spop_target_pfs[curr_group_patients], 
                                             spop_target_right_censored[curr_group_patients],
                                             max_all_t, 0).1;
                                             
      cond_spop_target_obs_cens_km_est[c] = estimate_kaplan_meier(spop_target_obs_cens_pfs[curr_group_patients],
                                             spop_target_obs_cens_right_censored[curr_group_patients],
                                             max_all_t, 0).1; 
                                             
      cond_spop_non_target_km_est[c] = estimate_kaplan_meier(spop_non_target_pfs[curr_group_patients],
                                             spop_non_target_right_censored[curr_group_patients],
                                             max_all_t, 0).1; 
      
      cond_spop_km_est[c] = estimate_kaplan_meier(spop_pfs[curr_group_patients],
                                             spop_right_censored[curr_group_patients],
                                             max_all_t, 0).1; 

      cond_median_pfs[c] = km_median(cond_sample_km_est[c]).1;

      for (n in 1:n_pfs_timepoints) {
        cond_forecast_target_pfs_n[c, n] = calc_km_pfs_n(cond_sample_km_est[c], months_to_weeks(pfs_timepoints[n])); 
      }

      for (m in 1:n_mature_cutoffs_calendar_days) {
        cutoff_cond_forecast_target_orr[c, m] = mean(cutoff_forecast_confirmed_response[curr_group_patients, m]);

        array[cond_group_size[c]] int curr_truncated_cond_spop_target_pfs, curr_truncated_conf_spop_target_right_censored;
        (curr_truncated_cond_spop_target_pfs, curr_truncated_conf_spop_target_right_censored) = truncate_at_max_time(
          spop_target_pfs[curr_group_patients], spop_target_right_censored[curr_group_patients], 
          patient_relative_week_at_cutoff[, m] 
        );

        cutoff_cond_spop_target_km_est[c, m] = estimate_kaplan_meier(
          curr_truncated_cond_spop_target_pfs, curr_truncated_conf_spop_target_right_censored, truncated_max_all_t).1;

        cutoff_cond_median_pfs[c, m] = km_median(cutoff_cond_spop_target_km_est[c, m]).1;

        for (n in 1:n_pfs_timepoints) {
          cutoff_cond_forecast_target_pfs_n[c, m, n] = calc_km_pfs_n(cutoff_cond_spop_target_km_est[c, m], months_to_weeks(pfs_timepoints[n])); 
        }
      }
    }
  }
  
  #include "sf-ssls-accuracy_gen_quant.stan"
}

