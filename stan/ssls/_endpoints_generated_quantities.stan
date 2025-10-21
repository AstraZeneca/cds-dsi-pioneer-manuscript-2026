
  real pop_log_decrease_frac = log_inv_logit(frac_logit_loc_pop);
  real pop_log_growth_frac   = log1m_inv_logit(frac_logit_loc_pop);
  real pop_log_decrease_rate = tr_loc_pop + pop_log_decrease_frac;
  real pop_log_growth_rate   = tr_loc_pop + pop_log_growth_frac;
  vector[n_trials] trial_log_total_rate = tr_loc_pop + tr_effect_trial_intercept;
  vector[n_trials] trial_log_decrease_rate = trial_log_total_rate + pop_log_decrease_frac;
  vector[n_trials] trial_log_growth_rate   = trial_log_total_rate + pop_log_growth_frac;
  vector[n_trials] trial_log_growth_rate_residual = trial_log_growth_rate - pop_log_growth_rate;
  vector[n_patients] patient_log_growth_rate_residual = patient_log_growth_rate - trial_log_growth_rate[patient_trial];
  vector[n_trials] trial_log_decrease_rate_residual = trial_log_decrease_rate - pop_log_decrease_rate;
  vector[n_patients] patient_log_decrease_rate_residual = patient_log_decrease_rate - trial_log_decrease_rate[patient_trial];

  vector<lower = 0, upper = 1>[max_all_t] all_growth_factor = get_growth_lag_factor(all_tumor_measure_t, exp(pop_log_growth_lag), exp(pop_log_growth_transition_rate));
  
  matrix[max_all_t, 2] all_scaled_process_sd = scale_process_sd(all_tumor_measure_t, pop_process_sd);

  // Latent states //////////////////////////////////////////
  
  matrix[n_total_visits_m1, 2] obs_patient_process_noise;
  matrix[n_total_forecast_visits, 2] forecast_patient_process_noise;
  
  matrix[n_total_forecast_visits, 2] forecast_patient_states;
  
  vector[sum(n_patient_visits)] mean_patient_log_sld, rep_patient_log_sld;
  vector[n_total_forecast_visits] forecast_mean_patient_log_sld, forecast_patient_log_sld;
  
  array[sum(n_patient_visits)] int<lower = CR, upper = PD + 1> rep_recist = rep_array(PD + 1, sum(n_patient_visits));
  array[n_total_forecast_visits] int<lower = CR, upper = PD> forecast_recist;
   
  // Endpoints (PFS, ORR, Median PFS, PFSn, ...) ////////////
 
  array[n_patients] int<lower = 0> sample_target_pfs, spop_target_pfs, spop_target_obs_cens_pfs; // spop_non_target_pfs, spop_pfs,  
  array[n_patients] int<lower = 0, upper = 1> 
    sample_target_right_censored, spop_target_right_censored, spop_target_obs_cens_right_censored; // spop_non_target_right_censored, spop_right_censored; 

  // Forecasting for right censored patients 
  array[n_right_censored_patients] int<lower = 0> forecast_target_pfs; //, forecast_non_target_pfs, forecast_pfs; 
  array[n_right_censored_patients] int<lower = 0, upper = 1> 
    forecast_target_right_censored; //, forecast_non_target_right_censored, forecast_right_censored; 
  
  array[n_trials] vector<lower = 0, upper = 1>[max_all_t + 1] sample_km_est, // sample_target_km_est, sample_non_target_km_est,
                                                              spop_target_km_est, spop_target_obs_cens_km_est; // spop_non_target_km_est, spop_km_est, 

  array[n_cond_group] vector<lower = 0, upper = 1>[max_all_t + 1] cond_sample_km_est, 
    cond_spop_target_km_est, cond_spop_target_obs_cens_km_est; // cond_spop_non_target_km_est, cond_spop_km_est, 

  array[n_trials] vector<lower = 0, upper = 1>[n_pfs_timepoints] forecast_target_pfs_n; //, sample_pfs_n; 
  array[n_cond_group] vector<lower = 0, upper = 1>[n_pfs_timepoints] cond_forecast_target_pfs_n; //, cond_sample_pfs_n;

  array[n_patients] int<lower = 0, upper = 1> forecast_confirmed_response;
  vector<lower = 0, upper = 1>[n_trials] forecast_target_orr;
  vector<lower = 0, upper = 1>[n_cond_group] cond_forecast_target_orr;

  vector<lower = 0>[n_trials] trial_median_pfs = zeros_vector(n_trials);
  vector<lower = 0>[n_cond_group] cond_median_pfs = zeros_vector(n_cond_group);
  
  profile("gen_quant") {
    int right_censored_idx = 1;
    
    for (i in 1:n_patients) {
  
      int visit_start, screening_visit_end, treat_visit_start, visit_end;
      (visit_start, screening_visit_end, treat_visit_start, visit_end) = get_visit_pos(
        patient_visit_pos, i, n_patient_screening_visits[i]);
      int visit_size = get_pos_size(patient_visit_pos, i);
      int n_obs_treat_visits = visit_size - n_patient_screening_visits[i];
  
      int visit_m1_start,  visit_m1_end;
      (visit_m1_start, visit_m1_end) = get_pos(patient_visit_m1_pos, i);
      int visit_m1_size = get_pos_size(patient_visit_m1_pos, i);
  
      int forecast_visit_start, forecast_visit_end;
      (forecast_visit_start, forecast_visit_end) = get_pos(forecast_visits_pos, i);
      int forecast_size = get_pos_size(forecast_visits_pos, i);
      
      array[n_patient_visits[i]] int curr_visits = get_int_sub_array(t_patient_visits, patient_visit_pos, i);
      // Treatment-only visits (exclude screening) for Convention B
      array[visit_size - n_patient_screening_visits[i]] int treat_curr_visits = curr_visits[(n_patient_screening_visits[i] + 1):];
  
      
      // Build forecast time WITH anchor (duplicate last observed week as element 1),
      // so size = n_future + 1. The anchor is needed for state propagation but
      // should be excluded from week mapping logic; hence we pass forecast_time[2:] there.
      array[n_patient_forecast_visits[i] + 1] int forecast_time = linspaced_int_array(
        n_patient_forecast_visits[i] + 1,
        patient_last_obs_visit[i],
        last_predict_visit);
      // Derived slice without anchor for mapping (progression / response) logic
      array[n_patient_forecast_visits[i] > 0 ? n_patient_forecast_visits[i] : 0] int forecast_time_wo_anchor =
        n_patient_forecast_visits[i] > 0 ? forecast_time[2:] : zeros_int_array(0);
      
      // Generate patient states using direct tuple assignment
      // Generate zeros for observed process noise (as it was hardcoded before)
      obs_patient_process_noise[visit_m1_start:visit_m1_end] = rep_matrix(0.0, n_patient_visits[i] - 1, 2);

      mean_patient_log_sld[visit_start:visit_end] = calc_log_sld_mean(states[visit_start:visit_end], sum_tumor_size[visit_start]);
      
      // Generate states, SLD replications, and forecasts
      (forecast_patient_states[forecast_visit_start:forecast_visit_end], 
       rep_patient_log_sld[visit_start:visit_end], 
       forecast_patient_log_sld[forecast_visit_start:forecast_visit_end]) = 
        generate_patient_states_rng(
          states[visit_start:visit_end],
          forecast_time,
          patient_log_decrease_rate[i], patient_log_growth_rate[i],
          sum_tumor_size[visit_start], 
          0.0001, 0.0001, // exp(patient_log_growth_lag[i]), exp(pop_log_growth_transition_rate),
          rep_matrix(0.0, forecast_size, 2), // Hardcode zeros for forecast process noise
          measure_sd
        );

      forecast_mean_patient_log_sld[forecast_visit_start:forecast_visit_end] = 
        calc_log_sld_mean(forecast_patient_states[forecast_visit_start:forecast_visit_end], sum_tumor_size[visit_start]);
      
      // Generate forecast process noise conditionally (as it was in the original code)
      if (independ_long_process_noise) {
        forecast_patient_process_noise[forecast_visit_start:forecast_visit_end] = multi_normal_rng(
          forecast_size, pop_process_sd, independ_cross_process_noise ? diag_matrix(ones_vector(2)) : L_process_corr
        );
      } else {
        forecast_patient_process_noise[forecast_visit_start:forecast_visit_end] = multi_normal_rng(
          obs_patient_process_noise[visit_m1_start:visit_m1_end],
          get_int_sub_array(t_patient_visits, patient_visit_pos, i)[2:],
          forecast_time[2:],
          exp(log_pop_tumor_gp_rho),
          pop_process_sd,
          independ_cross_process_noise ? diag_matrix(ones_vector(2)) : L_process_corr,
          delta
        );
      }      
      
    //   (spop_non_target_pfs[i], spop_non_target_right_censored[i]) = survival_time_rng(log_cond_prob_surv[1, i]);
     
      array[n_obs_treat_visits + forecast_size] int full_predict_recist = calculate_target_recist(
        // exp(append_row(rep_patient_log_sld[visit_start:visit_end], 
        //               forecast_patient_log_sld[forecast_visit_start:forecast_visit_end])) * 10,
        // This is a temporary fix to address the problem of very high measurement noise because of imperfect mechanistic model fit.
        exp(append_row(mean_patient_log_sld[visit_start:visit_end], 
                      forecast_mean_patient_log_sld[forecast_visit_start:forecast_visit_end])) * 10,
        n_patient_screening_visits[i]
      );
      
      rep_recist[treat_visit_start:visit_end] = full_predict_recist[:n_obs_treat_visits]; 
      forecast_recist[forecast_visit_start:forecast_visit_end] = full_predict_recist[(n_obs_treat_visits + 1):];
     
      int forecast_confirmed_response_week = 0, confirmed_response_censored = 0; 

      array[visit_size + forecast_size] int curr_recist;
      curr_recist[:visit_size] = recist[visit_start:visit_end]; // includes screening, will slice below when calling find_first_week

      sample_target_pfs[i] = pfs[i] + interval_censored[i] + 1; 
      sample_target_right_censored[i] = right_censored[i];
          
      if (n_patient_forecast_visits[i] > 0) {
        curr_recist[(visit_size + 1):] = forecast_recist[forecast_visit_start:forecast_visit_end];

        (forecast_confirmed_response_week, confirmed_response_censored) = 
          find_first_week(curr_recist[(n_patient_screening_visits[i] + 1):], { PR, CR }, 2, 
                          treat_curr_visits, forecast_time_wo_anchor, max_all_t);
         
        if (right_censored[i]) {
          (sample_target_pfs[i], sample_target_right_censored[i]) = 
            find_first_week(curr_recist[(n_patient_screening_visits[i] + 1):], { PD }, 1, treat_curr_visits, forecast_time_wo_anchor, max_all_t);

        //   (forecast_non_target_pfs[right_censored_idx], forecast_non_target_right_censored[right_censored_idx]) = survival_time_rng(
        //     log_cond_prob_surv[1, i], pfs[i] + interval_censored[i], right_censored[i], 0
        //   );

          (forecast_target_pfs[right_censored_idx], forecast_target_right_censored[right_censored_idx]) = 
            find_first_forecast_week(forecast_recist[forecast_visit_start:forecast_visit_end], { PD }, 1, forecast_time_wo_anchor, max_all_t);
          
        //   forecast_pfs[right_censored_idx] = min(forecast_target_pfs[right_censored_idx], forecast_non_target_pfs[right_censored_idx]);
        //   forecast_right_censored[right_censored_idx] = forecast_target_right_censored[right_censored_idx] && forecast_non_target_right_censored[right_censored_idx];
          
          right_censored_idx += 1;
        } 
        
        (spop_target_pfs[i], spop_target_right_censored[i]) = find_first_week(
          append_array(rep_recist[treat_visit_start:visit_end], 
                       forecast_recist[forecast_visit_start:forecast_visit_end]), 
          { PD }, 1, treat_curr_visits, forecast_time_wo_anchor, max_all_t); 
      } else {
        (spop_target_pfs[i], spop_target_right_censored[i]) = find_first_week(
          rep_recist[treat_visit_start:visit_end], { PD }, 1, treat_curr_visits, forecast_time_wo_anchor, max_all_t); 

        (forecast_confirmed_response_week, confirmed_response_censored) = find_first_week(
          curr_recist[(n_patient_screening_visits[i] + 1):], { PR, CR }, 2, treat_curr_visits, forecast_time_wo_anchor, max_all_t);
      }
     
      spop_target_obs_cens_right_censored[i] = right_censored[i]; 
      spop_target_obs_cens_pfs[i] = right_censored[i] ? min(spop_target_pfs[i], pfs[i]) : spop_target_pfs[i];
      
    //   spop_pfs[i] = min(spop_non_target_pfs[i] + 1, 
    //                             max(0, spop_target_pfs[i])); // BUG a couple of patients end up with negative weeks. We need to figure out why.
    //   spop_right_censored[i] = spop_target_right_censored[i] && spop_non_target_right_censored[i]; 

      forecast_confirmed_response[i] = 
        !confirmed_response_censored && forecast_confirmed_response_week < sample_target_pfs[i];
    }
      
    for (s in 1:n_trials) {
      if (get_pos_size(trial_patient_pos, s) > 0) {
        forecast_target_orr[s] = mean(get_int_sub_array(forecast_confirmed_response, trial_patient_pos, s));

        int n_curr_uncensored_obs = get_pos_size(trial_right_uncensored_pos, s);
        array[n_curr_uncensored_obs] int curr_uncensored_obs = get_int_sub_array(right_uncensored_patients, trial_right_uncensored_pos, s);
        int n_curr_right_censored = get_pos_size(trial_right_censored_pos, s);

        array[n_curr_uncensored_obs + n_curr_right_censored] int 
            curr_sample_pfs = get_int_sub_array(sample_target_pfs, trial_patient_pos, s),
            curr_sample_right_censored = get_int_sub_array(sample_target_right_censored, trial_patient_pos, s);

        sample_km_est[s] = estimate_kaplan_meier(curr_sample_pfs, curr_sample_right_censored, max_all_t).1; 

        array[n_curr_uncensored_obs + n_curr_right_censored] int 
            curr_spop_target_pfs = get_int_sub_array(spop_target_pfs, trial_patient_pos, s),
            curr_spop_target_right_censored = get_int_sub_array(spop_target_right_censored, trial_patient_pos, s); 
        
        spop_target_km_est[s] = estimate_kaplan_meier(curr_spop_target_pfs, curr_spop_target_right_censored, max_all_t, 0).1; 
                                               
        spop_target_obs_cens_km_est[s] = estimate_kaplan_meier(get_int_sub_array(spop_target_obs_cens_pfs, trial_patient_pos, s), 
                                               get_int_sub_array(spop_target_obs_cens_right_censored, trial_patient_pos, s), 
                                               max_all_t, 0).1; 
                                               
        // spop_non_target_km_est[s] = estimate_kaplan_meier(get_int_sub_array(spop_non_target_pfs, trial_patient_pos, s), 
        //                                        get_int_sub_array(spop_non_target_right_censored, trial_patient_pos, s), 
        //                                        max_all_t, 0).1; 
        
        // spop_km_est[s] = estimate_kaplan_meier(get_int_sub_array(spop_pfs, trial_patient_pos, s), 
        //                                        get_int_sub_array(spop_right_censored, trial_patient_pos, s), 
        //                                        max_all_t, 0).1;
                                               
        trial_median_pfs[s] = km_median(sample_km_est[s]).1;

        for (n in 1:n_pfs_timepoints) {
          forecast_target_pfs_n[s, n] = calc_km_pfs_n(sample_km_est[s], months_to_weeks(pfs_timepoints[n])); 
        }
      } else {
        spop_target_km_est[s] = zeros_vector(max_all_t + 1);
        spop_target_obs_cens_km_est[s] = zeros_vector(max_all_t + 1);
        // spop_non_target_km_est[s] = zeros_vector(max_all_t + 1);
        // spop_km_est[s] = zeros_vector(max_all_t + 1);
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
                                             
    //   cond_spop_non_target_km_est[c] = estimate_kaplan_meier(spop_non_target_pfs[curr_group_patients],
    //                                          spop_non_target_right_censored[curr_group_patients],
    //                                          max_all_t, 0).1; 
      
    //   cond_spop_km_est[c] = estimate_kaplan_meier(spop_pfs[curr_group_patients],
    //                                          spop_right_censored[curr_group_patients],
    //                                          max_all_t, 0).1; 

      cond_median_pfs[c] = km_median(cond_sample_km_est[c]).1;

      for (n in 1:n_pfs_timepoints) {
        cond_forecast_target_pfs_n[c, n] = calc_km_pfs_n(cond_sample_km_est[c], months_to_weeks(pfs_timepoints[n])); 
      }
    }
  }