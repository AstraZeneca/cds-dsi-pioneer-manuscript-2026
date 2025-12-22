functions {
  #include "../util.stan"
  #include "../pos.stan"
  #include "../gp.stan"
  #include "../pfs_functions.stan"
  #include "../lfo.stan"
  #include "_sf_functions.stan"
  #include "../recist.stanfunctions"
}  

data {
  #include "../base_data.stan"
  #include "../tumor/base_data.stan"
  #include "_sf_outcomes_info.stan"

  #include "modules/measurement/hyperparams.stan"
  #include "modules/other_events/data.stan"
  #include "modules/tr/hyperparams.stan"
  #include "modules/frac/hyperparams.stan"
  #include "modules/init/hyperparams.stan"
  #include "modules/other_events/hyperparams.stan"
  #include "modules/measurement/flags.stan"
  #include "modules/other_events/flags.stan"
  #include "modules/tr/flags.stan"
  #include "modules/frac/flags.stan"
  #include "modules/init/flags.stan"

  #include "_sf-ssls-lfo-data.stan"
  
  // --- LFO CV specific ---
  int<lower = 0, upper = 1> train_beyond_cutoff;
} 

transformed data {
  print("cutoff_calendar_day = ", cutoff_calendar_day);
  
  #include "../base_transformed_data.stan"
  #include "../tumor/tumor_transformed_data.stan"
  #include "modules/measurement/transformed_data.stan"
  #include "_sf_transformed_data.stan"
  #include "modules/other_events/transformed_data.stan"
  #include "_lfo_transformed_data.stan"
}

parameters {
  #include "modules/measurement/parameters.stan"
  #include "modules/other_events/parameters.stan"
  #include "modules/tr/parameters.stan"
  #include "modules/frac/parameters.stan"
  #include "modules/init/parameters.stan"
}

transformed parameters {
  #include "modules/measurement/transformed_parameters.stan"
  #include "modules/tr/transformed_parameters.stan"
  #include "modules/frac/transformed_parameters.stan"
  #include "modules/init/transformed_parameters.stan"
  #include "_sf_transformed_parameters.stan"
  #include "modules/other_events/transformed_parameters.stan"
}

model {
  #include "modules/measurement/priors.stan"
  #include "modules/other_events/priors.stan"
  #include "modules/tr/priors.stan"
  #include "modules/frac/priors.stan"
  #include "modules/init/priors.stan"

  if (fit_tumor_data) {
    // --- LFO CV specific ---
    for (i in 1:n_patients) {
      if (cutoff_last_visit_idx[i] > 0) {
        int visit_start, visit_end;
        (visit_start, visit_end) = get_pos(patient_visit_pos, i);

        int cutoff_idx = cutoff_last_visit_idx[i];

        normalized_sld[visit_start:cutoff_idx] ~ sf_log_space_obs(states[visit_start:cutoff_idx], measure_sd, log_lod - log_baseline_sld[i]);
      }
    }

    // Other events likelihood contribution (cutoff-aware)
    if (n_causes > 0) {
      matrix[n_cutoff_observed_patients, n_causes] patient_response_lp = rep_matrix(0, n_cutoff_observed_patients, n_causes);

      // Use cutoff-censored data and enforce time window to prevent data leakage
      patient_response_lp[, 1] = calc_pch_loglik(
        cutoff_ic_other_events_pfs, 
        cutoff_other_events_right_censored, 
        zeros_int_array(n_cutoff_observed_patients), // interval_censored
        0, 
        log_cond_prob_surv[1, cutoff_observed_patients],
        ones_int_array(n_cutoff_observed_patients), // start_from
        train_beyond_cutoff ? rep_array(max_all_t, n_cutoff_observed_patients) : cutoff_last_visit_week // end_at
      );

      target += sum(patient_response_lp);
    }
  }
}

generated quantities {
  // Include comprehensive endpoints that integrate other events with target RECIST
  #include "_lfo_endpoints_generated_quantities.stan"
  
  // Sentinel PD+1 not allowed by bound; assert below ensures no leakage
  array[sum(n_patient_testing_visits)] int<lower = CR, upper = PD> oos_recist = rep_array(PD + 1, sum(n_patient_testing_visits));   

  for (i in last_visit_calendar_day_sort_idx[testing_patient_idx[1]:]) {
    int visit_start, visit_screening_end, visit_treat_pos, visit_end;
    (visit_start, visit_screening_end, visit_treat_pos, visit_end) = get_visit_pos(patient_visit_pos, i, n_patient_screening_visits[i]);

    int cutoff_idx = cutoff_last_visit_idx[i];
    cutoff_idx = cutoff_idx > 0 ? cutoff_idx : visit_start;

    int visit_size = cutoff_idx - visit_start + 1;
    int treat_visit_size = max(0, cutoff_idx - visit_treat_pos + 1);
    int start_idx = testing_start_idx[1, i];
    // Only generate OOS predictions for patients who:
    // 1) Have post-cutoff visits at the first cutoff (start_idx > 0), AND
    // 2) Were observed before/at the cutoff (cutoff_observed_mask[i] == 1)
    // This excludes newly enrolled patients who entered the study after the cutoff.
    if (start_idx > 0 && cutoff_observed_mask[i] == 1) {
      int n_oos_visits = visit_end - start_idx + 1; 
    
      array[n_oos_visits + 1] int forecast_time = get_int_sub_array(t_patient_visits, patient_visit_pos, i)[visit_size:];      
      
      // Extract observed visit indices for this patient
      array[visit_size] int visit_indices = t_patient_visit_idx[visit_start:cutoff_idx];

      // Extract observed states and compute forecast states
      matrix[visit_size, 2] patient_states;
      matrix[n_oos_visits, 2] forecast_patient_states;

      if (enable_patient_process_noise_tr) {
        // Process noise ON: Extract from dense grid computed in transformed_parameters
        patient_states[, 1] = to_vector(states_full_grid[1][i, visit_indices]);
        patient_states[, 2] = to_vector(states_full_grid[2][i, visit_indices]);

        if (n_oos_visits > 0) {
          forecast_patient_states[, 1] = to_vector(states_full_grid[1][i, forecast_time[2:]]);
          forecast_patient_states[, 2] = to_vector(states_full_grid[2][i, forecast_time[2:]]);
        }
      } else {
        // Process noise OFF: Use states directly and compute forecast on-the-fly
        patient_states = states[visit_start:cutoff_idx];

        if (n_oos_visits > 0) {
          // Compute forecast states using constant rates
          matrix[n_oos_visits + 1, 2] full_forecast_expected;  // Unused but required by tuple return
          matrix[n_oos_visits + 1, 2] full_forecast;
          (full_forecast_expected, full_forecast) = sf_log_space_trajectory_ncp(
            patient_states[visit_size],  // Last observed state as initial
            forecast_time,
            exp(patient_log_decrease_rate[i, 1]),
            exp(patient_log_growth_rate[i, 1]),
            negative_infinity(),  // growth lag disabled
            1.0,  // growth transition
            rep_matrix(0.0, n_oos_visits, 2),  // No process noise
            0  // No debug
          );
          forecast_patient_states = full_forecast[2:];  // Skip anchor
        }
      }
      
      // Calculate replicated SLD for observed visits
      vector[visit_size] rep_patient_log_sld = zeros_vector(visit_size);
      rep_patient_log_sld[1] = log(sum_tumor_size[visit_start]);
      if (visit_size > 1) {
        rep_patient_log_sld[2:] = to_vector(normal_rng(
          calc_log_sld_mean(patient_states[2:], sum_tumor_size[visit_start]),
          rep_vector(measure_sd, visit_size - 1)
        ));
      }
      
      // Calculate mean log SLD (deterministic, no measurement noise)
      vector[visit_size] rep_mean_patient_log_sld = 
        calc_log_sld_mean(patient_states, sum_tumor_size[visit_start]);
      
      // Calculate forecast SLD with measurement noise
      vector[n_oos_visits] forecast_patient_log_sld = zeros_vector(n_oos_visits);
      vector[n_oos_visits] forecast_mean_patient_log_sld = zeros_vector(n_oos_visits);
      if (n_oos_visits > 0) {
        forecast_mean_patient_log_sld = 
          calc_log_sld_mean(forecast_patient_states, sum_tumor_size[visit_start]);
        forecast_patient_log_sld = to_vector(normal_rng(
          forecast_mean_patient_log_sld,
          rep_vector(measure_sd, n_oos_visits)
        ));
      }

      // calculate_target_recist returns RECIST for treatment visits only (screening dropped),
      // so length(full_predict_recist) == treat_visit_size + n_oos_visits.
      array[treat_visit_size + n_oos_visits] int full_predict_overall_recist = calculate_target_recist(
        exp(append_row(rep_mean_patient_log_sld, forecast_mean_patient_log_sld)) * 10,
        n_patient_screening_visits[i]
      );
      
      // Check if there was already a PD in the in-sample (observed) period
      // If so, all forecast visits must also be PD (overall RECIST remains PD once reached)
      int had_insample_pd = treat_visit_size > 0 && full_predict_overall_recist[treat_visit_size] == PD;
      
      // Mark all forecast visits as PD if:
      // 1) Patient had PD in the in-sample period (had_insample_pd == 1), OR
      // 2) Other events cause PD in the forecast period
      if (had_insample_pd == 1) {
        // All forecast visits are PD since patient already had PD before cutoff
        full_predict_overall_recist[(treat_visit_size + 1):] = rep_array(PD, n_oos_visits);
      } else {
        // Integrate other events PFS to mark RECIST as PD when other events cause progression
        // Since we only process cutoff-observed patients (cutoff_observed_mask[i] == 1),
        // we can always use the already-calculated sample_other_events_pfs from _lfo_endpoints_generated_quantities.stan
        int cutoff_patient_idx = patient_to_cutoff_idx[i];
        int forecast_other_events_pfs = sample_other_events_pfs[cutoff_patient_idx];
        int forecast_other_events_censored = sample_other_events_right_censored[cutoff_patient_idx];

        if (!forecast_other_events_censored) {
          // Other events PD occurs at week forecast_other_events_pfs
          // Find first forecast visit at or after other events PFS
          int forecast_other_events_visit_idx = 1;
          while (forecast_other_events_visit_idx <= n_oos_visits && forecast_time[forecast_other_events_visit_idx + 1] < forecast_other_events_pfs) {
            forecast_other_events_visit_idx += 1;
          }
          
          // Mark all subsequent forecast visits as PD (from the first visit >= other events PFS onward)
          if (forecast_other_events_visit_idx <= n_oos_visits) {
            full_predict_overall_recist[(treat_visit_size + forecast_other_events_visit_idx):] = rep_array(PD, n_oos_visits - forecast_other_events_visit_idx + 1);
          }
        }
      }

      int oos_recist_start, oos_recist_end;
      (oos_recist_start, oos_recist_end) = get_pos(testing_visit_pos, i);

      // We write only the out-of-sample part: drop the in-sample treatment visits (treat_visit_size)
      // and keep exactly n_oos_visits RECIST values.
      if (n_patient_testing_visits[i] > 0) {
        // Bounds/sanity checks for the per-patient OOS slice
        assert_equal(oos_recist_end - oos_recist_start + 1, n_oos_visits);
        oos_recist[oos_recist_start:oos_recist_end] = full_predict_overall_recist[(treat_visit_size + 1):];
      }
    }
  }

  // Separate tracking for interpretability
  array[n_cutoffs, n_cutoffs] vector[n_all_testing_patients] patient_log_lik_tumor;   // P(SLD | tumor model)
  array[n_cutoffs, n_cutoffs] vector[n_all_testing_patients] patient_log_lik_oe;      // P(OE PFS | OE model)  
  array[n_cutoffs, n_cutoffs] vector[n_all_testing_patients] patient_log_lik;         // P(SLD, OE PFS | joint model)

  array[n_cutoffs, n_cutoffs] matrix<lower = 0>[PD, PD] oos_recist_confusion_matrix; // rows = observed, cols = predicted

  // Index usage notes:
  //   start_idx = testing_start_idx[n, i] is first post-cutoff-n visit (0 if none yet)
  //   end_idx   = testing_end_idx[n, m+1, i] (inclusive) for horizon ending at cutoff m (< n_cutoffs), else patient's last visit
  //   first_start_idx = testing_start_idx[1, i] anchor for contiguous per-patient OOS slice
  //   test_start_offset = start_idx - first_start_idx (>=0) positions current window inside that slice
  for (n in 1:n_cutoffs) {
    int n_curr_patients = n_patients - testing_patient_idx[n] + 1; // How many patients after the current patient index
    int curr_first_testing_patient_idx = n_all_testing_patients - n_curr_patients + 1; 
    array[n_curr_patients] int curr_patients = last_visit_calendar_day_sort_idx[testing_patient_idx[n]:]; // Who are these patients
    
    for (m in 1:n_cutoffs) {
      patient_log_lik_tumor[n, m] = zeros_vector(n_all_testing_patients);
      patient_log_lik_oe[n, m] = zeros_vector(n_all_testing_patients);
      patient_log_lik[n, m] = zeros_vector(n_all_testing_patients);
      oos_recist_confusion_matrix[n, m] = rep_matrix(0, PD, PD);
    
      if (m >= n) {
        for (i_idx in 1:n_curr_patients) {
          // Note: i is the original patient ID (1-based index from input data), not a sort position.
          // curr_patients contains original patient IDs that were reordered by sorting on last_visit_calendar_day
          int i = curr_patients[i_idx];

          int visit_start, visit_end, first_start_idx = testing_start_idx[1, i];
          (visit_start, visit_end) = get_pos(patient_visit_pos, i); 

          int start_idx = testing_start_idx[n, i];
          int end_idx = m < n_cutoffs ? testing_end_idx[n, m + 1, i] : visit_end;

          // Only evaluate patients who were observed at cutoff (exclude newly enrolled patients)
          // cutoff_observed_mask[i] == 1 means patient had at least one visit before/at cutoff
          if (start_idx > 0 && end_idx >= start_idx && cutoff_observed_mask[i] == 1) {
            int patient_idx = curr_first_testing_patient_idx + i_idx - 1;
            
            // Component 1: Tumor model log-likelihood using observed SLD
            real tumor_ll = sf_log_space_obs_lpdf(
                normalized_sld[start_idx:end_idx] | states[start_idx:end_idx], 
                measure_sd, log_lod - log_baseline_sld[i]);
            
            patient_log_lik_tumor[n, m, patient_idx] = tumor_ll;
            
            // Component 2: Other events model log-likelihood using OBSERVED other events PFS
            real oe_ll = 0;
            int cutoff_patient_idx = patient_to_cutoff_idx[i];

            if (n_causes > 0 && cutoff_patient_idx > 0) {
              // Get test window boundaries in weeks
              int test_start_week = t_patient_visits[start_idx];
              int test_end_week = t_patient_visits[end_idx];
              
              // Use OBSERVED other events PFS from cutoff-censored data
              array[1] int obs_oe_pfs = {cutoff_ic_other_events_pfs[cutoff_patient_idx]};
              array[1] int obs_oe_censored = {cutoff_other_events_right_censored[cutoff_patient_idx]};
              array[1] int obs_oe_ic = {0}; // not interval censored
              array[1] int test_start = {test_start_week};
              array[1] int test_end = {test_end_week};
              
              // Extract single patient's survival probabilities as a 1-row matrix
              // log_cond_prob_surv is array[n_causes] matrix[n_patients, max_all_t]
              // We need matrix[1, max_all_t] for this single patient
              int patient_row = cutoff_observed_patients[cutoff_patient_idx];
              matrix[1, max_all_t] patient_log_surv = log_cond_prob_surv[1, patient_row:patient_row];
              
              // Calculate log-likelihood using the same function as in model block
              oe_ll = calc_pch_loglik(
                obs_oe_pfs,
                obs_oe_censored,
                obs_oe_ic,
                0, // ignore_interval_censoring
                patient_log_surv,
                test_start,
                test_end
              )[1]; // Extract single element from returned vector
            }
            
            patient_log_lik_oe[n, m, patient_idx] = oe_ll;
            
            // Joint log-likelihood: log P(SLD, OE PFS | θ) = log P(SLD | θ) + log P(OE PFS | θ)
            patient_log_lik[n, m, patient_idx] = tumor_ll + oe_ll;

            // Below part is for OOS RECIST confusion matrix calculation

            int oos_recist_start, oos_recist_end;
            (oos_recist_start, oos_recist_end) = get_pos(testing_visit_pos, i);

            int n_curr_testing_visits = end_idx - start_idx + 1;
            int test_start_offset = start_idx - first_start_idx;
            assert_greater_than_or_equal(test_start_offset, 0);

            // The patient slice in oos_recist spans all OOS visits from the first cutoff:
            assert_equal(oos_recist_start + n_patient_testing_visits[i] - 1, oos_recist_end);
            // For current (n, m) window, ensure we don't step past the end of that slice:
            assert_greater_than_or_equal(oos_recist_end, oos_recist_start + test_start_offset + n_curr_testing_visits - 1);

            for (t_idx in 1:n_curr_testing_visits) {
              // Predicted RECIST must be within 1..PD
              assert_less_or_equal(oos_recist[oos_recist_start + test_start_offset + t_idx - 1], PD);

              int obs_val = recist[start_idx + t_idx - 1];
              int pred_val = oos_recist[oos_recist_start + test_start_offset + t_idx - 1];

              if (obs_val <= PD) { // Ignore first pre-screening visits that don't have a response yet 
                oos_recist_confusion_matrix[n, m][obs_val, pred_val] += 1;
              }
            }
          }
        }
      } else {
        patient_log_lik_tumor[n, m] = rep_vector(negative_infinity(), n_all_testing_patients); 
        patient_log_lik_oe[n, m] = rep_vector(negative_infinity(), n_all_testing_patients);
        patient_log_lik[n, m] = rep_vector(negative_infinity(), n_all_testing_patients);
      }
    }
  } 
}
