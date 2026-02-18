functions {
  #include "util.stanfunctions"
  #include "pos.stanfunctions"
  #include "gp.stanfunctions"
  #include "pfs.stanfunctions"
  #include "lfo.stanfunctions"
  #include "multistate.stanfunctions"
  #include "modules/state_space/sf.stanfunctions"
  #include "modules/tumor/tumor.stanfunctions"
}

data {
  #include "_base_data.stan"
  #include "modules/tumor/data.stan"
  #include "modules/tumor/hyperparams.stan"
  #include "modules/state_space/data.stan"
  #include "modules/multistate/flags.stan"
  #include "modules/multistate/data.stan"
  #include "modules/multistate/hyperparams.stan"
  #include "modules/tr/hyperparams.stan"
  #include "modules/frac/hyperparams.stan"
  #include "modules/init/hyperparams.stan"
  #include "modules/tr/flags.stan"
  #include "modules/frac/flags.stan"
  #include "modules/init/flags.stan"

  int<lower = 0, upper = 1> fit_multistate_data;

  #include "modules/state_space/lfo_data.stan"
}

transformed data {
  print("cutoff_calendar_day = ", cutoff_calendar_day);

  #include "_base_transformed_data.stan"
  #include "modules/tumor/transformed_data.stan"
  #include "modules/tr/transformed_data.stan"
  #include "modules/frac/transformed_data.stan"
  #include "modules/init/transformed_data.stan"
  #include "modules/state_space/transformed_data.stan"
  #include "modules/multistate/transformed_data.stan"
  #include "_lfo_transformed_data.stan"
}

parameters {
  #include "modules/tumor/parameters.stan"
  #include "modules/multistate/parameters.stan"
  #include "modules/tr/parameters.stan"
  #include "modules/frac/parameters.stan"
  #include "modules/init/parameters.stan"
}

transformed parameters {
  #include "modules/tr/transformed_parameters.stan"
  #include "modules/frac/transformed_parameters.stan"
  #include "modules/init/transformed_parameters.stan"
  #include "modules/state_space/transformed_parameters.stan"
  #include "_ms_time_varying_covar.stan"
  #include "modules/multistate/transformed_parameters.stan"
}

model {
  #include "modules/tumor/priors.stan"
  #include "modules/multistate/priors.stan"
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

        normalized_sld[visit_start:cutoff_idx] ~ sf_log_space_obs(states[visit_start:cutoff_idx], measure_sd_sld, log_lod - log_baseline_sld[i]);
      }
    }

    // Multistate likelihood contribution (cutoff-aware)
    if (enable_ms_01) {
      target += sum(calc_ms_single_transition_loglik(
        cutoff_ms_time_01,
        cutoff_ms_censored_01,
        ms_log_cond_surv_01[cutoff_observed_patients]
      ));
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
          rep_vector(measure_sd_sld, visit_size - 1)
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
          rep_vector(measure_sd_sld, n_oos_visits)
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
        // Integrate multistate PFS to mark RECIST as PD when multistate events cause progression
        // Since we only process cutoff-observed patients (cutoff_observed_mask[i] == 1),
        // we can always use the already-calculated sample_ms_pfs from _lfo_endpoints_generated_quantities.stan
        int cutoff_patient_idx = patient_to_cutoff_idx[i];
        int forecast_ms_pfs = sample_ms_pfs[cutoff_patient_idx];
        int forecast_ms_censored = sample_ms_right_censored[cutoff_patient_idx];

        if (!forecast_ms_censored) {
          // Multistate PD occurs at week forecast_ms_pfs
          // Find first forecast visit at or after multistate PFS
          int forecast_ms_visit_idx = 1;
          while (forecast_ms_visit_idx <= n_oos_visits && forecast_time[forecast_ms_visit_idx + 1] < forecast_ms_pfs) {
            forecast_ms_visit_idx += 1;
          }

          // Mark all subsequent forecast visits as PD (from the first visit >= multistate PFS onward)
          if (forecast_ms_visit_idx <= n_oos_visits) {
            full_predict_overall_recist[(treat_visit_size + forecast_ms_visit_idx):] = rep_array(PD, n_oos_visits - forecast_ms_visit_idx + 1);
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
  // Array dimensions are configurable:
  //   - Exact LFO: max_n_rows=1, max_forecast_horizon=2 → [1, 2] arrays (minimal memory)
  //   - PSIS LFO: max_n_rows=n_cutoffs, max_forecast_horizon=n_cutoffs → [n, n] arrays (for approximation)
  array[max_n_rows, max_forecast_horizon] vector[n_all_testing_patients] patient_log_lik_tumor;   // P(SLD | tumor model)
  array[max_n_rows, max_forecast_horizon] vector[n_all_testing_patients] patient_log_lik_oe;      // P(OE PFS | OE model)
  array[max_n_rows, max_forecast_horizon] vector[n_all_testing_patients] patient_log_lik;         // P(SLD, OE PFS | joint model)

  array[max_n_rows, max_forecast_horizon] matrix<lower = 0>[PD, PD] oos_recist_confusion_matrix; // rows = observed, cols = predicted

  // Index usage notes:
  //   start_idx = testing_start_idx[n, i] is first post-cutoff-n visit (0 if none yet)
  //   end_idx   = testing_end_idx[n, m_abs+1, i] (inclusive) for horizon ending at cutoff m_abs (< n_cutoffs), else patient's last visit
  //   first_start_idx = testing_start_idx[1, i] anchor for contiguous per-patient OOS slice
  //   test_start_offset = start_idx - first_start_idx (>=0) positions current window inside that slice
  //
  // Loop bounds explanation:
  //   - max_n_rows controls how many training cutoffs to compute (1 for exact LFO, n_cutoffs for PSIS)
  //   - max_forecast_horizon controls the forecast window size (2 for exact LFO, n_cutoffs for PSIS)
  //   - m_rel is the relative column index for array storage (1, 2, ...)
  //   - m_abs is the absolute cutoff index for data access (n, n+1, ...)
  for (n in 1:max_n_rows) {
    int n_curr_patients = n_patients - testing_patient_idx[n] + 1; // How many patients after the current patient index
    int curr_first_testing_patient_idx = n_all_testing_patients - n_curr_patients + 1;
    array[n_curr_patients] int curr_patients = last_visit_calendar_day_sort_idx[testing_patient_idx[n]:]; // Who are these patients

    // Only compute for forecast window: m_abs in [n, min(n + max_forecast_horizon - 1, n_cutoffs)]
    int m_end = min(n + max_forecast_horizon - 1, n_cutoffs);
    for (m_abs in n:m_end) {
      // m_rel is the relative column index for array storage (1-based: 1, 2, ...)
      int m_rel = m_abs - n + 1;

      patient_log_lik_tumor[n, m_rel] = zeros_vector(n_all_testing_patients);
      patient_log_lik_oe[n, m_rel] = zeros_vector(n_all_testing_patients);
      patient_log_lik[n, m_rel] = zeros_vector(n_all_testing_patients);
      oos_recist_confusion_matrix[n, m_rel] = rep_matrix(0, PD, PD);

      for (i_idx in 1:n_curr_patients) {
        // Note: i is the original patient ID (1-based index from input data), not a sort position.
        // curr_patients contains original patient IDs that were reordered by sorting on last_visit_calendar_day
        int i = curr_patients[i_idx];

        int visit_start, visit_end, first_start_idx = testing_start_idx[1, i];
        (visit_start, visit_end) = get_pos(patient_visit_pos, i);

        int start_idx = testing_start_idx[n, i];
        int end_idx = m_abs < n_cutoffs ? testing_end_idx[n, m_abs + 1, i] : visit_end;

        // Only evaluate patients who were observed at cutoff (exclude newly enrolled patients)
        // cutoff_observed_mask[i] == 1 means patient had at least one visit before/at cutoff
        if (start_idx > 0 && end_idx >= start_idx && cutoff_observed_mask[i] == 1) {
          int patient_idx = curr_first_testing_patient_idx + i_idx - 1;

          // Component 1: Tumor model log-likelihood using observed SLD
          real tumor_ll = sf_log_space_obs_lpdf(
              normalized_sld[start_idx:end_idx] | states[start_idx:end_idx],
              measure_sd_sld, log_lod - log_baseline_sld[i]);

          patient_log_lik_tumor[n, m_rel, patient_idx] = tumor_ll;

          // Component 2: Multistate model log-likelihood using ACTUAL observed PFS
          // (not cutoff-censored - we want true OOS evaluation against real outcomes)
          real ms_ll = 0;

          if (enable_ms_01) {
            // Get test window boundaries in weeks
            int test_start_week = t_patient_visits[start_idx];
            int test_end_week = t_patient_visits[end_idx];

            // Use ACTUAL observed multistate PFS (not cutoff-censored)
            // This ensures proper OOS evaluation against real outcomes
            array[1] int obs_ms_time = {ms_time_01[i]};
            array[1] int obs_ms_censored = {ms_censored_01[i]};
            array[1] int test_start = {test_start_week};
            array[1] int test_end = {test_end_week};

            // Extract single patient's survival probabilities as a 1-row matrix
            matrix[1, max_all_t] patient_log_surv = ms_log_cond_surv_01[i:i];

            // Calculate log-likelihood using the same function as in model block
            ms_ll = calc_pch_loglik(
              obs_ms_time,
              obs_ms_censored,
              zeros_int_array(1), // no interval censoring
              0, // ignore_interval_censoring
              patient_log_surv,
              test_start,
              test_end
            )[1]; // Extract single element from returned vector
          }

          patient_log_lik_oe[n, m_rel, patient_idx] = ms_ll;

          // Joint log-likelihood: log P(SLD, OE PFS | θ) = log P(SLD | θ) + log P(OE PFS | θ)
          patient_log_lik[n, m_rel, patient_idx] = tumor_ll + ms_ll;

          // Below part is for OOS RECIST confusion matrix calculation

          int oos_recist_start, oos_recist_end;
          (oos_recist_start, oos_recist_end) = get_pos(testing_visit_pos, i);

          int n_curr_testing_visits = end_idx - start_idx + 1;
          int test_start_offset = start_idx - first_start_idx;
          assert_greater_than_or_equal(test_start_offset, 0);

          // The patient slice in oos_recist spans all OOS visits from the first cutoff:
          assert_equal(oos_recist_start + n_patient_testing_visits[i] - 1, oos_recist_end);
          // For current (n, m_abs) window, ensure we don't step past the end of that slice:
          assert_greater_than_or_equal(oos_recist_end, oos_recist_start + test_start_offset + n_curr_testing_visits - 1);

          for (t_idx in 1:n_curr_testing_visits) {
            // Predicted RECIST must be within 1..PD
            assert_less_or_equal(oos_recist[oos_recist_start + test_start_offset + t_idx - 1], PD);

            int obs_val = recist[start_idx + t_idx - 1];
            int pred_val = oos_recist[oos_recist_start + test_start_offset + t_idx - 1];

            if (obs_val <= PD) { // Ignore first pre-screening visits that don't have a response yet
              oos_recist_confusion_matrix[n, m_rel][obs_val, pred_val] += 1;
            }
          }
        }
      }
    }
  }
}
