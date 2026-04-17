// ============================================================================
// MS-Standalone LFO: Held-out log-likelihoods per transition
// ============================================================================
// Computes patient_log_lik[n, m, i] for compatibility with lfo_log_lik_rvar().
//
// For each cutoff window (n = training cutoff, m = evaluation horizon):
//   - Identify patients with events in [cutoff_n, cutoff_m]
//   - Compute conditional log-likelihood: P(event | survived to cutoff, theta)
//   - All three transitions (0→1, 0→2, 1→2) contribute
//
// Requires in scope:
//   n_cutoffs, max_n_rows, max_forecast_horizon, testing_patient_idx[]
//   last_visit_calendar_day_sort_idx[], cutoff_last_visit_week[]
//   n_all_testing_patients, all_testing_patients[]
//   ms_final_state[], ms_time_01[], ms_time_02[], ms_time_12[]
//   ms_censored_01[], ms_censored_02[], ms_censored_12[], ms_os_event_12[]
//   log_cond_surv_01, log_cond_surv_02, log_cond_surv_12_s, log_cond_surv_12_t
//   enable_ms_01, enable_ms_02, enable_ms_12, ms_time_scale_12

// Per-transition log-likelihoods for diagnostic breakdowns.
// The R-side lfo_log_lik_rvar() auto-discovers these via regex
// matches("^patient(_.+)?_log_lik"), so no R changes needed.
array[max_n_rows, max_forecast_horizon] vector[n_all_testing_patients] patient_log_lik_01;
array[max_n_rows, max_forecast_horizon] vector[n_all_testing_patients] patient_log_lik_02;
array[max_n_rows, max_forecast_horizon] vector[n_all_testing_patients] patient_log_lik_12;
array[max_n_rows, max_forecast_horizon] vector[n_all_testing_patients] patient_log_lik;

profile("lfo_gq") {
  for (n in 1:max_n_rows) {
    int n_curr_patients = n_patients - testing_patient_idx[n] + 1;
    int curr_first_testing_patient_idx = n_all_testing_patients - n_curr_patients + 1;
    array[n_curr_patients] int curr_patients =
      last_visit_calendar_day_sort_idx[testing_patient_idx[n]:];

    int m_end = min(n + max_forecast_horizon - 1, n_cutoffs);
    for (m_abs in n:m_end) {
      int m_rel = m_abs - n + 1;
      patient_log_lik_01[n, m_rel] = zeros_vector(n_all_testing_patients);
      patient_log_lik_02[n, m_rel] = zeros_vector(n_all_testing_patients);
      patient_log_lik_12[n, m_rel] = zeros_vector(n_all_testing_patients);
      patient_log_lik[n, m_rel] = zeros_vector(n_all_testing_patients);

      for (i_idx in 1:n_curr_patients) {
        int i = curr_patients[i_idx];
        int patient_idx = curr_first_testing_patient_idx + i_idx - 1;

        int visit_start, visit_end;
        (visit_start, visit_end) = get_pos(patient_visit_pos, i);

        int start_idx = testing_start_idx[n, i];
        int end_idx = m_abs < n_cutoffs ? testing_end_idx[n, m_abs + 1, i] : visit_end;

        if (start_idx <= 0 || end_idx < start_idx || cutoff_last_visit_week[i] <= 0) continue;

        int test_start_week = t_patient_visits[start_idx];
        int test_end_week = t_patient_visits[end_idx];

        real ll_01 = 0;
        real ll_02 = 0;
        real ll_12 = 0;

        // 0→1 transition: conditional log-likelihood for progression.
        // Uses zeros for interval censoring to avoid double-counting IC.
        if (enable_ms_01 && ms_time_01[i] > 0) {
          ll_01 = calc_pch_loglik(
            {ms_time_01[i]},
            {ms_censored_01[i]},
            zeros_int_array(1),
            0,
            log_cond_surv_01[i:i],
            {test_start_week},
            {test_end_week}
          )[1];
        }

        // 0→2 transition: conditional log-likelihood for death without progression.
        // Uses pre-computed ms_censored_02 from modules/multistate/transformed_data.stan.
        if (enable_ms_02 && ms_time_02[i] > 0) {
          ll_02 = calc_pch_loglik(
            {ms_time_02[i]},
            {ms_censored_02[i]},
            zeros_int_array(1),
            0,
            log_cond_surv_02[i:i],
            {test_start_week},
            {test_end_week}
          )[1];
        }

        // 1→2 transition: conditional log-likelihood for post-progression death
        if (enable_ms_12 && ms_final_state[i] == 2 && ms_censored_01[i] == 0) {
          if (ms_time_scale_12 == 0) {
            // Markov: clock-forward time, condition on survival to progression
            int death_or_censor_t = ms_time_12[i] > 0
              ? ms_time_01[i] + ms_time_12[i] : ms_time_01[i];
            ll_12 = calc_pch_loglik(
              {death_or_censor_t},
              {ms_censored_12[i]},
              zeros_int_array(1),
              0,
              log_cond_surv_12_t[i:i],
              {max(1, ms_time_01[i])},
              {test_end_week}
            )[1];
          } else {
            // Semi-Markov: sojourn time from progression
            ll_12 = calc_pch_loglik(
              {ms_time_12[i]},
              {ms_censored_12[i]},
              zeros_int_array(1),
              0,
              log_cond_surv_12_s[i:i],
              {1},
              {test_end_week}
            )[1];
          }
        }

        patient_log_lik_01[n, m_rel, patient_idx] = ll_01;
        patient_log_lik_02[n, m_rel, patient_idx] = ll_02;
        patient_log_lik_12[n, m_rel, patient_idx] = ll_12;
        patient_log_lik[n, m_rel, patient_idx] = ll_01 + ll_02 + ll_12;
      }
    }
  }
}
