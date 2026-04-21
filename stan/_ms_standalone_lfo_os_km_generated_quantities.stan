// ============================================================================
// MS-Standalone LFO: Posterior-predictive OS KM from re-censored events
// ============================================================================
// Wraps compute_ms_os_km_rng() with the lfo_ms_* arrays from
// _ms_standalone_lfo_transformed_data.stan. Post-cutoff events are hidden
// (re-censored), so sample_os_km_est at this cutoff represents the model's
// out-of-sample forecast for the horizon beyond the cutoff while preserving
// known pre-cutoff events.
//
// Requires in scope (declared in transformed_data / other GQ includes):
//   forecast_patient_idx, forecast_trial_patient_pos
//   n_total_visits, n_total_forecast_obs_visits, forecast_observation_interval
//   max_all_t, log_cond_surv_01/02/03/12_s/12_t/32
//   enable_ms_02/03/32/12, ms_time_scale_12
//   lfo_ms_final_state, lfo_ms_time_01, lfo_ms_censored_01,
//   lfo_ms_time_02, lfo_ms_time_12, lfo_ms_time_03, lfo_ms_time_32
//   ms_os_event_12 — unchanged (recensor_ms_at_cutoff does not return it)
//   patient_visit_pos, forecast_visits_pos, forecast_obs_visits_pos,
//   patient_last_obs_visit, last_predict_visit, t_patient_visits,
//   n_patient_screening_visits

// Derive re-censored 02/12/32 indicators from re-censored state/time arrays.
// Mirrors multistate/transformed_data.stan but against lfo_ms_*.
array[n_patients] int lfo_ms_censored_02;
array[n_patients] int lfo_ms_censored_12;
array[n_patients] int lfo_ms_censored_32;
(lfo_ms_censored_02, lfo_ms_censored_12, lfo_ms_censored_32) =
  derive_ms_censoring_indicators(
    lfo_ms_final_state, lfo_ms_time_01, lfo_ms_time_32
  );

// ms_os_event_12 is the exact calendar week of observed 1→2 deaths; for
// patients whose 1→2 event moved past the cutoff, lfo_ms_censored_12 is
// already 1, so derive_sample_os_rng ignores ms_os_event_12 for them.
array[n_trials] vector<lower=0, upper=1>[max_all_t + 1] sample_os_km_est =
  compute_ms_os_km_rng(
    forecast_patient_idx,
    forecast_trial_patient_pos,
    max_all_t,
    n_total_visits,
    n_total_forecast_obs_visits,
    forecast_observation_interval,
    log_cond_surv_01,
    enable_ms_02, enable_ms_03, enable_ms_32, enable_ms_12,
    ms_time_scale_12,
    log_cond_surv_02, log_cond_surv_03, log_cond_surv_32,
    log_cond_surv_12_s, log_cond_surv_12_t,
    lfo_ms_final_state,
    lfo_ms_time_01,
    lfo_ms_censored_01,
    lfo_ms_time_02,
    lfo_ms_censored_02,
    lfo_ms_time_12,
    lfo_ms_censored_12,
    ms_os_event_12,
    lfo_ms_time_03,
    lfo_ms_time_32,
    lfo_ms_censored_32,
    patient_visit_pos,
    forecast_visits_pos,
    forecast_obs_visits_pos,
    patient_last_obs_visit,
    last_predict_visit,
    t_patient_visits,
    n_patient_screening_visits
  );
