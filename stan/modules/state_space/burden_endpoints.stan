// ============================================================================
// BURDEN ENDPOINTS — SHARED COMPUTATION
// ============================================================================
// Included inside a local { } block. Contract interface variables and endpoint
// output variables must already be in scope.
//
// Contract interface (local, set by model before this include):
//   obs_biomarker_cat          — array[n_total_visits] int
//   rep_biomarker_cat          — array[n_total_visits] int
//   forecast_obs_biomarker_cat — array[n_total_forecast_obs_visits] int
//   burden_pfs                 — array[n_patients] int
//   burden_right_censored      — array[n_patients] int
//   burden_target_pfs          — array[n_patients] int
//   burden_target_right_censored — array[n_patients] int
//
// Output variables (GQ scope, declared by model before the { } block):
//   sample_target_pfs, spop_target_pfs, ..., spop_cif_01, etc.

profile("burden_endpoints") {
  (sample_target_pfs, sample_target_right_censored,
   spop_target_pfs, spop_target_right_censored,
   spop_target_obs_cens_pfs, spop_target_obs_cens_right_censored,
   sample_ms_pfs, sample_ms_right_censored,
   spop_ms_pfs, spop_ms_right_censored,
   sample_pfs, sample_right_censored,
   spop_pfs, spop_right_censored,
   sample_target_confirmed_response, sample_target_unconfirmed_response,
   spop_target_confirmed_response, spop_target_unconfirmed_response,
   forecast_target_pfs, forecast_target_right_censored,
   sample_os, sample_os_censored,
   spop_os, spop_os_censored,
   spop_is_dropout, sample_is_dropout
  ) = calculate_all_patients_endpoints_rng(
    hmc_patient_idx,
    obs_biomarker_cat,
    rep_biomarker_cat,
    forecast_obs_biomarker_cat,
    forecast_obs_visits_pos,
    forecast_observation_interval,
    log_cond_surv_01,
    enable_ms_02,
    enable_ms_03,
    enable_ms_32,
    enable_ms_12,
    ms_time_scale_12,
    log_cond_surv_02,
    log_cond_surv_03,
    log_cond_surv_32,
    log_cond_surv_12_s,
    log_cond_surv_12_t,
    burden_pfs,
    burden_right_censored,
    burden_target_pfs,
    burden_target_right_censored,
    ms_censored_01,
    ms_final_state,
    ms_time_02,
    ms_censored_02,
    ms_censored_12,
    ms_time_12,
    ms_os_event_12,
    ms_time_03,
    ms_time_32,
    ms_censored_32,
    patient_visit_pos,
    forecast_visits_pos,
    patient_last_obs_visit,
    last_predict_visit,
    t_patient_visits,
    max_all_t,
    n_patient_screening_visits
  );
}

// Aggregate to trial-level metrics
(sample_target_orr, spop_target_orr,
 sample_target_km_est, spop_target_km_est, spop_target_obs_cens_km_est,
 sample_ms_km_est, spop_ms_km_est,
 sample_km_est, spop_km_est,
 sample_target_quant_pfs, spop_target_quant_pfs,
 sample_target_quant_pfs_exceeds_max, spop_target_quant_pfs_exceeds_max,
 sample_ms_quant_pfs, spop_ms_quant_pfs,
 sample_ms_quant_pfs_exceeds_max, spop_ms_quant_pfs_exceeds_max,
 sample_quant_pfs, spop_quant_pfs,
 sample_quant_pfs_exceeds_max, spop_quant_pfs_exceeds_max,
 sample_target_pfs_n, spop_target_pfs_n,
 sample_ms_pfs_n, spop_ms_pfs_n,
 sample_pfs_n, spop_pfs_n,
 sample_os_km_est, spop_os_km_est,
 sample_os_quant, spop_os_quant,
 sample_os_quant_exceeds_max, spop_os_quant_exceeds_max,
 sample_os_n, spop_os_n) =
  aggregate_trial_metrics(
    sample_target_confirmed_response,
    spop_target_confirmed_response,
    sample_target_pfs,
    sample_target_right_censored,
    spop_target_pfs,
    spop_target_right_censored,
    spop_target_obs_cens_pfs,
    spop_target_obs_cens_right_censored,
    sample_ms_pfs,
    sample_ms_right_censored,
    spop_ms_pfs,
    spop_ms_right_censored,
    sample_pfs,
    sample_right_censored,
    spop_pfs,
    spop_right_censored,
    sample_os,
    sample_os_censored,
    spop_os,
    spop_os_censored,
    hmc_trial_patient_pos,
    max_all_t,
    pfs_quantiles,
    pfs_timepoints
  );

// Aggregate to conditional group-level metrics
(cond_sample_target_orr, cond_spop_target_orr,
 cond_sample_target_km_est, cond_spop_target_km_est, cond_spop_target_obs_cens_km_est,
 cond_sample_ms_km_est, cond_spop_ms_km_est,
 cond_sample_km_est, cond_spop_km_est,
 cond_sample_target_quant_pfs, cond_spop_target_quant_pfs,
 cond_sample_target_quant_pfs_exceeds_max, cond_spop_target_quant_pfs_exceeds_max,
 cond_sample_ms_quant_pfs, cond_spop_ms_quant_pfs,
 cond_sample_ms_quant_pfs_exceeds_max, cond_spop_ms_quant_pfs_exceeds_max,
 cond_sample_quant_pfs, cond_spop_quant_pfs,
 cond_sample_quant_pfs_exceeds_max, cond_spop_quant_pfs_exceeds_max,
 cond_sample_target_pfs_n, cond_spop_target_pfs_n,
 cond_sample_ms_pfs_n, cond_spop_ms_pfs_n,
 cond_sample_pfs_n, cond_spop_pfs_n,
 cond_sample_os_km_est, cond_spop_os_km_est,
 cond_sample_os_quant, cond_spop_os_quant,
 cond_sample_os_quant_exceeds_max, cond_spop_os_quant_exceeds_max,
 cond_sample_os_n, cond_spop_os_n) =
  aggregate_conditional_group_metrics(
    sample_target_confirmed_response,
    spop_target_confirmed_response,
    sample_target_pfs,
    sample_target_right_censored,
    spop_target_pfs,
    spop_target_right_censored,
    spop_target_obs_cens_pfs,
    spop_target_obs_cens_right_censored,
    sample_ms_pfs,
    sample_ms_right_censored,
    spop_ms_pfs,
    spop_ms_right_censored,
    sample_pfs,
    sample_right_censored,
    spop_pfs,
    spop_right_censored,
    sample_os,
    sample_os_censored,
    spop_os,
    spop_os_censored,
    cond_group,
    cond_group_pos,
    max_all_t,
    pfs_quantiles,
    pfs_timepoints
  );

// Competing Risks CIF (per-trial, empirical subdistribution)
// Uses hmc_trial_patient_pos — endpoint arrays are n_hmc_patients-sized,
// not n_patients-sized (differs when RWD patients are Laplace-marginalized).
for (s in 1:n_trials) {
  int n_tr = get_pos_size(hmc_trial_patient_pos, s);
  if (n_tr > 0) {
    int tr_start; int tr_end;
    (tr_start, tr_end) = get_pos(hmc_trial_patient_pos, s);

    (spop_cif_01[s], spop_cif_02[s], spop_cif_03[s]) = compute_trial_cif(
      spop_pfs[tr_start:tr_end], spop_right_censored[tr_start:tr_end],
      spop_is_dropout[tr_start:tr_end],
      spop_os[tr_start:tr_end], spop_os_censored[tr_start:tr_end],
      max_all_t
    );

    (sample_cif_01[s], sample_cif_02[s], sample_cif_03[s]) = compute_trial_cif(
      sample_pfs[tr_start:tr_end], sample_right_censored[tr_start:tr_end],
      sample_is_dropout[tr_start:tr_end],
      sample_os[tr_start:tr_end], sample_os_censored[tr_start:tr_end],
      max_all_t
    );
  }
}
