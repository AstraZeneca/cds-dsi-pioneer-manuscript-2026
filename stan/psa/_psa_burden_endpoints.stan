// ============================================================================
// BURDEN ENDPOINTS — PSA FLAT AGGREGATION
// ============================================================================
// Flat (no trial dimension) version of modules/state_space/burden_endpoints.stan.
// Used by pioneer.stan only. All other models use burden_endpoints.stan.
//
// Contract interface: identical to burden_endpoints.stan.
// GQ output variables in scope must be flat (not array[n_trials]):
//   sample_target_orr, spop_target_orr                 — real
//   sample_target_km_est, spop_target_km_est, ...      — vector[max_all_t + 1]
//   sample_target_pfs_quant, ...                       — vector[n_pfs_quantiles]
//   sample_target_pfs_quant_exceeds_max, ...           — array[n_pfs_quantiles] int
//   sample_target_pfs_n, ...                           — vector[n_pfs_timepoints]
//   spop_cif_01, spop_cif_02, spop_cif_03              — vector[max_all_t + 1]
//   sample_cif_01, sample_cif_02, sample_cif_03        — vector[max_all_t + 1]

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
    forecast_patient_idx,
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
    n_patient_screening_visits,
    burden_enable_ms_visit_gated_01,
    burden_tv_coef_01_val,
    burden_forecast_obs_log_psa,
    burden_median_log_psa_obs,
    burden_iqr_log_psa_obs
  );
}

array[n_forecast_patients] int orr_sample_response =
  burden_orr_use_confirmed_response
    ? sample_target_confirmed_response
    : sample_target_unconfirmed_response;
array[n_forecast_patients] int orr_spop_response =
  burden_orr_use_confirmed_response
    ? spop_target_confirmed_response
    : spop_target_unconfirmed_response;

// ── Flat aggregation over all forecast patients ──────────────────────────────
profile("flat_aggregate_metrics") {
  sample_target_orr = mean(orr_sample_response);
  spop_target_orr   = mean(orr_spop_response);

  sample_target_km_est        = estimate_kaplan_meier(sample_target_pfs, sample_target_right_censored, max_all_t).1;
  spop_target_km_est          = estimate_kaplan_meier(spop_target_pfs, spop_target_right_censored, max_all_t, 0).1;
  spop_target_obs_cens_km_est = estimate_kaplan_meier(spop_target_obs_cens_pfs, spop_target_obs_cens_right_censored, max_all_t, 0).1;
  sample_ms_pfs_km_est        = estimate_kaplan_meier(sample_ms_pfs, sample_ms_right_censored, max_all_t, 0).1;
  spop_ms_pfs_km_est          = estimate_kaplan_meier(spop_ms_pfs, spop_ms_right_censored, max_all_t, 0).1;
  sample_pfs_km_est           = estimate_kaplan_meier(sample_pfs, sample_right_censored, max_all_t, 0).1;
  spop_pfs_km_est             = estimate_kaplan_meier(spop_pfs, spop_right_censored, max_all_t, 0).1;
  sample_os_km_est            = estimate_kaplan_meier(sample_os, sample_os_censored, max_all_t, 0).1;
  spop_os_km_est              = estimate_kaplan_meier(spop_os, spop_os_censored, max_all_t, 0).1;

  (sample_target_pfs_quant, sample_target_pfs_quant_exceeds_max) = km_quantiles(sample_target_km_est, pfs_quantiles);
  (spop_target_pfs_quant,   spop_target_pfs_quant_exceeds_max)   = km_quantiles(spop_target_km_est,   pfs_quantiles);
  (sample_ms_pfs_quant,     sample_ms_pfs_quant_exceeds_max)     = km_quantiles(sample_ms_pfs_km_est, pfs_quantiles);
  (spop_ms_pfs_quant,       spop_ms_pfs_quant_exceeds_max)       = km_quantiles(spop_ms_pfs_km_est,   pfs_quantiles);
  (sample_pfs_quant,        sample_pfs_quant_exceeds_max)        = km_quantiles(sample_pfs_km_est,    pfs_quantiles);
  (spop_pfs_quant,          spop_pfs_quant_exceeds_max)          = km_quantiles(spop_pfs_km_est,      pfs_quantiles);
  (sample_os_quant,         sample_os_quant_exceeds_max)         = km_quantiles(sample_os_km_est,     pfs_quantiles);
  (spop_os_quant,           spop_os_quant_exceeds_max)           = km_quantiles(spop_os_km_est,       pfs_quantiles);

  for (n in 1:n_pfs_timepoints) {
    sample_target_pfs_n[n] = calc_km_pfs_n(sample_target_km_est, months_to_weeks(pfs_timepoints[n]));
    spop_target_pfs_n[n]   = calc_km_pfs_n(spop_target_km_est,   months_to_weeks(pfs_timepoints[n]));
    sample_ms_pfs_n[n]     = calc_km_pfs_n(sample_ms_pfs_km_est, months_to_weeks(pfs_timepoints[n]));
    spop_ms_pfs_n[n]       = calc_km_pfs_n(spop_ms_pfs_km_est,   months_to_weeks(pfs_timepoints[n]));
    sample_pfs_n[n]        = calc_km_pfs_n(sample_pfs_km_est,    months_to_weeks(pfs_timepoints[n]));
    spop_pfs_n[n]          = calc_km_pfs_n(spop_pfs_km_est,      months_to_weeks(pfs_timepoints[n]));
    sample_os_n[n]         = calc_km_pfs_n(sample_os_km_est,     months_to_weeks(pfs_timepoints[n]));
    spop_os_n[n]           = calc_km_pfs_n(spop_os_km_est,       months_to_weeks(pfs_timepoints[n]));
  }
}

// ── Conditional group aggregation (unchanged from burden_endpoints.stan) ─────
(cond_sample_target_orr, cond_spop_target_orr,
 cond_sample_target_km_est, cond_spop_target_km_est, cond_spop_target_obs_cens_km_est,
 cond_sample_ms_pfs_km_est, cond_spop_ms_pfs_km_est,
 cond_sample_pfs_km_est, cond_spop_pfs_km_est,
 cond_sample_target_pfs_quant, cond_spop_target_pfs_quant,
 cond_sample_target_pfs_quant_exceeds_max, cond_spop_target_pfs_quant_exceeds_max,
 cond_sample_ms_pfs_quant, cond_spop_ms_pfs_quant,
 cond_sample_ms_pfs_quant_exceeds_max, cond_spop_ms_pfs_quant_exceeds_max,
 cond_sample_pfs_quant, cond_spop_pfs_quant,
 cond_sample_pfs_quant_exceeds_max, cond_spop_pfs_quant_exceeds_max,
 cond_sample_target_pfs_n, cond_spop_target_pfs_n,
 cond_sample_ms_pfs_n, cond_spop_ms_pfs_n,
 cond_sample_pfs_n, cond_spop_pfs_n,
 cond_sample_os_km_est, cond_spop_os_km_est,
 cond_sample_os_quant, cond_spop_os_quant,
 cond_sample_os_quant_exceeds_max, cond_spop_os_quant_exceeds_max,
 cond_sample_os_n, cond_spop_os_n) =
  aggregate_conditional_group_metrics(
    orr_sample_response,
    orr_spop_response,
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

// ── Flat CIF (empirical subdistribution for all forecast patients) ────────────
(spop_cif_01,   spop_cif_02,   spop_cif_03) = compute_trial_cif(
  spop_pfs,   spop_right_censored,   spop_is_dropout,
  spop_os,   spop_os_censored,   max_all_t);
(sample_cif_01, sample_cif_02, sample_cif_03) = compute_trial_cif(
  sample_pfs, sample_right_censored, sample_is_dropout,
  sample_os, sample_os_censored, max_all_t);
