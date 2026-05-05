// ============================================================================
// FLAT AGGREGATION INJECT
// ============================================================================
// Included AFTER burden_endpoints.stan by pioneer.stan.
// Aggregates over all n_forecast_patients without a trial dimension.
//
// Prerequisites (in scope, set by burden_endpoints.stan before this include):
//   orr_sample_response, orr_spop_response
//   sample_target_pfs, ..., spop_os_censored — patient-level endpoint arrays
//   spop_is_dropout, sample_is_dropout
// Also required (from model data):
//   max_all_t, pfs_quantiles, pfs_timepoints, n_pfs_timepoints
//
// Output variables (declared by the model in GQ scope, flat shaped):
//   sample_target_orr, spop_target_orr               — real
//   sample_target_km_est, ..., spop_os_n              — vector[...] (no trial dim)
//   spop_cif_01, spop_cif_02, spop_cif_03            — vector[max_all_t + 1]
//   sample_cif_01, sample_cif_02, sample_cif_03      — vector[max_all_t + 1]

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

(spop_cif_01,   spop_cif_02,   spop_cif_03) = compute_trial_cif(
  spop_pfs,   spop_right_censored,   spop_is_dropout,
  spop_os,   spop_os_censored,   max_all_t);
(sample_cif_01, sample_cif_02, sample_cif_03) = compute_trial_cif(
  sample_pfs, sample_right_censored, sample_is_dropout,
  sample_os, sample_os_censored, max_all_t);
