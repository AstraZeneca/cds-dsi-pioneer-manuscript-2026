// ============================================================================
// TRIAL-INDEXED AGGREGATION INJECT
// ============================================================================
// Included AFTER burden_endpoints.stan by tumor and ms-standalone models.
// Provides trial-level aggregate metrics using aggregate_trial_metrics() and
// per-trial CIF via forecast_trial_patient_pos.
//
// Prerequisites (in scope, set by burden_endpoints.stan before this include):
//   orr_sample_response, orr_spop_response
//   sample_target_pfs, ..., spop_os_censored — patient-level endpoint arrays
//   spop_is_dropout, sample_is_dropout
// Also required (from model data/transformed data):
//   forecast_trial_patient_pos, n_trials, max_all_t, pfs_quantiles, pfs_timepoints
//
// Output variables (declared by the model in GQ scope, array[n_trials] shaped):
//   sample_target_orr, spop_target_orr               — vector[n_trials]
//   sample_target_km_est, ..., spop_os_n              — array[n_trials] vector[...]
//   spop_cif_01, spop_cif_02, spop_cif_03            — array[n_trials] vector[...]
//   sample_cif_01, sample_cif_02, sample_cif_03      — array[n_trials] vector[...]

(sample_target_orr, spop_target_orr,
 sample_target_km_est, spop_target_km_est, spop_target_obs_cens_km_est,
 sample_ms_pfs_km_est, spop_ms_pfs_km_est,
 sample_pfs_km_est, spop_pfs_km_est,
 sample_target_pfs_quant, spop_target_pfs_quant,
 sample_target_pfs_quant_exceeds_max, spop_target_pfs_quant_exceeds_max,
 sample_ms_pfs_quant, spop_ms_pfs_quant,
 sample_ms_pfs_quant_exceeds_max, spop_ms_pfs_quant_exceeds_max,
 sample_pfs_quant, spop_pfs_quant,
 sample_pfs_quant_exceeds_max, spop_pfs_quant_exceeds_max,
 sample_target_pfs_n, spop_target_pfs_n,
 sample_ms_pfs_n, spop_ms_pfs_n,
 sample_pfs_n, spop_pfs_n,
 sample_os_km_est, spop_os_km_est,
 sample_os_quant, spop_os_quant,
 sample_os_quant_exceeds_max, spop_os_quant_exceeds_max,
 sample_os_n, spop_os_n) =
  aggregate_trial_metrics(
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
    forecast_trial_patient_pos,
    max_all_t,
    pfs_quantiles,
    pfs_timepoints
  );

// Competing Risks CIF (per-trial, empirical subdistribution)
for (s in 1:n_trials) {
  int n_tr = get_pos_size(forecast_trial_patient_pos, s);
  if (n_tr > 0) {
    int tr_start; int tr_end;
    (tr_start, tr_end) = get_pos(forecast_trial_patient_pos, s);

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
