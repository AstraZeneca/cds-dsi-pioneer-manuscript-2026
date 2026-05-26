// ============================================================================
// OBSERVED-PSA COVARIATE FOR VISIT-GATED 0->1 TRANSITION
// ============================================================================
// Populates ms_obs_visit_covar_flat (declared in modules/multistate/transformed_data.stan)
// from OBSERVED PSA values (data), not from the modeled state-space trajectory.
//
// No carry-forward: the value at unmeasured visits is 0 (placeholder), but
// unmeasured visits are skipped by the likelihood (psa_measured[v] check).
//
// MUST be included AFTER: modules/psa/transformed_data.stan
//   (needs log_psa_values, psa_measured, median_log_psa_obs, iqr_log_psa_obs)
// MUST be included AFTER: modules/multistate/transformed_data.stan
//   (declares ms_obs_visit_covar_flat)

if (enable_ms_visit_gated_01 && !enable_ms_visit_gated_latent_01) {
  for (v in 1:sum(n_patient_visits)) {
    ms_obs_visit_covar_flat[v] = psa_measured[v]
      ? (log_psa_values[v] - median_log_psa_obs) / iqr_log_psa_obs
      : 0.0;
  }
}
