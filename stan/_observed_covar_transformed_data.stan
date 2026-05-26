// ============================================================================
// OBSERVED-BURDEN COVARIATE FOR VISIT-GATED 0->1 TRANSITION (BURDEN-AGNOSTIC)
// ============================================================================
// Populates ms_obs_visit_covar_flat (declared in modules/multistate/transformed_data.stan)
// from observed disease-burden values (data), not from a modeled trajectory.
//
// The caller MUST have populated, before including this file:
//   - vector[sum(n_patient_visits)] log_burden_obs;
//   - real median_log_burden_obs;
//   - real iqr_log_burden_obs;
//   - array[sum(n_patient_visits)] int burden_measured;
//
// PSA models alias these from log_psa_values / median_log_psa_obs /
// iqr_log_psa_obs / psa_measured. Tumor models alias from log_sum_tumor_size /
// median_log_sld_obs / iqr_log_sld_obs and a unit mask (every visit measured).
//
// MUST be included AFTER: modules/multistate/transformed_data.stan
//   (declares ms_obs_visit_covar_flat)
//
// No carry-forward: at unmeasured visits the covariate is 0 (placeholder), and
// those visits are skipped by the visit-gated likelihood loop anyway.

if (enable_ms_visit_gated_01 && !enable_ms_visit_gated_latent_01) {
  for (v in 1:sum(n_patient_visits)) {
    ms_obs_visit_covar_flat[v] = burden_measured[v]
      ? (log_burden_obs[v] - median_log_burden_obs) / iqr_log_burden_obs
      : 0.0;
  }
}
