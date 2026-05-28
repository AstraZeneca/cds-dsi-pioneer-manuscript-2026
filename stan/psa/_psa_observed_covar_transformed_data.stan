// ============================================================================
// PSA → SHARED OBSERVED-BURDEN COVARIATE ADAPTER
// ============================================================================
// Two responsibilities (see _tumor_observed_covar_transformed_data.stan for
// the symmetric tumor adapter):
//
//   (1) Declare the burden-agnostic aliases that downstream shared includes
//       expect — `log_burden_obs`, `median_log_burden_obs`, `iqr_log_burden_obs`,
//       `burden_measured` — using the PSA-side data values.
//
//   (2) When `enable_ms_visit_gated_01 && !enable_ms_visit_gated_latent_01`
//       (observed visit-gated mode), populate `ms_obs_visit_covar_flat`
//       from the standardized observed log-PSA values. In latent visit-
//       gated mode the populator below is a no-op.
//
// MUST be included AFTER: modules/psa/transformed_data.stan
//   (provides log_psa_values, psa_measured, median_log_psa_obs, iqr_log_psa_obs)
// MUST be included AFTER: modules/multistate/transformed_data.stan
//   (declares ms_obs_visit_covar_flat)

vector[sum(n_patient_visits)] log_burden_obs = log_psa_values;
real median_log_burden_obs = median_log_psa_obs;
real iqr_log_burden_obs = iqr_log_psa_obs;
array[sum(n_patient_visits)] int burden_measured = psa_measured;

#include "../_observed_covar_transformed_data.stan"
