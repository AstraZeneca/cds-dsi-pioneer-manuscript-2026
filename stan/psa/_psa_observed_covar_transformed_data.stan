// ============================================================================
// PSA → SHARED OBSERVED-BURDEN COVARIATE ADAPTER
// ============================================================================
// Aliases the PSA-side names (log_psa_values, median_log_psa_obs,
// iqr_log_psa_obs, psa_measured) onto the burden-agnostic names expected by
// stan/_observed_covar_transformed_data.stan, then includes the shared
// populator.
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
