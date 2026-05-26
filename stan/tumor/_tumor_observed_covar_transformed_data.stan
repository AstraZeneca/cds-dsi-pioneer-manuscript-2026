// ============================================================================
// TUMOR (SLD) → SHARED OBSERVED-BURDEN COVARIATE ADAPTER
// ============================================================================
// Aliases the tumor-side names (log_sum_tumor_size, median_log_sld_obs,
// iqr_log_sld_obs) onto the burden-agnostic names expected by
// stan/_observed_covar_transformed_data.stan, then includes the shared
// populator.
//
// Tumor SLD: a visit is "measured" when sum_tumor_size > 0. Zero-valued
// visits represent missing/below-LOD readings — log(0) is -inf, so feeding
// those into the standardization would propagate non-finite values into
// the log probability and stall sampling at initialization.
//
// MUST be included AFTER: modules/tumor/transformed_data.stan
//   (provides log_sum_tumor_size, median_log_sld_obs, iqr_log_sld_obs,
//   sum_tumor_size)
// MUST be included AFTER: modules/multistate/transformed_data.stan
//   (declares ms_obs_visit_covar_flat)

vector[sum(n_patient_visits)] log_burden_obs = log_sum_tumor_size;
real median_log_burden_obs = median_log_sld_obs;
real iqr_log_burden_obs = iqr_log_sld_obs;
array[sum(n_patient_visits)] int burden_measured;
for (v in 1:sum(n_patient_visits)) {
  burden_measured[v] = sum_tumor_size[v] > 0 ? 1 : 0;
}

#include "../_observed_covar_transformed_data.stan"
