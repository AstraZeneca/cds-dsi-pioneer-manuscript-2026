// ============================================================================
// TUMOR (SLD) → SHARED OBSERVED-BURDEN COVARIATE ADAPTER
// ============================================================================
// Two responsibilities:
//
//   (1) Declare the burden-agnostic aliases that downstream shared includes
//       expect — `median_log_burden_obs`, `iqr_log_burden_obs` — using the
//       tumor-side data values (`median_log_sld_obs`, `iqr_log_sld_obs`).
//       These aliases are referenced from `transformed parameters` by
//       `_ms_burden_tv_covar.stan` even when the observed populator below
//       is gated off, so this file is REQUIRED for the tumor model whether
//       visit-gating runs in observed or latent mode.
//
//   (2) When `enable_ms_visit_gated_01 && !enable_ms_visit_gated_latent_01`
//       (observed visit-gated mode), populate `ms_obs_visit_covar_flat`
//       from the standardized observed log-SLD values. In latent visit-
//       gated mode this populator is a no-op — the shared file's `if`
//       gate skips it — and the multistate likelihood reads the latent
//       trajectory from `ms_time_varying_covar_01` instead.
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

// Build burden_measured mask + a finite log_burden_obs vector. We can't
// reuse log_sum_tumor_size directly because it contains log(0) = -inf at
// zero / below-LOD visits — even though those entries are masked out in
// the standardization, Stan's autodiff still evaluates the masked branch
// of the ternary and propagates -inf into the log probability gradient.
vector[sum(n_patient_visits)] log_burden_obs;
array[sum(n_patient_visits)] int burden_measured;
for (v in 1:sum(n_patient_visits)) {
  if (sum_tumor_size[v] > 0) {
    burden_measured[v] = 1;
    log_burden_obs[v] = log(sum_tumor_size[v]);
  } else {
    burden_measured[v] = 0;
    log_burden_obs[v] = 0.0;  // placeholder; never read because mask is 0
  }
}
real median_log_burden_obs = median_log_sld_obs;
real iqr_log_burden_obs = iqr_log_sld_obs;
// Velocity standardization aliases for the (level, velocity) coupling basis.
// The burden-agnostic names are read by _ms_burden_tv_covar.stan /
// _ms_burden_inline_tv_covar.stan when enable_ms_velocity_basis = 1. PSA does
// not alias these (it stays on the legacy 3-feature basis).
real median_velocity_burden_obs = median_velocity_obs;
real iqr_velocity_burden_obs = iqr_velocity_obs;

#include "../_observed_covar_transformed_data.stan"
