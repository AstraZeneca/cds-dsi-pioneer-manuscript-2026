// frac/transformed_data.stan
// Compute enabled group counts, position arrays, and pre-computed indices for
// efficient vectorized operations in transformed_parameters.

int n_enabled_groups_frac_intercept;
array[n_levels + 1] int enabled_level_pos_frac_intercept;
array[n_patients, n_levels] int patient_frac_intercept_flat_idx;
int n_enabled_groups_frac_slope;
array[n_levels + 1] int enabled_level_pos_frac_slope;
array[n_patients, n_levels] int patient_frac_slope_flat_idx;
(n_enabled_groups_frac_intercept, enabled_level_pos_frac_intercept, patient_frac_intercept_flat_idx,
 n_enabled_groups_frac_slope, enabled_level_pos_frac_slope, patient_frac_slope_flat_idx) =
  compute_level_module_flags(n_patients, n_levels, n_forecast_patients,
    n_forecast_groups_per_level, patient_level_groups,
    enable_level_intercept_frac, enable_level_cov_frac);
