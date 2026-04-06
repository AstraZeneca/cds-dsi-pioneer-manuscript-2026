// init/transformed_data.stan
// Compute enabled group counts, position arrays, and pre-computed indices for
// efficient vectorized operations in transformed_parameters.

int n_enabled_groups_init_intercept;
array[n_levels + 1] int enabled_level_pos_init_intercept;
array[n_patients, n_levels] int patient_init_intercept_flat_idx;
int n_enabled_groups_init_slope;
array[n_levels + 1] int enabled_level_pos_init_slope;
array[n_patients, n_levels] int patient_init_slope_flat_idx;
(n_enabled_groups_init_intercept, enabled_level_pos_init_intercept, patient_init_intercept_flat_idx,
 n_enabled_groups_init_slope, enabled_level_pos_init_slope, patient_init_slope_flat_idx) =
  compute_level_module_flags(n_patients, n_levels, n_forecast_patients,
    n_forecast_groups_per_level, patient_level_groups,
    enable_level_intercept_init, enable_level_cov_init);

// Count RE-only levels for right-sizing the SD parameter array
int n_re_levels_init_intercept = 0;
for (lv in 1:n_levels)
  if (enable_level_intercept_init[lv] == LEVEL_MODE_RE) n_re_levels_init_intercept += 1;
