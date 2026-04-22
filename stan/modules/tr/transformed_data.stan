// tr/transformed_data.stan
// Compute enabled group counts, position arrays, and pre-computed indices for
// efficient vectorized operations in transformed_parameters.

int n_enabled_groups_tr_intercept;
array[n_levels + 1] int enabled_level_pos_tr_intercept;
array[n_patients, n_levels] int patient_tr_intercept_flat_idx;
int n_enabled_groups_tr_slope;
array[n_levels + 1] int enabled_level_pos_tr_slope;
array[n_patients, n_levels] int patient_tr_slope_flat_idx;
(n_enabled_groups_tr_intercept, enabled_level_pos_tr_intercept, patient_tr_intercept_flat_idx,
 n_enabled_groups_tr_slope, enabled_level_pos_tr_slope, patient_tr_slope_flat_idx) =
  compute_level_module_flags(n_patients, n_levels, n_forecast_patients,
    n_forecast_groups_per_level, patient_level_groups,
    enable_level_intercept_tr, enable_level_cov_tr);

// Count RE-only levels for right-sizing the SD parameter array
int n_re_levels_tr_intercept = 0;
for (lv in 1:n_levels)
  if (enable_level_intercept_tr[lv] == LEVEL_MODE_RE) n_re_levels_tr_intercept += 1;
