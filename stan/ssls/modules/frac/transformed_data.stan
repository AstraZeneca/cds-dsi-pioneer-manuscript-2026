// frac/transformed_data.stan
// Compute enabled group counts, position arrays, and pre-computed indices for
// efficient vectorized operations in transformed_parameters.

// Enabled group count for frac intercepts
int n_enabled_groups_frac_intercept = compute_n_enabled_groups(
  n_groups_per_level, enable_level_intercept_frac
);

// Position array for enabled intercept levels only
array[n_levels + 1] int enabled_level_pos_frac_intercept = create_enabled_pos(
  n_groups_per_level, enable_level_intercept_frac
);

// Pre-computed flat indices for intercepts
array[n_patients, n_levels] int patient_frac_intercept_flat_idx;
{
  for (i in 1:n_patients) {
    for (lv in 1:n_levels) {
      if (enable_level_intercept_frac[lv]) {
        patient_frac_intercept_flat_idx[i, lv] =
          get_global_group_idx(enabled_level_pos_frac_intercept, lv, patient_level_groups[i, lv]);
      } else {
        patient_frac_intercept_flat_idx[i, lv] = 1;
      }
    }
  }
}

// Enabled group count for frac slopes
int n_enabled_groups_frac_slope = compute_n_enabled_groups(
  n_groups_per_level, enable_level_cov_frac
);

// Position array for enabled slope levels only
array[n_levels + 1] int enabled_level_pos_frac_slope = create_enabled_pos(
  n_groups_per_level, enable_level_cov_frac
);

// Pre-computed flat indices for slopes
array[n_patients, n_levels] int patient_frac_slope_flat_idx;
{
  for (i in 1:n_patients) {
    for (lv in 1:n_levels) {
      if (enable_level_cov_frac[lv]) {
        patient_frac_slope_flat_idx[i, lv] =
          get_global_group_idx(enabled_level_pos_frac_slope, lv, patient_level_groups[i, lv]);
      } else {
        patient_frac_slope_flat_idx[i, lv] = 1;
      }
    }
  }
}
