// init/transformed_data.stan
// Compute enabled group counts, position arrays, and pre-computed indices for
// efficient vectorized operations in transformed_parameters.

// Enabled group count for init intercepts (used in parameters block)
int n_enabled_groups_init_intercept = compute_n_enabled_groups(
  n_hmc_groups_per_level, enable_level_intercept_init
);

// Position array for enabled intercept levels only
array[n_levels + 1] int enabled_level_pos_init_intercept = create_enabled_pos(
  n_hmc_groups_per_level, enable_level_intercept_init
);

// Pre-computed flat indices for intercepts
array[n_patients, n_levels] int patient_init_intercept_flat_idx;
{
  for (i in 1:n_patients) {
    for (lv in 1:n_levels) {
      if (enable_level_intercept_init[lv]) {
        patient_init_intercept_flat_idx[i, lv] =
          (lv == n_levels && patient_level_groups[i, lv] > n_hmc_patients) ? 1
          : get_global_group_idx(enabled_level_pos_init_intercept, lv, patient_level_groups[i, lv]);
      } else {
        patient_init_intercept_flat_idx[i, lv] = 1;
      }
    }
  }
}

// Enabled group count for init slopes
int n_enabled_groups_init_slope = compute_n_enabled_groups(
  n_hmc_groups_per_level, enable_level_cov_init
);

// Position array for enabled slope levels only
array[n_levels + 1] int enabled_level_pos_init_slope = create_enabled_pos(
  n_hmc_groups_per_level, enable_level_cov_init
);

// Pre-computed flat indices for slopes
array[n_patients, n_levels] int patient_init_slope_flat_idx;
{
  for (i in 1:n_patients) {
    for (lv in 1:n_levels) {
      if (enable_level_cov_init[lv]) {
        patient_init_slope_flat_idx[i, lv] =
          (lv == n_levels && patient_level_groups[i, lv] > n_hmc_patients) ? 1
          : get_global_group_idx(enabled_level_pos_init_slope, lv, patient_level_groups[i, lv]);
      } else {
        patient_init_slope_flat_idx[i, lv] = 1;
      }
    }
  }
}
