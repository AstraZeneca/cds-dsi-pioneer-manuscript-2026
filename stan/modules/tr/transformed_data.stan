// tr/transformed_data.stan
// Compute enabled group counts, position arrays, and pre-computed indices for
// efficient vectorized operations in transformed_parameters.

// Enabled group count for tr intercepts
int n_enabled_groups_tr_intercept = compute_n_enabled_groups(
  n_forecast_groups_per_level, enable_level_intercept_tr
);

// Position array for enabled intercept levels only
array[n_levels + 1] int enabled_level_pos_tr_intercept = create_enabled_pos(
  n_forecast_groups_per_level, enable_level_intercept_tr
);

// Pre-computed flat indices: patient_tr_intercept_flat_idx[i, lv] gives the
// direct index into tr_raw_level_intercept for patient i at level lv.
// This avoids array slicing in transformed_parameters.
array[n_patients, n_levels] int patient_tr_intercept_flat_idx;
{
  for (i in 1:n_patients) {
    for (lv in 1:n_levels) {
      if (enable_level_intercept_tr[lv]) {
        int gid = patient_level_groups[i, lv];
        // Patient level (lv == n_levels): RWD patients have gid > n_forecast_patients.
        // Their patient-level entry is dummy=1; patient level is skipped in Laplace.
        patient_tr_intercept_flat_idx[i, lv] =
          (lv == n_levels && gid > n_forecast_patients) ? 1
          : get_global_group_idx(enabled_level_pos_tr_intercept, lv, gid);
      } else {
        patient_tr_intercept_flat_idx[i, lv] = 1;
      }
    }
  }
}

// Enabled group count for tr slopes
int n_enabled_groups_tr_slope = compute_n_enabled_groups(
  n_forecast_groups_per_level, enable_level_cov_tr
);

// Position array for enabled slope levels only
array[n_levels + 1] int enabled_level_pos_tr_slope = create_enabled_pos(
  n_forecast_groups_per_level, enable_level_cov_tr
);

// Pre-computed flat indices for slopes
array[n_patients, n_levels] int patient_tr_slope_flat_idx;
{
  for (i in 1:n_patients) {
    for (lv in 1:n_levels) {
      if (enable_level_cov_tr[lv]) {
        int gid = patient_level_groups[i, lv];
        patient_tr_slope_flat_idx[i, lv] =
          (lv == n_levels && gid > n_forecast_patients) ? 1
          : get_global_group_idx(enabled_level_pos_tr_slope, lv, gid);
      } else {
        patient_tr_slope_flat_idx[i, lv] = 1;  // dummy (level disabled)
      }
    }
  }
}
