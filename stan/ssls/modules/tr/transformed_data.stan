// tr/transformed_data.stan
// Compute enabled group counts and position arrays for total rate module.
// These allow parameters to be sized exactly for enabled levels only.

// Enabled group count for tr intercepts
int n_enabled_groups_tr_intercept = compute_n_enabled_groups(
  n_groups_per_level, enable_level_intercept_tr
);

// Position array for enabled intercept levels only
array[n_levels + 1] int enabled_level_pos_tr_intercept = create_enabled_pos(
  n_groups_per_level, enable_level_intercept_tr
);

// Enabled group count for tr slopes
int n_enabled_groups_tr_slope = compute_n_enabled_groups(
  n_groups_per_level, enable_level_cov_tr
);

// Position array for enabled slope levels only
array[n_levels + 1] int enabled_level_pos_tr_slope = create_enabled_pos(
  n_groups_per_level, enable_level_cov_tr
);
