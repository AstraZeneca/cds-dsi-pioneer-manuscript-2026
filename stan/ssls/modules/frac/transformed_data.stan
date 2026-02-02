// frac/transformed_data.stan
// Compute enabled group counts and position arrays for fraction module.
// These allow parameters to be sized exactly for enabled levels only.

// Enabled group count for frac intercepts
int n_enabled_groups_frac_intercept = compute_n_enabled_groups(
  n_groups_per_level, enable_level_intercept_frac
);

// Position array for enabled intercept levels only
array[n_levels + 1] int enabled_level_pos_frac_intercept = create_enabled_pos(
  n_groups_per_level, enable_level_intercept_frac
);

// Enabled group count for frac slopes
int n_enabled_groups_frac_slope = compute_n_enabled_groups(
  n_groups_per_level, enable_level_cov_frac
);

// Position array for enabled slope levels only
array[n_levels + 1] int enabled_level_pos_frac_slope = create_enabled_pos(
  n_groups_per_level, enable_level_cov_frac
);
