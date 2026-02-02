// init/transformed_data.stan
// Compute enabled group counts and position arrays for initial state module.
// These allow parameters to be sized exactly for enabled levels only.

// Enabled group count for init intercepts
int n_enabled_groups_init_intercept = compute_n_enabled_groups(
  n_groups_per_level, enable_level_intercept_init
);

// Position array for enabled intercept levels only
array[n_levels + 1] int enabled_level_pos_init_intercept = create_enabled_pos(
  n_groups_per_level, enable_level_intercept_init
);

// Enabled group count for init slopes
int n_enabled_groups_init_slope = compute_n_enabled_groups(
  n_groups_per_level, enable_level_cov_init
);

// Position array for enabled slope levels only
array[n_levels + 1] int enabled_level_pos_init_slope = create_enabled_pos(
  n_groups_per_level, enable_level_cov_init
);
