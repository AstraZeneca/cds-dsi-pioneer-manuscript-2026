// frac/flags.stan
// Fraction mix module gating flags.
// Multi-level hierarchy: flags are arrays indexed by level (1..n_levels)

int<lower=0,upper=1> enable_pop_cov_frac;

// Per-level flags for hierarchical intercepts and slopes
array[n_levels] int<lower=0,upper=4> enable_level_intercept_frac;
array[n_levels] int<lower=0,upper=1> enable_level_cov_frac;
