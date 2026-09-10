// init/flags.stan
// Initial state module gating flags.
// Multi-level hierarchy: flags are arrays indexed by level (1..n_levels)

int<lower=0,upper=1> enable_pop_cov_init;

// Static compartment: when 1, baseline burden splits decrease/static/growth
// via a nested conditional logit (init_logit_loc unchanged: decrease vs. rest).
int<lower=0,upper=1> enable_static_init;

// Per-level flags for hierarchical intercepts and slopes
array[n_levels] int<lower=0,upper=4> enable_level_intercept_init;
array[n_levels] int<lower=0,upper=1> enable_level_cov_init;
