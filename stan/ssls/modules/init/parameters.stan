// init/parameters.stan
// Activated Initial proportion (init) module parameter declarations.
// Multi-level hierarchy: parameters use unified level-indexed structure.

// Population intercept (always on)
real init_logit_loc_pop;

// Population covariate coefficients (QR space) — length 0 if disabled
vector[enable_pop_cov_init ? n_covar : 0] init_coef_qr_pop;

// ===== UNIFIED LEVEL STRUCTURE =====

// Intercept SD hyperparameters - one per level
array[n_levels] real<lower=0> init_sd_level_intercept;

// Raw standard normal draws for intercepts - flattened across all levels
vector[n_total_groups] init_raw_level_intercept;

// Slope SD hyperparameters - one vector per level
array[n_levels] vector<lower=0>[n_covar] init_sd_level_slope;

// Raw standard normal draws for slopes - flattened across all levels
matrix[n_total_groups, n_covar] init_raw_level_slope;
