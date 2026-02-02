// frac/parameters.stan
// Activated Fraction (frac) module parameter declarations.
// Multi-level hierarchy: parameters use unified level-indexed structure.

// Population intercept (always on)
real frac_logit_loc_pop;

// Population covariate coefficients (QR space) — length 0 if disabled
vector[enable_pop_cov_frac ? n_covar : 0] frac_coef_qr_pop;

// ===== UNIFIED LEVEL STRUCTURE =====

// Intercept SD hyperparameters - one per level
array[n_levels] real<lower=0> frac_sd_level_intercept;

// Raw standard normal draws for intercepts - flattened across all levels
vector[n_total_groups] frac_raw_level_intercept;

// Slope SD hyperparameters - one vector per level
array[n_levels] vector<lower=0>[n_covar] frac_sd_level_slope;

// Raw standard normal draws for slopes - flattened across all levels
matrix[n_total_groups, n_covar] frac_raw_level_slope;
