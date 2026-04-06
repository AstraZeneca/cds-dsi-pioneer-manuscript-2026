// init/parameters.stan
// Activated Initial proportion (init) module parameter declarations.
// Multi-level hierarchy: parameters use unified level-indexed structure.
// Parameters are sized by enabled group counts to avoid wasted sampling.

// Population intercept (always on)
real init_logit_loc_pop;

// Population covariate coefficients (QR space) — length 0 if disabled
vector[enable_pop_cov_init ? n_covar : 0] init_coef_qr_pop;

// ===== UNIFIED LEVEL STRUCTURE =====

// Intercept SD free parameters - sized to RE levels only (mode=2)
// FE levels (mode=1) use init_fe_sd_level_intercept from hyperparams; disabled levels use 0.
// Full n_levels array init_sd_level_intercept is assembled in transformed_parameters.
array[n_re_levels_init_intercept] real<lower=0> init_sd_level_intercept_raw;

// Raw standard normal draws for intercepts - sized by ENABLED groups only
vector[n_enabled_groups_init_intercept] init_raw_level_intercept;

// Slope SD hyperparameters - one vector per level (always n_levels for simplicity)
array[n_levels] vector<lower=0>[n_covar] init_sd_level_slope;

// Raw standard normal draws for slopes - sized by ENABLED groups only
matrix[n_enabled_groups_init_slope, n_covar] init_raw_level_slope;
