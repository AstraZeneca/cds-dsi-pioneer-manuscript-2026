// frac/parameters.stan
// Activated Fraction (frac) module parameter declarations.
// Multi-level hierarchy: parameters use unified level-indexed structure.
// Parameters are sized by enabled group counts to avoid wasted sampling.

// Population intercept (always on)
real frac_logit_loc_pop;

// Population covariate coefficients (QR space) — length 0 if disabled
vector[enable_pop_cov_frac ? n_covar : 0] frac_coef_qr_pop;

// ===== UNIFIED LEVEL STRUCTURE =====

// Intercept SD free parameters - sized to RE levels only (mode=2)
// FE levels (mode=1) use frac_fe_sd_level_intercept from hyperparams; disabled levels use 0.
// Full n_levels array frac_sd_level_intercept is assembled in transformed_parameters.
array[n_re_levels_frac_intercept] real<lower=0> frac_sd_level_intercept_raw;

// Raw NCP and CP draws for intercepts — split by RE vs RE_CP mode
vector[n_raw_groups_frac_intercept] frac_raw_level_intercept;
vector[n_cp_groups_frac_intercept]  frac_cp_level_intercept;

// Slope SD hyperparameters - one vector per level (always n_levels for simplicity)
array[n_levels] vector<lower=0>[n_covar] frac_sd_level_slope;

// Raw NCP and CP draws for slopes — split by RE vs RE_CP mode
matrix[n_raw_groups_frac_slope, n_covar] frac_raw_level_slope;
matrix[n_cp_groups_frac_slope,  n_covar] frac_cp_level_slope;

// Student-t hierarchy: degrees of freedom per level (size 0 when disabled)
array[enable_student_t_hierarchy ? n_levels : 0] real<lower=2> frac_nu_level;
