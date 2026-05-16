// frac/hyperparams.stan
// Hyperparameter (data) declarations for Fraction Mix (frac) module.
// Multi-level hierarchy: hyperparameters are arrays indexed by level (1..n_levels)

// Population intercept (logit scale of decrease fraction)
real frac_logit_loc_pop_mean;
real<lower=0> frac_logit_loc_pop_sd;

// Hierarchical intercept prior scale hyperparameters - one per level
array[n_levels] real<lower=0> frac_sd_level_intercept_sd;

// Fixed-effect SD hyperparameter for FE levels (mode=1); unused for RE/disabled levels
array[n_levels] real<lower=0> frac_fe_sd_level_intercept;

// Student-t hierarchy: nu prior hyperparameters (one per level)
array[n_levels] real<lower=0> frac_nu_level_prior_alpha;
array[n_levels] real<lower=0> frac_nu_level_prior_beta;

// QR-space coefficient hyperparameters (applied in model block)
vector[n_covar] frac_coef_qr_pop_mean;
vector<lower=0>[n_covar] frac_coef_qr_pop_sd;

// Hierarchical slope SD hyperpriors - one vector per level
array[n_levels] row_vector<lower=0>[n_covar] frac_sd_level_slope_sd;
