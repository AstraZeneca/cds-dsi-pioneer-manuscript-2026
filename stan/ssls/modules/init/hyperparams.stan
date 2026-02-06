// init/hyperparams.stan
// Hyperparameter (data) declarations for Initial State Proportion (init) module.
// Multi-level hierarchy: hyperparameters are arrays indexed by level (1..n_levels)

real init_logit_loc_pop_mean;
real<lower=0> init_logit_loc_pop_sd;

// Hierarchical intercept prior scale hyperparameters - one per level
array[n_levels] real<lower=0> init_sd_level_intercept_sd;

// QR-space coefficient hyperparameters (applied in model block)
vector[n_covar] init_coef_qr_pop_mean;
vector<lower=0>[n_covar] init_coef_qr_pop_sd;

// Hierarchical slope SD hyperpriors - one vector per level
array[n_levels] row_vector<lower=0>[n_covar] init_sd_level_slope_sd;
