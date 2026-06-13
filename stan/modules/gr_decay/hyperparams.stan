// gr_decay/hyperparams.stan
// Hyperparameter (data) declarations for the Gompertz decay (gr_decay) module.
// Always declared with flag-independent shape, unused when off (cf. init/hyperparams.stan).

// Population intercept prior on log(kappa). Mean is NEGATIVE (e.g. log(0.02) ~ -3.9),
// so this is unconstrained real (NOT <lower=0>).
real gr_decay_log_loc_pop_mean;
real<lower=0> gr_decay_log_loc_pop_sd;

// QR-space covariate coefficient hyperparameters (applied in model block when enabled).
vector[n_covar] gr_decay_coef_qr_pop_mean;
vector<lower=0>[n_covar] gr_decay_coef_qr_pop_sd;

// Hierarchical intercept prior scale hyperparameters — one per level (RE/RE_CP).
array[n_levels] real<lower=0> gr_decay_sd_level_intercept_sd;

// Fixed-effect SD for FE levels (mode=1); unused for RE/disabled levels.
array[n_levels] real<lower=0> gr_decay_fe_sd_level_intercept;

// Hierarchical slope SD hyperpriors — one row_vector per level.
array[n_levels] row_vector<lower=0>[n_covar] gr_decay_sd_level_slope_sd;
