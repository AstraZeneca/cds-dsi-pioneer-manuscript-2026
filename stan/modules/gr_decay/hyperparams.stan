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
