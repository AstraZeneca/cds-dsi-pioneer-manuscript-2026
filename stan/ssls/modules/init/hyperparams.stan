// init/hyperparams.stan
// Hyperparameter (data) declarations for Initial State Proportion (init) module.
// Mapping old legacy names -> new:
//   pop_decrease_prop_logis_mean      -> init_logit_loc_pop_mean
//   pop_decrease_prop_logis_sd        -> init_logit_loc_pop_sd
//   trial_decrease_prop_logis_sd_sd   -> init_sd_trial_intercept_sd
//   patient_decrease_prop_logis_sd_sd -> init_sd_patient_intercept_sd
//   pop_decrease_prop_logis_coef_mean -> init_coef_pop_mean
//   pop_decrease_prop_logis_coef_sd   -> init_coef_pop_sd
// (No existing per-covariate trial/patient slope SD hyperparams for init; predeclare for symmetry.)

real init_logit_loc_pop_mean;
real<lower=0> init_logit_loc_pop_sd;

real<lower=0> init_sd_trial_intercept_sd;     // trial intercept SD prior scale
real<lower=0> init_sd_patient_intercept_sd;   // patient intercept SD prior scale

// QR-space coefficient hyperparameters (applied in model block)
vector[n_covar] init_coef_qr_pop_mean;                   // mean for QR coefficients (typically 0)
vector<lower=0>[n_covar] init_coef_qr_pop_sd;            // sd for QR coefficients (typically 1)

row_vector<lower=0>[n_covar] init_sd_trial_slope_sd;    // placeholder (may be zero-length if disabled in R)
row_vector<lower=0>[n_covar] init_sd_patient_slope_sd;  // placeholder
