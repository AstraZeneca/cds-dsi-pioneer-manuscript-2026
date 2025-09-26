// frac/hyperparams.stan
// Hyperparameter (data) declarations for Fraction Mix (frac) module.
// Mapping old legacy names -> new:
//   pop_decrease_frac_logit_mean  -> frac_logit_loc_pop_mean
//   pop_decrease_frac_logit_sd    -> frac_logit_loc_pop_sd
//   trial_decrease_frac_logit_coef_sd_sd -> frac_sd_trial_slope_sd (per-covariate of slope SD)
//   patient_decrease_frac_logit_coef_sd_sd -> frac_sd_patient_slope_sd
//   patient_decrease_frac_logit_sd_sd -> frac_sd_patient_intercept_sd
//   (trial intercept SD prior previously borrowed total-rate naming; introduce explicit:) trial_decrease_frac_logit_sd_sd -> frac_sd_trial_intercept_sd (if exists later)
// Covariate coefficient hyperparams:
//   pop_decrease_frac_logit_coef_mean -> frac_coef_pop_mean
//   pop_decrease_frac_logit_coef_sd   -> frac_coef_pop_sd

// Population intercept (logit scale of decrease fraction)
real frac_logit_loc_pop_mean;
real<lower=0> frac_logit_loc_pop_sd;

// Intercept hierarchy SD prior scales
real<lower=0> frac_sd_trial_intercept_sd;   // may be repurposed from existing trial frac intercept hyperparam (define in R later if missing)
real<lower=0> frac_sd_patient_intercept_sd; // from patient_decrease_frac_logit_sd_sd

// Population covariate coefficient hyperparameters
vector[n_covar] frac_coef_pop_mean;
vector<lower=0>[n_covar] frac_coef_pop_sd;

// Per-covariate slope SD hyperpriors
row_vector<lower=0>[n_covar] frac_sd_trial_slope_sd;
row_vector<lower=0>[n_covar] frac_sd_patient_slope_sd;
