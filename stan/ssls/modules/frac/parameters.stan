// frac/parameters.stan
// Activated Fraction (frac) module parameter declarations.
// Covariate slope components (QR) are scaffolded but currently gated off unless
// enable_*_cov_frac flags are set from data. We reuse existing n_covar design size.

// Population intercept (always on)
real frac_logit_loc_pop; // population intercept (logit scale)

// Population covariate coefficients (QR space) — length 0 if disabled
vector[enable_pop_cov_frac ? n_covar : 0] frac_coef_qr_pop;

// Trial-level random intercept hierarchy
real<lower=0> frac_sd_trial_intercept; // prior scale hyperparam
vector[enable_trial_intercept_frac ? n_trials : 0] frac_raw_trial_intercept; // std normal draws

// Patient-level random intercept hierarchy
real<lower=0> frac_sd_patient_intercept;
vector[enable_patient_intercept_frac ? n_patients : 0] frac_raw_patient_intercept;

// Trial-level covariate slope deviations (originally on QR scale)
vector<lower=0>[enable_trial_cov_frac ? n_covar : 0] frac_sd_trial_slope;
matrix[enable_trial_cov_frac ? n_trials : 0, enable_trial_cov_frac ? n_covar : 0] frac_raw_trial_slope;

// Patient-level covariate slope deviations
vector<lower=0>[enable_patient_cov_frac ? n_covar : 0] frac_sd_patient_slope;
matrix[enable_patient_cov_frac ? n_patients : 0, enable_patient_cov_frac ? n_covar : 0] frac_raw_patient_slope;
