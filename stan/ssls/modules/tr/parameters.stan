// tr/parameters.stan
// Activated Total Rate (tr) module parameter declarations.
// For compatibility with legacy downstream code we expose pop_log_total_rate,
// trial_log_total_rate_effect, patient_log_total_rate_effect, patient_log_total_rate
// (the latter two are derived in transformed parameters but their raw components
// are declared here as raw standard normal draws + scales).
// Covariate slope components (QR) are scaffolded but currently gated off unless
// enable_*_cov_tr flags are set from data. We reuse existing n_covar design size.

// Population intercept (always on)
real pop_log_total_rate;

// Population covariate coefficients (QR space) — length 0 if disabled
vector[enable_pop_cov_tr ? n_covar : 0] tr_coef_qr_pop;

// Trial-level random intercept hierarchy
real<lower=0> tr_sd_trial_intercept; // prior scale hyperparam
vector[enable_trial_intercept_tr ? n_train_trials : 0] tr_raw_trial_intercept; // std normal draws

// Patient-level random intercept hierarchy
real<lower=0> tr_sd_patient_intercept;
vector[enable_patient_intercept_tr ? n_train_patients : 0] tr_raw_patient_intercept;

// Trial-level covariate slope deviations (originally on QR scale)
vector<lower=0>[enable_trial_cov_tr ? n_covar : 0] tr_sd_trial_slope;
matrix[enable_trial_cov_tr ? n_train_trials : 0, enable_trial_cov_tr ? n_covar : 0] tr_raw_trial_slope;

// Patient-level covariate slope deviations
vector<lower=0>[enable_patient_cov_tr ? n_covar : 0] tr_sd_patient_slope;
matrix[enable_patient_cov_tr ? n_train_patients : 0, enable_patient_cov_tr ? n_covar : 0] tr_raw_patient_slope;
