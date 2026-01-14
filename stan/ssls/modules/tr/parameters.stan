// tr/parameters.stan
// Activated Total Rate (tr) module parameter declarations.
// Covariate slope components (QR) are scaffolded but currently gated off unless
// enable_*_cov_tr flags are set from data. We reuse existing n_covar design size.

// Population intercept (always on)
real tr_loc_pop;

// Population covariate coefficients (QR space) — length 0 if disabled
vector[enable_pop_cov_tr ? n_covar : 0] tr_coef_qr_pop;

// Trial-level random intercept hierarchy
real<lower=0> tr_sd_trial_intercept; // prior scale hyperparam
vector[enable_trial_intercept_tr ? n_trials : 0] tr_raw_trial_intercept; // std normal draws

// Patient-level random intercept hierarchy
real<lower=0> tr_sd_patient_intercept;
vector[enable_patient_intercept_tr ? n_patients : 0] tr_raw_patient_intercept;

// Trial-level covariate slope deviations (originally on QR scale)
vector<lower=0>[enable_trial_cov_tr ? n_covar : 0] tr_sd_trial_slope;
matrix[enable_trial_cov_tr ? n_trials : 0, enable_trial_cov_tr ? n_covar : 0] tr_raw_trial_slope;

// Patient-level covariate slope deviations
vector<lower=0>[enable_patient_cov_tr ? n_covar : 0] tr_sd_patient_slope;
matrix[enable_patient_cov_tr ? n_patients : 0, enable_patient_cov_tr ? n_covar : 0] tr_raw_patient_slope;

// Patient-level process noise parameters - only declared when feature is enabled
matrix[enable_patient_process_noise_tr ? n_patients : 0, max_t_width] tr_raw_patient_process_noise;
array[enable_patient_process_noise_tr ? 1 : 0] real tr_log_sd_pop_process_noise;
array[enable_patient_process_noise_tr ? 1 : 0] real<lower=0> tr_sd_patient_log_sd_process_noise;
vector[enable_patient_process_noise_sd_tr ? n_patients : 0] tr_raw_patient_log_sd_process_noise;
array[enable_patient_process_noise_tr ? 1 : 0] real tr_logit_phi_pop_process_noise;
array[enable_patient_process_noise_tr ? 1 : 0] real<lower=0> tr_sd_patient_phi_process_noise;
vector[enable_patient_process_noise_phi_tr ? n_patients : 0] tr_raw_patient_phi_process_noise;

// Population-level time-varying process noise parameters (shared AR(1) across all patients)
array[enable_pop_process_noise_tr ? 1 : 0] real tr_log_sd_pop_process_noise_pop;
array[enable_pop_process_noise_tr ? 1 : 0] real tr_logit_phi_pop_process_noise_pop;
row_vector[enable_pop_process_noise_tr ? max_t_width : 0] tr_raw_pop_process_noise;