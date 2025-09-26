// frac/parameters.stan — active Fraction (decrease share) module
real pop_decrease_frac_logit; // population intercept (legacy name preserved for now)
vector[enable_pop_cov_frac ? n_covar : 0] frac_coef_qr_pop;

// Trial intercept hierarchy
real<lower=0> frac_sd_trial_intercept;
vector[enable_trial_intercept_frac ? n_train_trials : 0] frac_raw_trial_intercept;

// Patient intercept hierarchy
real<lower=0> frac_sd_patient_intercept;
vector[enable_patient_intercept_frac ? n_train_patients : 0] frac_raw_patient_intercept;

// Trial slopes
vector<lower=0>[enable_trial_cov_frac ? n_covar : 0] frac_sd_trial_slope; 
matrix[enable_trial_cov_frac ? n_train_trials : 0, enable_trial_cov_frac ? n_covar : 0] frac_raw_trial_slope;

// Patient slopes
vector<lower=0>[enable_patient_cov_frac ? n_covar : 0] frac_sd_patient_slope;
matrix[enable_patient_cov_frac ? n_train_patients : 0, enable_patient_cov_frac ? n_covar : 0] frac_raw_patient_slope;
