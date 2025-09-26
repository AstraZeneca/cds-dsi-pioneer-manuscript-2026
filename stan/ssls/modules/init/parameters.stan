// init/parameters.stan — active Initial proportion module
real pop_decrease_prop_logis; // population intercept (legacy name preserved)
vector[enable_pop_cov_init ? n_covar : 0] init_coef_qr_pop;

// Trial intercept hierarchy
real<lower=0> init_sd_trial_intercept;
vector[enable_trial_intercept_init ? n_train_trials : 0] init_raw_trial_intercept;

// Patient intercept hierarchy
real<lower=0> init_sd_patient_intercept;
vector[enable_patient_intercept_init ? n_train_patients : 0] init_raw_patient_intercept;

// Trial slopes
vector<lower=0>[enable_trial_cov_init ? n_covar : 0] init_sd_trial_slope;
matrix[enable_trial_cov_init ? n_train_trials : 0, enable_trial_cov_init ? n_covar : 0] init_raw_trial_slope;

// Patient slopes
vector<lower=0>[enable_patient_cov_init ? n_covar : 0] init_sd_patient_slope;
matrix[enable_patient_cov_init ? n_train_patients : 0, enable_patient_cov_init ? n_covar : 0] init_raw_patient_slope;
