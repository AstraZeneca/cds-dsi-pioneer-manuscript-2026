// MIGRATED (total rate) -> see modules/tr/parameters.stan
// (Removed fraction & init parameters migrated to modules/frac and modules/init)
// (Retained process, GP, lag, noise related legacy parameters below)
real pop_log_growth_lag;
real pop_log_growth_transition_rate;

real<lower = 0> patient_log_growth_lag_sd;
vector[pop_growth_lag_param_only ? 0 : n_patients] raw_patient_log_growth_lag;

// Patient-level GP
// row_vector<lower = 0>[2] pop_tumor_gp_alpha;
real log_pop_tumor_gp_rho;

// Hierarchical length-scale parameter
// vector<lower = 0>[2] log_trial_tumor_gp_rho_sd;
// vector[add_trial_level_tumor_gp_param ? n_trials : 0] raw_log_trial_tumor_gp_rho_effect; 

real<lower = 0> log_patient_tumor_gp_rho_sd;
vector[independ_long_process_noise || pop_rho_param_only ? 0 : n_patients] raw_log_patient_tumor_gp_rho_effect; 

matrix[n_total_train_visits_m1, 2] raw_patient_process_noise;

// real<lower=0> pop_decrease_process_sd;
// real<lower=0> pop_growth_process_sd;
vector<lower = 0>[2] pop_process_sd;
cholesky_factor_corr[independ_cross_process_noise ? 0 : 2] L_process_corr;

real<lower=0> measure_sd;        // Measurement noise variance
