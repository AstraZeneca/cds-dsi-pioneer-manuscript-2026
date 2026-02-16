// tr/parameters.stan
// Activated Total Rate (tr) module parameter declarations.
// Multi-level hierarchy: parameters use unified level-indexed structure.
// Parameters are sized by enabled group counts to avoid wasted sampling.

// Population intercept (always on)
real tr_loc_pop;

// Population covariate coefficients (QR space) — length 0 if disabled
vector[enable_pop_cov_tr ? n_covar : 0] tr_coef_qr_pop;

// ===== UNIFIED LEVEL STRUCTURE =====

// Intercept SD hyperparameters - one per level (always n_levels for simplicity)
array[n_levels] real<lower=0> tr_sd_level_intercept;

// Raw standard normal draws for intercepts - sized by ENABLED groups only
// Size: n_enabled_groups_tr_intercept (sum of groups at enabled levels)
vector[n_enabled_groups_tr_intercept] tr_raw_level_intercept;

// Slope SD hyperparameters - one vector per level (always n_levels for simplicity)
array[n_levels] vector<lower=0>[n_covar] tr_sd_level_slope;

// Raw standard normal draws for slopes - sized by ENABLED groups only
matrix[n_enabled_groups_tr_slope, n_covar] tr_raw_level_slope;

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