// ============================================================================
// Other Events Model Parameters
// ============================================================================

// --- Baseline Hazard GP (Population Level) ---
row_vector[n_causes] log_lambda_gp_pop_intercept;
vector<lower=0>[n_causes] log_lambda_gp_pop_alpha;
vector<lower=0>[n_causes] log_lambda_gp_pop_rho;
array[n_causes] row_vector[max_all_t] log_lambda_gp_pop_eta;

// --- Baseline Hazard GP (Trial Level) ---
vector<lower=0>[oe_enable_trial_baseline_hazard ? n_trials : 0] log_lambda_gp_trial_alpha;
vector<lower=0>[oe_enable_trial_baseline_hazard ? n_trials : 0] log_lambda_gp_trial_rho;
array[n_causes] matrix[oe_enable_trial_baseline_hazard ? n_trials : 0, max_all_t] log_lambda_gp_trial_eta;
array[n_causes] vector[oe_enable_trial_baseline_hazard ? n_trials : 0] raw_log_lambda_gp_trial_intercept;
vector<lower=0>[oe_enable_trial_baseline_hazard ? n_causes : 0] log_lambda_gp_trial_intercept_sd;

// --- Proportional Hazard: Population-level Coefficients (QR space) ---
array[n_causes] vector[oe_enable_pop_tumor_cov ? n_tumor_covar : 0] oe_tumor_coef_qr_pop;
array[n_causes] vector[oe_enable_pop_cov ? n_covar : 0] oe_covar_coef_qr_pop;

// --- Proportional Hazard: Trial-level Random Slopes (non-centered) ---
array[n_causes] vector<lower=0>[oe_enable_trial_tumor_cov ? n_tumor_covar : 0] oe_sd_trial_tumor_slope;
array[n_causes] matrix[oe_enable_trial_tumor_cov ? n_trials : 0, oe_enable_trial_tumor_cov ? n_tumor_covar : 0] oe_raw_trial_tumor_slope;

array[n_causes] vector<lower=0>[oe_enable_trial_cov ? n_covar : 0] oe_sd_trial_slope;
array[n_causes] matrix[oe_enable_trial_cov ? n_trials : 0, oe_enable_trial_cov ? n_covar : 0] oe_raw_trial_slope;
