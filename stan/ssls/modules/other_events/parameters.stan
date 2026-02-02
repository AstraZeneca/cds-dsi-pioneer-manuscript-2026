// ============================================================================
// Other Events Model Parameters
// ============================================================================

// --- Baseline Hazard GP (Population Level) ---
row_vector[n_causes] log_lambda_gp_pop_intercept;
vector<lower=0>[n_causes] log_lambda_gp_pop_alpha;
vector<lower=0>[n_causes] log_lambda_gp_pop_rho;
array[n_causes] row_vector[max_all_t] log_lambda_gp_pop_eta;

// --- Baseline Hazard GP (Trial Level) ---
vector<lower=0>[oe_enable_trial_baseline_hazard ? n_causes : 0] log_lambda_gp_trial_alpha;
vector<lower=0>[oe_enable_trial_baseline_hazard ? n_causes : 0] log_lambda_gp_trial_rho;
array[n_causes] matrix[oe_enable_trial_baseline_hazard ? n_trials : 0, max_all_t] log_lambda_gp_trial_eta;
array[n_causes] vector[oe_enable_trial_baseline_hazard ? n_trials : 0] raw_log_lambda_gp_trial_intercept;
vector<lower=0>[oe_enable_trial_baseline_hazard ? n_causes : 0] log_lambda_gp_trial_intercept_sd;

// --- Proportional Hazard: Population-level Coefficients ---
// Tumor-derived coefficients (size determined by hyperparameters):
//   Minimum 1: log(SLD) effect
//   If size >= 3: also includes log(decrease rate) and log(growth rate) effects
array[n_causes] vector[oe_enable_pop_tumor_cov ? n_tumor_covar : 0] oe_tumor_coef_pop;

// Non-tumor covariate coefficients (QR space)
array[n_causes] vector[oe_enable_pop_cov ? n_covar : 0] oe_covar_coef_qr_pop;

// --- Proportional Hazard: Multi-level Random Slopes (non-centered) ---
// Non-tumor covariate random slopes - unified level structure
// SD hyperparameters: one vector per cause per level
array[n_causes, n_levels] vector<lower=0>[n_covar] oe_sd_level_slope;
// Raw effects: flattened across all levels (n_total_groups rows)
array[n_causes] matrix[n_total_groups, n_covar] oe_raw_level_slope;
