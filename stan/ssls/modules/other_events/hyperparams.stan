// ============================================================================
// Other Events Model Hyperparameters
// ============================================================================

// --- Baseline Hazard GP: Population-level ---
vector<lower=0>[n_causes] oe_log_lambda_gp_pop_alpha_alpha;
vector<lower=0>[n_causes] oe_log_lambda_gp_pop_alpha_beta;
vector<lower=0>[n_causes] oe_log_lambda_gp_pop_rho_alpha;
vector<lower=0>[n_causes] oe_log_lambda_gp_pop_rho_beta;
vector[n_causes] oe_log_lambda_gp_pop_intercept_mean;
vector<lower=0>[n_causes] oe_log_lambda_gp_pop_intercept_sd;

// --- Baseline Hazard GP: Trial-level ---
vector<lower=0>[n_causes] oe_log_lambda_gp_trial_alpha_alpha;
vector<lower=0>[n_causes] oe_log_lambda_gp_trial_alpha_beta;
vector<lower=0>[n_causes] oe_log_lambda_gp_trial_rho_alpha;
vector<lower=0>[n_causes] oe_log_lambda_gp_trial_rho_beta;
vector<lower=0>[n_causes] oe_log_lambda_gp_trial_intercept_sd_sd;

// --- Proportional Hazard: Population-level Covariate Coefficients ---
// Time-varying tumor coefficients (vector per cause, currently just log(SLD))
array[n_causes] vector[n_tumor_covar] oe_tumor_coef_pop_mean;
array[n_causes] vector<lower=0>[n_tumor_covar] oe_tumor_coef_pop_sd;

// Non-tumor covariates (QR space, vector per cause)
array[n_causes] vector[n_covar] oe_covar_coef_qr_pop_mean;
array[n_causes] vector<lower=0>[n_covar] oe_covar_coef_qr_pop_sd;

// --- Proportional Hazard: Multi-level Random Slope SDs ---
// Non-tumor covariate random slopes - one vector per cause per level
array[n_causes, n_levels] row_vector<lower=0>[n_covar] oe_sd_level_slope_sd;
