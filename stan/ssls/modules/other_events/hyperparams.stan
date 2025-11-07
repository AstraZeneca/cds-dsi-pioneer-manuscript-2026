// ============================================================================
// Other Events Model Hyperparameters
// ============================================================================

int<lower=1> n_causes; // Number of competing risk causes (e.g., death, non-target PD)

// --- Baseline Hazard GP: Population-level ---
vector<lower=0>[n_causes] oe_log_lambda_gp_pop_alpha_sd;
vector<lower=0>[n_causes] oe_log_lambda_gp_pop_rho_alpha;
vector<lower=0>[n_causes] oe_log_lambda_gp_pop_rho_beta;
vector[n_causes] oe_log_lambda_gp_pop_intercept_mean;
vector<lower=0>[n_causes] oe_log_lambda_gp_pop_intercept_sd;

// --- Baseline Hazard GP: Trial-level ---
vector<lower=0>[n_causes] oe_log_lambda_gp_trial_alpha_sd;
vector<lower=0>[n_causes] oe_log_lambda_gp_trial_rho_alpha;
vector<lower=0>[n_causes] oe_log_lambda_gp_trial_rho_beta;
vector<lower=0>[n_causes] oe_log_lambda_gp_trial_intercept_sd_sd;

// --- Proportional Hazard: Population-level Covariate Coefficients (QR space) ---
array[n_causes] vector[n_tumor_covar] oe_tumor_coef_qr_pop_mean;
array[n_causes] vector<lower=0>[n_tumor_covar] oe_tumor_coef_qr_pop_sd;
array[n_causes] vector[n_covar] oe_covar_coef_qr_pop_mean;
array[n_causes] vector<lower=0>[n_covar] oe_covar_coef_qr_pop_sd;

// --- Proportional Hazard: Trial-level Random Slope SDs ---
array[n_causes] row_vector<lower=0>[n_tumor_covar] oe_sd_trial_tumor_slope_sd;
array[n_causes] row_vector<lower=0>[n_covar] oe_sd_trial_slope_sd;
