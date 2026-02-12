// ============================================================================
// Multistate Hazard Model Hyperparameters
// ============================================================================

// --- Baseline Hazard GP: Population-level ---
// 0→1 transition
real<lower=0> ms_log_lambda_gp_01_pop_alpha_alpha;
real<lower=0> ms_log_lambda_gp_01_pop_alpha_beta;
real<lower=0> ms_log_lambda_gp_01_pop_rho_alpha;
real<lower=0> ms_log_lambda_gp_01_pop_rho_beta;
real ms_log_lambda_gp_01_pop_intercept_mean;
real<lower=0> ms_log_lambda_gp_01_pop_intercept_sd;

// 0→2 transition
real<lower=0> ms_log_lambda_gp_02_pop_alpha_alpha;
real<lower=0> ms_log_lambda_gp_02_pop_alpha_beta;
real<lower=0> ms_log_lambda_gp_02_pop_rho_alpha;
real<lower=0> ms_log_lambda_gp_02_pop_rho_beta;
real ms_log_lambda_gp_02_pop_intercept_mean;
real<lower=0> ms_log_lambda_gp_02_pop_intercept_sd;

// 1→2 transition (sojourn time GP)
real<lower=0> ms_log_lambda_gp_12_s_pop_alpha_alpha;
real<lower=0> ms_log_lambda_gp_12_s_pop_alpha_beta;
real<lower=0> ms_log_lambda_gp_12_s_pop_rho_alpha;
real<lower=0> ms_log_lambda_gp_12_s_pop_rho_beta;
real ms_log_lambda_gp_12_s_pop_intercept_mean;
real<lower=0> ms_log_lambda_gp_12_s_pop_intercept_sd;

// 1→2 transition (clock-forward time GP, for Markov/extended)
real<lower=0> ms_log_lambda_gp_12_t_pop_alpha_alpha;
real<lower=0> ms_log_lambda_gp_12_t_pop_alpha_beta;
real<lower=0> ms_log_lambda_gp_12_t_pop_rho_alpha;
real<lower=0> ms_log_lambda_gp_12_t_pop_rho_beta;
real ms_log_lambda_gp_12_t_pop_intercept_mean;
real<lower=0> ms_log_lambda_gp_12_t_pop_intercept_sd;

// --- Baseline Hazard GP: Level-level (N-level hierarchy) ---
// SD for level intercepts (per level)
array[n_levels] real<lower=0> ms_log_lambda_gp_01_level_intercept_sd_sd;
array[n_levels] real<lower=0> ms_log_lambda_gp_02_level_intercept_sd_sd;
array[n_levels] real<lower=0> ms_log_lambda_gp_12_s_level_intercept_sd_sd;
array[n_levels] real<lower=0> ms_log_lambda_gp_12_t_level_intercept_sd_sd;

// Alpha and rho hyperparameters for level GPs
array[n_levels] real<lower=0> ms_log_lambda_gp_01_level_alpha_alpha;
array[n_levels] real<lower=0> ms_log_lambda_gp_01_level_alpha_beta;
array[n_levels] real<lower=0> ms_log_lambda_gp_01_level_rho_alpha;
array[n_levels] real<lower=0> ms_log_lambda_gp_01_level_rho_beta;

array[n_levels] real<lower=0> ms_log_lambda_gp_02_level_alpha_alpha;
array[n_levels] real<lower=0> ms_log_lambda_gp_02_level_alpha_beta;
array[n_levels] real<lower=0> ms_log_lambda_gp_02_level_rho_alpha;
array[n_levels] real<lower=0> ms_log_lambda_gp_02_level_rho_beta;

array[n_levels] real<lower=0> ms_log_lambda_gp_12_s_level_alpha_alpha;
array[n_levels] real<lower=0> ms_log_lambda_gp_12_s_level_alpha_beta;
array[n_levels] real<lower=0> ms_log_lambda_gp_12_s_level_rho_alpha;
array[n_levels] real<lower=0> ms_log_lambda_gp_12_s_level_rho_beta;

array[n_levels] real<lower=0> ms_log_lambda_gp_12_t_level_alpha_alpha;
array[n_levels] real<lower=0> ms_log_lambda_gp_12_t_level_alpha_beta;
array[n_levels] real<lower=0> ms_log_lambda_gp_12_t_level_rho_alpha;
array[n_levels] real<lower=0> ms_log_lambda_gp_12_t_level_rho_beta;

// --- Covariate Coefficient Hyperparameters ---
// Time-varying coefficients (population-level)
vector[n_time_varying_covar] ms_time_varying_coef_01_mean;
vector<lower=0>[n_time_varying_covar] ms_time_varying_coef_01_sd;
vector[n_time_varying_covar] ms_time_varying_coef_02_mean;
vector<lower=0>[n_time_varying_covar] ms_time_varying_coef_02_sd;
vector[n_time_varying_covar] ms_time_varying_coef_12_mean;
vector<lower=0>[n_time_varying_covar] ms_time_varying_coef_12_sd;

// Time-invariant coefficients (population-level, QR space)
vector[n_time_invariant_covar] ms_time_invariant_coef_01_mean;
vector<lower=0>[n_time_invariant_covar] ms_time_invariant_coef_01_sd;
vector[n_time_invariant_covar] ms_time_invariant_coef_02_mean;
vector<lower=0>[n_time_invariant_covar] ms_time_invariant_coef_02_sd;
vector[n_time_invariant_covar] ms_time_invariant_coef_12_mean;
vector<lower=0>[n_time_invariant_covar] ms_time_invariant_coef_12_sd;

// Multi-level random slope SDs (per level)
array[n_levels] row_vector<lower=0>[n_time_invariant_covar] ms_sd_level_slope_01_sd;
array[n_levels] row_vector<lower=0>[n_time_invariant_covar] ms_sd_level_slope_02_sd;
array[n_levels] row_vector<lower=0>[n_time_invariant_covar] ms_sd_level_slope_12_sd;
