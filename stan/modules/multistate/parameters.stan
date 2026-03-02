// ============================================================================
// Multistate Hazard Model Parameters
// ============================================================================
// Parameters are sized conditionally based on enabled transitions

// ============================================================================
// 0→1 TRANSITION PARAMETERS (only if enable_ms_01=1)
// ============================================================================

// --- Population-level Baseline Hazard GP ---
array[enable_ms_01 ? 1 : 0] real log_lambda_gp_01_pop_intercept;
array[enable_ms_01 ? 1 : 0] real<lower=0> log_lambda_gp_01_pop_alpha;
array[enable_ms_01 ? 1 : 0] real<lower=0> log_lambda_gp_01_pop_rho;
row_vector[enable_ms_01 ? max_all_t : 0] log_lambda_gp_01_pop_eta;

// --- Level-level Baseline Hazard GP (N-level hierarchy) ---
array[enable_ms_01 ? n_levels : 0] real<lower=0> log_lambda_gp_01_level_alpha;
array[enable_ms_01 ? n_levels : 0] real<lower=0> log_lambda_gp_01_level_rho;
array[enable_ms_01 ? n_levels : 0] real<lower=0> log_lambda_gp_01_level_intercept_sd;
matrix[n_enabled_groups_ms_baseline_01, enable_ms_01 ? max_all_t : 0] log_lambda_gp_01_level_eta;
vector[n_enabled_groups_ms_baseline_01] raw_log_lambda_gp_01_level_intercept;

// --- Population-level Time-varying Covariates ---
vector[enable_ms_01 && enable_ms_pop_time_varying_cov ? n_time_varying_covar : 0] time_varying_coef_01;

// --- Population-level Time-invariant Covariates (QR space) ---
vector[enable_ms_01 && enable_ms_pop_time_invariant_cov ? n_time_invariant_covar : 0] time_invariant_coef_qr_01;

// --- Multi-level Random Slopes for Time-invariant Covariates ---
array[n_levels] vector<lower=0>[enable_ms_01 ? n_time_invariant_covar : 0] sd_level_slope_01;
matrix[enable_ms_01 ? n_enabled_groups_ms_slope : 0, n_time_invariant_covar] raw_level_slope_01;

// ============================================================================
// 0→2 TRANSITION PARAMETERS (only if enable_ms_02=1)
// ============================================================================

// --- Population-level Baseline Hazard GP ---
array[enable_ms_02 ? 1 : 0] real log_lambda_gp_02_pop_intercept;
array[enable_ms_02 ? 1 : 0] real<lower=0> log_lambda_gp_02_pop_alpha;
array[enable_ms_02 ? 1 : 0] real<lower=0> log_lambda_gp_02_pop_rho;
row_vector[enable_ms_02 ? max_all_t : 0] log_lambda_gp_02_pop_eta;

// --- Level-level Baseline Hazard GP ---
// Conditional on transition being enabled to avoid improper posteriors
array[enable_ms_02 ? n_levels : 0] real<lower=0> log_lambda_gp_02_level_alpha;
array[enable_ms_02 ? n_levels : 0] real<lower=0> log_lambda_gp_02_level_rho;
array[enable_ms_02 ? n_levels : 0] real<lower=0> log_lambda_gp_02_level_intercept_sd;
matrix[n_enabled_groups_ms_baseline_02, enable_ms_02 ? max_all_t : 0] log_lambda_gp_02_level_eta;
vector[n_enabled_groups_ms_baseline_02] raw_log_lambda_gp_02_level_intercept;

// --- Population-level Time-varying Covariates ---
vector[enable_ms_02 && enable_ms_pop_time_varying_cov ? n_time_varying_covar : 0] time_varying_coef_02;

// --- Population-level Time-invariant Covariates (QR space) ---
vector[enable_ms_02 && enable_ms_pop_time_invariant_cov ? n_time_invariant_covar : 0] time_invariant_coef_qr_02;

// --- Multi-level Random Slopes ---
array[n_levels] vector<lower=0>[enable_ms_02 ? n_time_invariant_covar : 0] sd_level_slope_02;
matrix[enable_ms_02 ? n_enabled_groups_ms_slope : 0, n_time_invariant_covar] raw_level_slope_02;

// ============================================================================
// 1→2 TRANSITION PARAMETERS (only if enable_ms_12=1)
// ============================================================================

// --- Sojourn Time GP (semi-Markov or extended) ---
array[need_12_s_gp ? 1 : 0] real log_lambda_gp_12_s_pop_intercept;
array[need_12_s_gp ? 1 : 0] real<lower=0> log_lambda_gp_12_s_pop_alpha;
array[need_12_s_gp ? 1 : 0] real<lower=0> log_lambda_gp_12_s_pop_rho;
row_vector[need_12_s_gp ? ms_max_sojourn_t : 0] log_lambda_gp_12_s_pop_eta;

// Level hierarchy for sojourn GP (conditional on need_12_s_gp)
array[need_12_s_gp ? n_levels : 0] real<lower=0> log_lambda_gp_12_s_level_alpha;
array[need_12_s_gp ? n_levels : 0] real<lower=0> log_lambda_gp_12_s_level_rho;
array[need_12_s_gp ? n_levels : 0] real<lower=0> log_lambda_gp_12_s_level_intercept_sd;
matrix[n_enabled_groups_ms_baseline_12_s, need_12_s_gp ? ms_max_sojourn_t : 0] log_lambda_gp_12_s_level_eta;
vector[n_enabled_groups_ms_baseline_12_s] raw_log_lambda_gp_12_s_level_intercept;

// --- Clock-forward Time GP (Markov or extended) ---
array[need_12_t_gp ? 1 : 0] real log_lambda_gp_12_t_pop_intercept;
array[need_12_t_gp ? 1 : 0] real<lower=0> log_lambda_gp_12_t_pop_alpha;
array[need_12_t_gp ? 1 : 0] real<lower=0> log_lambda_gp_12_t_pop_rho;
row_vector[need_12_t_gp ? max_all_t : 0] log_lambda_gp_12_t_pop_eta;

// Level hierarchy for clock-forward GP (conditional on need_12_t_gp)
array[need_12_t_gp ? n_levels : 0] real<lower=0> log_lambda_gp_12_t_level_alpha;
array[need_12_t_gp ? n_levels : 0] real<lower=0> log_lambda_gp_12_t_level_rho;
array[need_12_t_gp ? n_levels : 0] real<lower=0> log_lambda_gp_12_t_level_intercept_sd;
matrix[n_enabled_groups_ms_baseline_12_t, need_12_t_gp ? max_all_t : 0] log_lambda_gp_12_t_level_eta;
vector[n_enabled_groups_ms_baseline_12_t] raw_log_lambda_gp_12_t_level_intercept;

// --- Population-level Covariates for 1→2 ---
vector[enable_ms_12 && enable_ms_pop_time_varying_cov ? n_time_varying_covar : 0] time_varying_coef_12;
vector[enable_ms_12 && enable_ms_pop_time_invariant_cov ? n_time_invariant_covar : 0] time_invariant_coef_qr_12;

// --- Multi-level Random Slopes for 1→2 ---
array[n_levels] vector<lower=0>[enable_ms_12 ? n_time_invariant_covar : 0] sd_level_slope_12;
matrix[enable_ms_12 ? n_enabled_groups_ms_slope : 0, n_time_invariant_covar] raw_level_slope_12;

// ============================================================================
// 0→3 TRANSITION PARAMETERS (GP baseline hazard, N-level hierarchy)
// ============================================================================
array[enable_ms_03 ? 1 : 0] real log_lambda_gp_03_pop_intercept;
array[enable_ms_03 ? 1 : 0] real<lower=0> log_lambda_gp_03_pop_alpha;
array[enable_ms_03 ? 1 : 0] real<lower=0> log_lambda_gp_03_pop_rho;
row_vector[enable_ms_03 ? max_all_t : 0] log_lambda_gp_03_pop_eta;
array[enable_ms_03 ? n_levels : 0] real<lower=0> log_lambda_gp_03_level_alpha;
array[enable_ms_03 ? n_levels : 0] real<lower=0> log_lambda_gp_03_level_rho;
array[enable_ms_03 ? n_levels : 0] real<lower=0> log_lambda_gp_03_level_intercept_sd;
matrix[n_enabled_groups_ms_baseline_03, enable_ms_03 ? max_all_t : 0] log_lambda_gp_03_level_eta;
vector[n_enabled_groups_ms_baseline_03] raw_log_lambda_gp_03_level_intercept;

// ============================================================================
// 3→2 TRANSITION PARAMETERS (Sojourn time GP baseline hazard, semi-Markov)
// ============================================================================

// --- Population-level Baseline Hazard GP ---
array[enable_ms_32 ? 1 : 0] real log_lambda_gp_32_s_pop_intercept;
array[enable_ms_32 ? 1 : 0] real<lower=0> log_lambda_gp_32_s_pop_alpha;
array[enable_ms_32 ? 1 : 0] real<lower=0> log_lambda_gp_32_s_pop_rho;
row_vector[enable_ms_32 ? ms_max_sojourn_t_32 : 0] log_lambda_gp_32_s_pop_eta;

// --- Level-level Baseline Hazard GP ---
array[enable_ms_32 ? n_levels : 0] real<lower=0> log_lambda_gp_32_s_level_alpha;
array[enable_ms_32 ? n_levels : 0] real<lower=0> log_lambda_gp_32_s_level_rho;
array[enable_ms_32 ? n_levels : 0] real<lower=0> log_lambda_gp_32_s_level_intercept_sd;
matrix[n_enabled_groups_ms_baseline_32, enable_ms_32 ? ms_max_sojourn_t_32 : 0] log_lambda_gp_32_s_level_eta;
vector[n_enabled_groups_ms_baseline_32] raw_log_lambda_gp_32_s_level_intercept;
