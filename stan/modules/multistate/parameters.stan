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
row_vector[enable_ms_01 ? n_ms_gp_cal_knots : 0] log_lambda_gp_01_pop_eta;

// --- Level-level Baseline Hazard GP (N-level hierarchy) ---
array[enable_ms_01 ? n_levels : 0] real<lower=0> log_lambda_gp_01_level_alpha;
array[enable_ms_01 ? n_levels : 0] real<lower=0> log_lambda_gp_01_level_rho;
array[enable_ms_01 && any_re_level ? n_levels : 0] real<lower=0> log_lambda_gp_01_level_intercept_sd;
matrix[n_gp_groups_ms_baseline_01, enable_ms_01 ? n_ms_gp_cal_knots : 0] log_lambda_gp_01_level_eta;
vector[n_raw_groups_ms_baseline_01] raw_log_lambda_gp_01_level_intercept;
vector[n_cp_groups_ms_baseline_01]  cp_log_lambda_gp_01_level_intercept;

// --- Population-level Time-varying Covariates ---
// 0->1 dimensioning rules:
//   continuous mode (visit_gated_01 = 0): n_time_varying_covar features (full
//     SLD + decrease + growth covariate set, evaluated continuously over weeks).
//   latent visit-gated mode (visit_gated_01 = 1, latent = 1): n_time_varying_covar
//     features (same modeled covariates, just evaluated at visit weeks).
//   observed visit-gated mode (visit_gated_01 = 1, latent = 0): single feature
//     (one observed biomarker via ms_obs_visit_covar_flat per visit).
vector[enable_ms_01 && enable_ms_pop_time_varying_cov
    ? (enable_ms_visit_gated_01 && !enable_ms_visit_gated_latent_01 ? 1 : n_time_varying_covar) : 0] time_varying_coef_01;

// --- Population-level Time-invariant Covariates (QR space) ---
vector[enable_ms_01 && enable_ms_pop_time_invariant_cov ? n_time_invariant_covar : 0] time_invariant_coef_qr_01;

// --- Multi-level Random Slopes for Time-invariant Covariates ---
array[n_levels] vector<lower=0>[enable_ms_01 ? n_time_invariant_covar : 0] sd_level_slope_01;
matrix[enable_ms_01 ? n_raw_groups_ms_slope_shared : 0, n_time_invariant_covar] raw_level_slope_01;
matrix[enable_ms_01 ? n_cp_groups_ms_slope_shared  : 0, n_time_invariant_covar] cp_level_slope_01;

// ============================================================================
// 0→2 TRANSITION PARAMETERS (only if enable_ms_02=1)
// ============================================================================

// --- Population-level Baseline Hazard GP ---
array[enable_ms_02 ? 1 : 0] real log_lambda_gp_02_pop_intercept;
array[enable_ms_02 && !share_dead_gp_shape ? 1 : 0] real<lower=0> log_lambda_gp_02_pop_alpha;
array[enable_ms_02 && !share_dead_gp_shape ? 1 : 0] real<lower=0> log_lambda_gp_02_pop_rho;
row_vector[enable_ms_02 && !share_dead_gp_shape ? n_ms_gp_cal_knots : 0] log_lambda_gp_02_pop_eta;

// --- Level-level Baseline Hazard GP ---
// Conditional on transition being enabled to avoid improper posteriors
array[enable_ms_02 ? n_levels : 0] real<lower=0> log_lambda_gp_02_level_alpha;
array[enable_ms_02 ? n_levels : 0] real<lower=0> log_lambda_gp_02_level_rho;
array[enable_ms_02 && any_re_level ? n_levels : 0] real<lower=0> log_lambda_gp_02_level_intercept_sd;
matrix[n_gp_groups_ms_baseline_02, enable_ms_02 ? n_ms_gp_cal_knots : 0] log_lambda_gp_02_level_eta;
vector[n_raw_groups_ms_baseline_02] raw_log_lambda_gp_02_level_intercept;
vector[n_cp_groups_ms_baseline_02]  cp_log_lambda_gp_02_level_intercept;

// --- Population-level Time-varying Covariates ---
// 0->2: n_time_varying_covar features if enabled, 0 otherwise
vector[enable_ms_02 && enable_ms_pop_time_varying_cov && enable_ms_02_time_varying_cov
    ? n_time_varying_covar : 0] time_varying_coef_02;

// --- Population-level Time-invariant Covariates (QR space) ---
vector[enable_ms_02 && enable_ms_pop_time_invariant_cov ? n_time_invariant_covar : 0] time_invariant_coef_qr_02;

// --- Multi-level Random Slopes ---
array[n_levels] vector<lower=0>[enable_ms_02 ? n_time_invariant_covar : 0] sd_level_slope_02;
matrix[enable_ms_02 ? n_raw_groups_ms_slope_shared : 0, n_time_invariant_covar] raw_level_slope_02;
matrix[enable_ms_02 ? n_cp_groups_ms_slope_shared  : 0, n_time_invariant_covar] cp_level_slope_02;

// ============================================================================
// 1→2 TRANSITION PARAMETERS (only if enable_ms_12=1)
// ============================================================================

// --- Sojourn Time GP (semi-Markov or extended) ---
array[need_12_s_gp ? 1 : 0] real log_lambda_gp_12_s_pop_intercept;
array[need_12_s_gp ? 1 : 0] real<lower=0> log_lambda_gp_12_s_pop_alpha;
array[need_12_s_gp ? 1 : 0] real<lower=0> log_lambda_gp_12_s_pop_rho;
row_vector[need_12_s_gp ? n_ms_gp_sojourn_knots : 0] log_lambda_gp_12_s_pop_eta;

// Level hierarchy for sojourn GP (conditional on need_12_s_gp)
array[need_12_s_gp ? n_levels : 0] real<lower=0> log_lambda_gp_12_s_level_alpha;
array[need_12_s_gp ? n_levels : 0] real<lower=0> log_lambda_gp_12_s_level_rho;
array[need_12_s_gp && any_re_level ? n_levels : 0] real<lower=0> log_lambda_gp_12_s_level_intercept_sd;
matrix[n_gp_groups_ms_baseline_12_s, need_12_s_gp ? n_ms_gp_sojourn_knots : 0] log_lambda_gp_12_s_level_eta;
vector[n_raw_groups_ms_baseline_12_s] raw_log_lambda_gp_12_s_level_intercept;
vector[n_cp_groups_ms_baseline_12_s]  cp_log_lambda_gp_12_s_level_intercept;

// ============================================================================
// SHARED "DEAD" GP SHAPE (0→2 + 1→2 clock-forward, when share_dead_gp_shape=1)
// ============================================================================
array[share_dead_gp_shape ? 1 : 0] real<lower=0> log_lambda_gp_dead_pop_alpha;
array[share_dead_gp_shape ? 1 : 0] real<lower=0> log_lambda_gp_dead_pop_rho;
row_vector[share_dead_gp_shape ? n_ms_gp_cal_knots : 0] log_lambda_gp_dead_pop_eta;

// --- Clock-forward Time GP (Markov or extended) ---
array[need_12_t_gp ? 1 : 0] real log_lambda_gp_12_t_pop_intercept;
array[need_12_t_gp && !share_dead_gp_shape ? 1 : 0] real<lower=0> log_lambda_gp_12_t_pop_alpha;
array[need_12_t_gp && !share_dead_gp_shape ? 1 : 0] real<lower=0> log_lambda_gp_12_t_pop_rho;
row_vector[need_12_t_gp && !share_dead_gp_shape ? n_ms_gp_cal_knots : 0] log_lambda_gp_12_t_pop_eta;

// Level hierarchy for clock-forward GP (conditional on need_12_t_gp)
array[need_12_t_gp ? n_levels : 0] real<lower=0> log_lambda_gp_12_t_level_alpha;
array[need_12_t_gp ? n_levels : 0] real<lower=0> log_lambda_gp_12_t_level_rho;
array[need_12_t_gp && any_re_level ? n_levels : 0] real<lower=0> log_lambda_gp_12_t_level_intercept_sd;
matrix[n_gp_groups_ms_baseline_12_t, need_12_t_gp ? n_ms_gp_cal_knots : 0] log_lambda_gp_12_t_level_eta;
vector[n_raw_groups_ms_baseline_12_t] raw_log_lambda_gp_12_t_level_intercept;
vector[n_cp_groups_ms_baseline_12_t]  cp_log_lambda_gp_12_t_level_intercept;

// --- Population-level Covariates for 1→2 ---
vector[enable_ms_12 && enable_ms_pop_time_varying_cov ? n_time_varying_covar : 0] time_varying_coef_12;
vector[enable_ms_12 && enable_ms_pop_time_invariant_cov ? n_time_invariant_covar : 0] time_invariant_coef_qr_12;

// --- Multi-level Random Slopes for 1→2 ---
array[n_levels] vector<lower=0>[enable_ms_12 ? n_time_invariant_covar : 0] sd_level_slope_12;
matrix[enable_ms_12 ? n_raw_groups_ms_slope_shared : 0, n_time_invariant_covar] raw_level_slope_12;
matrix[enable_ms_12 ? n_cp_groups_ms_slope_shared  : 0, n_time_invariant_covar] cp_level_slope_12;

// ============================================================================
// 0→3 TRANSITION PARAMETERS (GP baseline hazard, N-level hierarchy)
// ============================================================================
array[enable_ms_03 ? 1 : 0] real log_lambda_gp_03_pop_intercept;
array[enable_ms_03 ? 1 : 0] real<lower=0> log_lambda_gp_03_pop_alpha;
array[enable_ms_03 ? 1 : 0] real<lower=0> log_lambda_gp_03_pop_rho;
row_vector[enable_ms_03 ? n_ms_gp_cal_knots : 0] log_lambda_gp_03_pop_eta;
array[enable_ms_03 ? n_levels : 0] real<lower=0> log_lambda_gp_03_level_alpha;
array[enable_ms_03 ? n_levels : 0] real<lower=0> log_lambda_gp_03_level_rho;
array[enable_ms_03 && any_re_level ? n_levels : 0] real<lower=0> log_lambda_gp_03_level_intercept_sd;
matrix[n_gp_groups_ms_baseline_03, enable_ms_03 ? n_ms_gp_cal_knots : 0] log_lambda_gp_03_level_eta;
vector[n_raw_groups_ms_baseline_03] raw_log_lambda_gp_03_level_intercept;
vector[n_cp_groups_ms_baseline_03]  cp_log_lambda_gp_03_level_intercept;

// --- Population-level Time-varying Covariates for 0->3 (tumor bridge) ---
vector[enable_ms_03 && enable_ms_pop_time_varying_cov && enable_ms_03_time_varying_cov
    ? n_time_varying_covar : 0] time_varying_coef_03;

// --- Population-level Time-invariant Covariates for 0->3 (QR space) ---
vector[enable_ms_03 && enable_ms_pop_time_invariant_cov && enable_ms_03_time_invariant_cov
    ? n_time_invariant_covar : 0] time_invariant_coef_qr_03;

// --- Multi-level Random Slopes for 0->3 ---
array[n_levels] vector<lower=0>[enable_ms_03 && enable_ms_03_time_invariant_cov ? n_time_invariant_covar : 0] sd_level_slope_03;
matrix[enable_ms_03 && enable_ms_03_time_invariant_cov ? n_raw_groups_ms_slope_shared : 0, n_time_invariant_covar] raw_level_slope_03;
matrix[enable_ms_03 && enable_ms_03_time_invariant_cov ? n_cp_groups_ms_slope_shared  : 0, n_time_invariant_covar] cp_level_slope_03;

// ============================================================================
// STUDENT-T HIERARCHY: DEGREES OF FREEDOM (size 0 when disabled)
// ============================================================================
array[enable_student_t_hierarchy ? n_levels : 0] real<lower=2> ms_nu_baseline_level;
array[enable_student_t_hierarchy ? n_levels : 0] real<lower=2> ms_nu_slope_level;

// ============================================================================
// PSA-AT-ENTRY COVARIATE COEFFICIENTS
// ============================================================================
// Scalar coefficient for standardized log-PSA at state entry.
// Shifts the entire sojourn hazard up/down per patient based on PSA burden.
array[enable_ms_12 && enable_ms_12_entry_covar ? 1 : 0] real coef_log_entry_covar_12;
array[enable_ms_32 && enable_ms_32_entry_covar ? 1 : 0] real coef_log_entry_covar_32;

// ============================================================================
// 3→2 TRANSITION PARAMETERS (Sojourn time GP baseline hazard, semi-Markov)
// ============================================================================

// --- Population-level Baseline Hazard GP ---
array[enable_ms_32 ? 1 : 0] real log_lambda_gp_32_s_pop_intercept;
array[enable_ms_32 ? 1 : 0] real<lower=0> log_lambda_gp_32_s_pop_alpha;
array[enable_ms_32 ? 1 : 0] real<lower=0> log_lambda_gp_32_s_pop_rho;
row_vector[enable_ms_32 ? n_ms_gp_sojourn_32_knots : 0] log_lambda_gp_32_s_pop_eta;

// --- Level-level Baseline Hazard GP ---
array[enable_ms_32 ? n_levels : 0] real<lower=0> log_lambda_gp_32_s_level_alpha;
array[enable_ms_32 ? n_levels : 0] real<lower=0> log_lambda_gp_32_s_level_rho;
array[enable_ms_32 && any_re_level ? n_levels : 0] real<lower=0> log_lambda_gp_32_s_level_intercept_sd;
matrix[n_gp_groups_ms_baseline_32, enable_ms_32 ? n_ms_gp_sojourn_32_knots : 0] log_lambda_gp_32_s_level_eta;
vector[n_raw_groups_ms_baseline_32] raw_log_lambda_gp_32_s_level_intercept;
vector[n_cp_groups_ms_baseline_32]  cp_log_lambda_gp_32_s_level_intercept;

// --- Population-level Time-invariant Covariates for 3->2 (QR space) ---
vector[enable_ms_32 && enable_ms_pop_time_invariant_cov && enable_ms_32_time_invariant_cov
    ? n_time_invariant_covar : 0] time_invariant_coef_qr_32;

// --- Multi-level Random Slopes for 3->2 ---
array[n_levels] vector<lower=0>[enable_ms_32 && enable_ms_32_time_invariant_cov ? n_time_invariant_covar : 0] sd_level_slope_32;
matrix[enable_ms_32 && enable_ms_32_time_invariant_cov ? n_raw_groups_ms_slope_shared : 0, n_time_invariant_covar] raw_level_slope_32;
matrix[enable_ms_32 && enable_ms_32_time_invariant_cov ? n_cp_groups_ms_slope_shared  : 0, n_time_invariant_covar] cp_level_slope_32;
