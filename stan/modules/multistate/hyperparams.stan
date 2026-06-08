// ============================================================================
// Multistate Hazard Model Hyperparameters
// ============================================================================

// --- Baseline Hazard GP: Population-level ---
// 0→1 transition
real<lower=0> log_lambda_gp_01_pop_alpha_alpha;
real<lower=0> log_lambda_gp_01_pop_alpha_beta;
real<lower=0> log_lambda_gp_01_pop_rho_alpha;
real<lower=0> log_lambda_gp_01_pop_rho_beta;
real log_lambda_gp_01_pop_intercept_mean;
real<lower=0> log_lambda_gp_01_pop_intercept_sd;

// 0→2 transition
real<lower=0> log_lambda_gp_02_pop_alpha_alpha;
real<lower=0> log_lambda_gp_02_pop_alpha_beta;
real<lower=0> log_lambda_gp_02_pop_rho_alpha;
real<lower=0> log_lambda_gp_02_pop_rho_beta;
real log_lambda_gp_02_pop_intercept_mean;
real<lower=0> log_lambda_gp_02_pop_intercept_sd;

// 1→2 transition (sojourn time GP)
real<lower=0> log_lambda_gp_12_s_pop_alpha_alpha;
real<lower=0> log_lambda_gp_12_s_pop_alpha_beta;
real<lower=0> log_lambda_gp_12_s_pop_rho_alpha;
real<lower=0> log_lambda_gp_12_s_pop_rho_beta;
real log_lambda_gp_12_s_pop_intercept_mean;
real<lower=0> log_lambda_gp_12_s_pop_intercept_sd;

// 1→2 transition (clock-forward time GP, for Markov/extended)
real<lower=0> log_lambda_gp_12_t_pop_alpha_alpha;
real<lower=0> log_lambda_gp_12_t_pop_alpha_beta;
real<lower=0> log_lambda_gp_12_t_pop_rho_alpha;
real<lower=0> log_lambda_gp_12_t_pop_rho_beta;
real log_lambda_gp_12_t_pop_intercept_mean;
real<lower=0> log_lambda_gp_12_t_pop_intercept_sd;

// Shared dead GP shape hyperparameters (used when share_dead_gp_shape=1)
real<lower=0> log_lambda_gp_dead_pop_alpha_alpha;
real<lower=0> log_lambda_gp_dead_pop_alpha_beta;
real<lower=0> log_lambda_gp_dead_pop_rho_alpha;
real<lower=0> log_lambda_gp_dead_pop_rho_beta;

// --- Baseline Hazard GP: Level-level (N-level hierarchy) ---
// SD for level intercepts (per level)
array[n_levels] real<lower=0> log_lambda_gp_01_level_intercept_sd_sd;
array[n_levels] real<lower=0> log_lambda_gp_02_level_intercept_sd_sd;
array[n_levels] real<lower=0> log_lambda_gp_12_s_level_intercept_sd_sd;
array[n_levels] real<lower=0> log_lambda_gp_12_t_level_intercept_sd_sd;

// Alpha and rho hyperparameters for level GPs
array[n_levels] real<lower=0> log_lambda_gp_01_level_alpha_alpha;
array[n_levels] real<lower=0> log_lambda_gp_01_level_alpha_beta;
array[n_levels] real<lower=0> log_lambda_gp_01_level_rho_alpha;
array[n_levels] real<lower=0> log_lambda_gp_01_level_rho_beta;

array[n_levels] real<lower=0> log_lambda_gp_02_level_alpha_alpha;
array[n_levels] real<lower=0> log_lambda_gp_02_level_alpha_beta;
array[n_levels] real<lower=0> log_lambda_gp_02_level_rho_alpha;
array[n_levels] real<lower=0> log_lambda_gp_02_level_rho_beta;

array[n_levels] real<lower=0> log_lambda_gp_12_s_level_alpha_alpha;
array[n_levels] real<lower=0> log_lambda_gp_12_s_level_alpha_beta;
array[n_levels] real<lower=0> log_lambda_gp_12_s_level_rho_alpha;
array[n_levels] real<lower=0> log_lambda_gp_12_s_level_rho_beta;

array[n_levels] real<lower=0> log_lambda_gp_12_t_level_alpha_alpha;
array[n_levels] real<lower=0> log_lambda_gp_12_t_level_alpha_beta;
array[n_levels] real<lower=0> log_lambda_gp_12_t_level_rho_alpha;
array[n_levels] real<lower=0> log_lambda_gp_12_t_level_rho_beta;

// --- 0→3 Dropout GP: Population + N-level hierarchy ---
real log_lambda_gp_03_pop_intercept_mean;
real<lower=0> log_lambda_gp_03_pop_intercept_sd;
real<lower=0> log_lambda_gp_03_pop_alpha_alpha;
real<lower=0> log_lambda_gp_03_pop_alpha_beta;
real<lower=0> log_lambda_gp_03_pop_rho_alpha;
real<lower=0> log_lambda_gp_03_pop_rho_beta;
array[n_levels] real<lower=0> log_lambda_gp_03_level_intercept_sd_sd;
array[n_levels] real<lower=0> log_lambda_gp_03_level_alpha_alpha;
array[n_levels] real<lower=0> log_lambda_gp_03_level_alpha_beta;
array[n_levels] real<lower=0> log_lambda_gp_03_level_rho_alpha;
array[n_levels] real<lower=0> log_lambda_gp_03_level_rho_beta;

// --- 3→2 Baseline Hazard GP: Population-level (sojourn time, semi-Markov) ---
real<lower=0> log_lambda_gp_32_s_pop_alpha_alpha;
real<lower=0> log_lambda_gp_32_s_pop_alpha_beta;
real<lower=0> log_lambda_gp_32_s_pop_rho_alpha;
real<lower=0> log_lambda_gp_32_s_pop_rho_beta;
real log_lambda_gp_32_s_pop_intercept_mean;
real<lower=0> log_lambda_gp_32_s_pop_intercept_sd;

// --- 3→2 Baseline Hazard GP: Level-level ---
array[n_levels] real<lower=0> log_lambda_gp_32_s_level_intercept_sd_sd;
array[n_levels] real<lower=0> log_lambda_gp_32_s_level_alpha_alpha;
array[n_levels] real<lower=0> log_lambda_gp_32_s_level_alpha_beta;
array[n_levels] real<lower=0> log_lambda_gp_32_s_level_rho_alpha;
array[n_levels] real<lower=0> log_lambda_gp_32_s_level_rho_beta;

// --- Fixed-Effect Prior SD for Level Intercepts (used when flag == 1) ---
array[n_levels] real<lower=0> fe_log_lambda_gp_01_level_intercept_sd;
array[n_levels] real<lower=0> fe_log_lambda_gp_02_level_intercept_sd;
array[n_levels] real<lower=0> fe_log_lambda_gp_12_s_level_intercept_sd;
array[n_levels] real<lower=0> fe_log_lambda_gp_12_t_level_intercept_sd;
array[n_levels] real<lower=0> fe_log_lambda_gp_03_level_intercept_sd;
array[n_levels] real<lower=0> fe_log_lambda_gp_32_s_level_intercept_sd;

// --- Covariate Coefficient Hyperparameters ---
// Time-varying coefficients (population-level)
// 0->1 dimensioning matches parameters.stan:
//   single coefficient only in observed visit-gated mode; full
//   n_time_varying_covar in continuous mode and latent visit-gated mode.
vector[enable_ms_pop_time_varying_cov
    ? (enable_ms_visit_gated_01 && !enable_ms_visit_gated_latent_01 ? 1 : n_time_varying_covar) : 0] time_varying_coef_01_mean;
vector<lower=0>[enable_ms_pop_time_varying_cov
    ? (enable_ms_visit_gated_01 && !enable_ms_visit_gated_latent_01 ? 1 : n_time_varying_covar) : 0] time_varying_coef_01_sd;
vector[enable_ms_pop_time_varying_cov && enable_ms_02_time_varying_cov
    ? n_time_varying_covar : 0] time_varying_coef_02_mean;
vector<lower=0>[enable_ms_pop_time_varying_cov && enable_ms_02_time_varying_cov
    ? n_time_varying_covar : 0] time_varying_coef_02_sd;
vector[enable_ms_pop_time_varying_cov && enable_ms_03_time_varying_cov
    ? n_time_varying_covar : 0] time_varying_coef_03_mean;
vector<lower=0>[enable_ms_pop_time_varying_cov && enable_ms_03_time_varying_cov
    ? n_time_varying_covar : 0] time_varying_coef_03_sd;
vector[n_time_varying_covar] time_varying_coef_12_mean;
vector<lower=0>[n_time_varying_covar] time_varying_coef_12_sd;

// Time-invariant coefficients (population-level, QR space)
// 0->3 / 3->2 dimensioning matches parameters.stan: gated by per-transition
// flag so length-0 input from priors.R (when flag off) round-trips cleanly.
// 0->1 / 0->2 / 1->2 do not have per-transition TI flags, so always full length.
vector[n_time_invariant_covar] time_invariant_coef_01_mean;
vector<lower=0>[n_time_invariant_covar] time_invariant_coef_01_sd;
vector[n_time_invariant_covar] time_invariant_coef_02_mean;
vector<lower=0>[n_time_invariant_covar] time_invariant_coef_02_sd;
vector[enable_ms_03 && enable_ms_pop_time_invariant_cov && enable_ms_03_time_invariant_cov
    ? n_time_invariant_covar : 0] time_invariant_coef_03_mean;
vector<lower=0>[enable_ms_03 && enable_ms_pop_time_invariant_cov && enable_ms_03_time_invariant_cov
    ? n_time_invariant_covar : 0] time_invariant_coef_03_sd;
vector[n_time_invariant_covar] time_invariant_coef_12_mean;
vector<lower=0>[n_time_invariant_covar] time_invariant_coef_12_sd;
vector[enable_ms_32 && enable_ms_pop_time_invariant_cov && enable_ms_32_time_invariant_cov
    ? n_time_invariant_covar : 0] time_invariant_coef_32_mean;
vector<lower=0>[enable_ms_32 && enable_ms_pop_time_invariant_cov && enable_ms_32_time_invariant_cov
    ? n_time_invariant_covar : 0] time_invariant_coef_32_sd;

// Multi-level random slope SDs (per level)
array[n_levels] row_vector<lower=0>[n_time_invariant_covar] sd_level_slope_01_sd;
array[n_levels] row_vector<lower=0>[n_time_invariant_covar] sd_level_slope_02_sd;
array[n_levels] row_vector<lower=0>[n_time_invariant_covar] sd_level_slope_03_sd;
array[n_levels] row_vector<lower=0>[n_time_invariant_covar] sd_level_slope_12_sd;
array[n_levels] row_vector<lower=0>[n_time_invariant_covar] sd_level_slope_32_sd;

// --- Student-t hierarchy: nu prior hyperparameters for multistate ---
// Baseline intercepts (one per level)
array[n_levels] real<lower=0> ms_nu_baseline_level_prior_alpha;
array[n_levels] real<lower=0> ms_nu_baseline_level_prior_beta;
// Slopes (one per level)
array[n_levels] real<lower=0> ms_nu_slope_level_prior_alpha;
array[n_levels] real<lower=0> ms_nu_slope_level_prior_beta;

// --- Burden-at-State-Entry Coefficient Hyperparameters ---
// Prior for the scalar log-burden coefficient on 1→2 and 3→2 sojourn hazards.
// (PSA in pioneer; symbol generic for any future per-trial burden marker.)
real coef_log_entry_covar_12_mean;
real<lower=0> coef_log_entry_covar_12_sd;
real coef_log_entry_covar_32_mean;
real<lower=0> coef_log_entry_covar_32_sd;

// --- 0->1 Baseline Log-Time Trend Hyperparameters ---
// Priors for the hierarchical log-time trend on the 0->1 baseline hazard,
// mirroring the baseline GP hierarchy (population slope + per-level slope
// deviation). Only consumed when enable_ms_baseline_trend_01 == 1.
real log_lambda_trend_01_pop_mean;          // prior mean, population slope (e.g. 0)
real<lower=0> log_lambda_trend_01_pop_sd;   // prior sd,  population slope (e.g. 0.5)
array[n_levels] real<lower=0> log_lambda_trend_01_level_sd_sd;  // half-normal scale for per-level slope SD (e.g. 0.25)

// --- Correlated Intercept Block Hyperparameter (Phase 2) ---
// LKJ shape for the cross-transition frailty correlation Cholesky factors.
// eta = 2 (default): gently concentrates toward the identity, symmetric about
// zero. Single scalar shared by every correlated block. Unused (but always in
// scope) when n_ms_corr_blocks == 0.
real<lower=0> ms_intercept_corr_eta;
