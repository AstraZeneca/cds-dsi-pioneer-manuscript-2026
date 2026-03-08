# nolint start: object_usage_linter

# Shared multistate initialization helper (prior-sampling version)
#
# Called by both the SCLC initializer (create_tumor_ssls_initializer) and
# the pioneer initializer (create_pioneer_initializer).
# Must be called inside a with(stan_data, { ... }) block so that all
# multistate fields (flags, hyperparameters, dimensions) are in scope.
#
# Compare with ms_init_values_fixed() in initializers_fixed.R which uses
# deterministic starting values instead of prior draws.

ms_init_values <- function(env) {
  with(env, {
    # Derived flags for 1→2 GPs
    need_12_s_gp <- enable_ms_12 && (ms_time_scale_12 == 1 || ms_time_scale_12 == 2)
    need_12_t_gp <- enable_ms_12 && (ms_time_scale_12 == 0 || ms_time_scale_12 == 2)

    # Multistate enabled group counts - SEPARATE for each transition
    n_enabled_groups_ms_baseline_01 <- if (enable_ms_01) sum(n_groups_per_level[enable_ms_level_baseline_hazard == 1]) else 0L
    n_enabled_groups_ms_baseline_02 <- if (enable_ms_02) sum(n_groups_per_level[enable_ms_level_baseline_hazard == 1]) else 0L
    n_enabled_groups_ms_baseline_12_s <- if (need_12_s_gp) sum(n_groups_per_level[enable_ms_level_baseline_hazard == 1]) else 0L
    n_enabled_groups_ms_baseline_12_t <- if (need_12_t_gp) sum(n_groups_per_level[enable_ms_level_baseline_hazard == 1]) else 0L
    n_enabled_groups_ms_slope <- sum(n_groups_per_level[enable_ms_level_cov == 1])

    tibble::lst(
      # --- 0→1 Transition (Progression / PFS event) ---
      log_lambda_gp_01_pop_intercept = if (enable_ms_01) {
        array(rnorm(1, log_lambda_gp_01_pop_intercept_mean, log_lambda_gp_01_pop_intercept_sd), dim = 1)
      },
      log_lambda_gp_01_pop_alpha = if (enable_ms_01) array(1.0, dim = 1),
      log_lambda_gp_01_pop_rho = if (enable_ms_01) {
        array(invgamma::rinvgamma(1, log_lambda_gp_01_pop_rho_alpha, log_lambda_gp_01_pop_rho_beta), dim = 1)
      },
      log_lambda_gp_01_pop_eta = if (enable_ms_01) rnorm(max_all_t),
      log_lambda_gp_01_level_alpha = rep(1.0, n_levels),
      log_lambda_gp_01_level_rho = invgamma::rinvgamma(n_levels, log_lambda_gp_01_level_rho_alpha, log_lambda_gp_01_level_rho_beta),
      log_lambda_gp_01_level_intercept_sd = abs(rnorm(n_levels, sd = log_lambda_gp_01_level_intercept_sd_sd)),
      log_lambda_gp_01_level_eta = if (enable_ms_01 && n_enabled_groups_ms_baseline_01 > 0) {
        matrix(rnorm(n_enabled_groups_ms_baseline_01 * max_all_t), nrow = n_enabled_groups_ms_baseline_01, ncol = max_all_t)
      },
      raw_log_lambda_gp_01_level_intercept = if (n_enabled_groups_ms_baseline_01 > 0) {
        rnorm(n_enabled_groups_ms_baseline_01)
      },

      # --- 0→2 Transition (Death without progression) ---
      log_lambda_gp_02_pop_intercept = if (enable_ms_02) {
        array(rnorm(1, log_lambda_gp_02_pop_intercept_mean, log_lambda_gp_02_pop_intercept_sd), dim = 1)
      },
      log_lambda_gp_02_pop_alpha = if (enable_ms_02) array(1.0, dim = 1),
      log_lambda_gp_02_pop_rho = if (enable_ms_02) {
        array(invgamma::rinvgamma(1, log_lambda_gp_02_pop_rho_alpha, log_lambda_gp_02_pop_rho_beta), dim = 1)
      },
      log_lambda_gp_02_pop_eta = if (enable_ms_02) rnorm(max_all_t),
      log_lambda_gp_02_level_alpha = if (enable_ms_02) rep(1.0, n_levels) else numeric(0),
      log_lambda_gp_02_level_rho = if (enable_ms_02) invgamma::rinvgamma(n_levels, log_lambda_gp_02_level_rho_alpha, log_lambda_gp_02_level_rho_beta) else numeric(0),
      log_lambda_gp_02_level_intercept_sd = if (enable_ms_02) abs(rnorm(n_levels, sd = log_lambda_gp_02_level_intercept_sd_sd)) else numeric(0),
      log_lambda_gp_02_level_eta = if (enable_ms_02 && n_enabled_groups_ms_baseline_02 > 0) {
        matrix(rnorm(n_enabled_groups_ms_baseline_02 * max_all_t), nrow = n_enabled_groups_ms_baseline_02, ncol = max_all_t)
      },
      raw_log_lambda_gp_02_level_intercept = if (n_enabled_groups_ms_baseline_02 > 0) {
        rnorm(n_enabled_groups_ms_baseline_02)
      },

      # --- 1→2 Transition: Sojourn Time GP ---
      log_lambda_gp_12_s_pop_intercept = if (need_12_s_gp) {
        array(rnorm(1, log_lambda_gp_12_s_pop_intercept_mean, log_lambda_gp_12_s_pop_intercept_sd), dim = 1)
      },
      log_lambda_gp_12_s_pop_alpha = if (need_12_s_gp) array(1.0, dim = 1),
      log_lambda_gp_12_s_pop_rho = if (need_12_s_gp) {
        array(invgamma::rinvgamma(1, log_lambda_gp_12_s_pop_rho_alpha, log_lambda_gp_12_s_pop_rho_beta), dim = 1)
      },
      log_lambda_gp_12_s_pop_eta = if (need_12_s_gp) rnorm(ms_max_sojourn_t),
      log_lambda_gp_12_s_level_alpha = if (need_12_s_gp) rep(1.0, n_levels) else numeric(0),
      log_lambda_gp_12_s_level_rho = if (need_12_s_gp) invgamma::rinvgamma(n_levels, log_lambda_gp_12_s_level_rho_alpha, log_lambda_gp_12_s_level_rho_beta) else numeric(0),
      log_lambda_gp_12_s_level_intercept_sd = if (need_12_s_gp) abs(rnorm(n_levels, sd = log_lambda_gp_12_s_level_intercept_sd_sd)) else numeric(0),
      log_lambda_gp_12_s_level_eta = if (need_12_s_gp && n_enabled_groups_ms_baseline_12_s > 0) {
        matrix(rnorm(n_enabled_groups_ms_baseline_12_s * ms_max_sojourn_t), nrow = n_enabled_groups_ms_baseline_12_s, ncol = ms_max_sojourn_t)
      },
      raw_log_lambda_gp_12_s_level_intercept = if (n_enabled_groups_ms_baseline_12_s > 0) {
        rnorm(n_enabled_groups_ms_baseline_12_s)
      },

      # --- 1→2 Transition: Clock-forward Time GP ---
      log_lambda_gp_12_t_pop_intercept = if (need_12_t_gp) {
        array(rnorm(1, log_lambda_gp_12_t_pop_intercept_mean, log_lambda_gp_12_t_pop_intercept_sd), dim = 1)
      },
      log_lambda_gp_12_t_pop_alpha = if (need_12_t_gp) array(1.0, dim = 1),
      log_lambda_gp_12_t_pop_rho = if (need_12_t_gp) {
        array(invgamma::rinvgamma(1, log_lambda_gp_12_t_pop_rho_alpha, log_lambda_gp_12_t_pop_rho_beta), dim = 1)
      },
      log_lambda_gp_12_t_pop_eta = if (need_12_t_gp) rnorm(max_all_t),
      log_lambda_gp_12_t_level_alpha = if (need_12_t_gp) rep(1.0, n_levels) else numeric(0),
      log_lambda_gp_12_t_level_rho = if (need_12_t_gp) invgamma::rinvgamma(n_levels, log_lambda_gp_12_t_level_rho_alpha, log_lambda_gp_12_t_level_rho_beta) else numeric(0),
      log_lambda_gp_12_t_level_intercept_sd = if (need_12_t_gp) abs(rnorm(n_levels, sd = log_lambda_gp_12_t_level_intercept_sd_sd)) else numeric(0),
      log_lambda_gp_12_t_level_eta = if (need_12_t_gp && n_enabled_groups_ms_baseline_12_t > 0) {
        matrix(rnorm(n_enabled_groups_ms_baseline_12_t * max_all_t), nrow = n_enabled_groups_ms_baseline_12_t, ncol = max_all_t)
      },
      raw_log_lambda_gp_12_t_level_intercept = if (n_enabled_groups_ms_baseline_12_t > 0) {
        rnorm(n_enabled_groups_ms_baseline_12_t)
      },

      # --- Time-varying Covariate Coefficients ---
      time_varying_coef_01 = if (enable_ms_01 && enable_ms_pop_time_varying_cov && n_time_varying_covar > 0) {
        rnorm(n_time_varying_covar, time_varying_coef_01_mean, time_varying_coef_01_sd)
      },
      time_varying_coef_02 = if (enable_ms_02 && enable_ms_pop_time_varying_cov && n_time_varying_covar > 0) {
        rnorm(n_time_varying_covar, time_varying_coef_02_mean, time_varying_coef_02_sd)
      },
      time_varying_coef_12 = if (enable_ms_12 && enable_ms_pop_time_varying_cov && n_time_varying_covar > 0) {
        rnorm(n_time_varying_covar, time_varying_coef_12_mean, time_varying_coef_12_sd)
      },

      # --- Time-invariant Covariate Coefficients (QR space) ---
      time_invariant_coef_qr_01 = if (enable_ms_01 && enable_ms_pop_time_invariant_cov && n_time_invariant_covar > 0) {
        rnorm(n_time_invariant_covar, time_invariant_coef_01_mean, time_invariant_coef_01_sd)
      },
      time_invariant_coef_qr_02 = if (enable_ms_02 && enable_ms_pop_time_invariant_cov && n_time_invariant_covar > 0) {
        rnorm(n_time_invariant_covar, time_invariant_coef_02_mean, time_invariant_coef_02_sd)
      },
      time_invariant_coef_qr_12 = if (enable_ms_12 && enable_ms_pop_time_invariant_cov && n_time_invariant_covar > 0) {
        rnorm(n_time_invariant_covar, time_invariant_coef_12_mean, time_invariant_coef_12_sd)
      },

      # --- Multi-level Random Slope SDs ---
      sd_level_slope_01 = if (enable_ms_01 && n_time_invariant_covar > 0) {
        lapply(seq_len(n_levels), function(lv) abs(rnorm(n_time_invariant_covar, sd = sd_level_slope_01_sd[[lv]])))
      },
      sd_level_slope_02 = if (enable_ms_02 && n_time_invariant_covar > 0) {
        lapply(seq_len(n_levels), function(lv) abs(rnorm(n_time_invariant_covar, sd = sd_level_slope_02_sd[[lv]])))
      },
      sd_level_slope_12 = if (enable_ms_12 && n_time_invariant_covar > 0) {
        lapply(seq_len(n_levels), function(lv) abs(rnorm(n_time_invariant_covar, sd = sd_level_slope_12_sd[[lv]])))
      },

      # --- Multi-level Raw Random Slopes ---
      raw_level_slope_01 = if (enable_ms_01 && n_time_invariant_covar > 0) {
        matrix(rnorm(n_enabled_groups_ms_slope * n_time_invariant_covar, sd = 0.5),
               nrow = n_enabled_groups_ms_slope, ncol = n_time_invariant_covar)
      },
      raw_level_slope_02 = if (enable_ms_02 && n_time_invariant_covar > 0) {
        matrix(rnorm(n_enabled_groups_ms_slope * n_time_invariant_covar, sd = 0.5),
               nrow = n_enabled_groups_ms_slope, ncol = n_time_invariant_covar)
      },
      raw_level_slope_12 = if (enable_ms_12 && n_time_invariant_covar > 0) {
        matrix(rnorm(n_enabled_groups_ms_slope * n_time_invariant_covar, sd = 0.5),
               nrow = n_enabled_groups_ms_slope, ncol = n_time_invariant_covar)
      },
    )
  })
}

# nolint end: object_usage_linter
