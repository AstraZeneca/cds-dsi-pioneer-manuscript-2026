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

    # Enabled group counts: truthy (> 0) — both intercept-only (1) and GP (2) count
    n_enabled_groups_ms_baseline <- sum(n_groups_per_level[enable_ms_level_baseline_hazard > 0])
    n_enabled_groups_ms_baseline_01 <- if (enable_ms_01) n_enabled_groups_ms_baseline else 0L
    n_enabled_groups_ms_baseline_02 <- if (enable_ms_02) n_enabled_groups_ms_baseline else 0L
    n_enabled_groups_ms_baseline_12_s <- if (need_12_s_gp) n_enabled_groups_ms_baseline else 0L
    n_enabled_groups_ms_baseline_12_t <- if (need_12_t_gp) n_enabled_groups_ms_baseline else 0L
    n_enabled_groups_ms_baseline_03 <- if (enable_ms_03) n_enabled_groups_ms_baseline else 0L
    n_enabled_groups_ms_baseline_32 <- if (enable_ms_32) n_enabled_groups_ms_baseline else 0L
    n_enabled_groups_ms_slope <- sum(n_groups_per_level[enable_ms_level_cov == 1])
    any_re_level <- any(enable_ms_level_baseline_hazard >= 2L)

    # GP-only group counts (mode == 3): used for _level_eta matrix sizing
    n_gp_groups_ms_baseline <- sum(n_groups_per_level[enable_ms_level_baseline_hazard == 3L])
    n_gp_groups_ms_baseline_01 <- if (enable_ms_01) n_gp_groups_ms_baseline else 0L
    n_gp_groups_ms_baseline_02 <- if (enable_ms_02) n_gp_groups_ms_baseline else 0L
    n_gp_groups_ms_baseline_12_s <- if (need_12_s_gp) n_gp_groups_ms_baseline else 0L
    n_gp_groups_ms_baseline_12_t <- if (need_12_t_gp) n_gp_groups_ms_baseline else 0L
    n_gp_groups_ms_baseline_03 <- if (enable_ms_03) n_gp_groups_ms_baseline else 0L
    n_gp_groups_ms_baseline_32 <- if (enable_ms_32) n_gp_groups_ms_baseline else 0L

    # GP knot counts (coarse grid, matches Stan transformed_data ceiling division)
    n_ms_gp_cal_knots        <- ceiling(max_all_t / ms_gp_grid_step)
    n_ms_gp_sojourn_knots    <- ceiling(ms_max_sojourn_t / ms_gp_grid_step)
    n_ms_gp_sojourn_32_knots <- ceiling(ms_max_sojourn_t_32 / ms_gp_grid_step)

    tibble::lst(
      # --- 0→1 Transition (Progression / PFS event) ---
      log_lambda_gp_01_pop_intercept = if (enable_ms_01) {
        array(rnorm(1, log_lambda_gp_01_pop_intercept_mean, log_lambda_gp_01_pop_intercept_sd), dim = 1)
      },
      log_lambda_gp_01_pop_alpha = if (enable_ms_01) array(1.0, dim = 1),
      log_lambda_gp_01_pop_rho = if (enable_ms_01) {
        array(max(invgamma::rinvgamma(1, log_lambda_gp_01_pop_rho_alpha, log_lambda_gp_01_pop_rho_beta), ms_gp_grid_step), dim = 1)
      },
      log_lambda_gp_01_pop_eta = if (enable_ms_01) rnorm(n_ms_gp_cal_knots),
      log_lambda_gp_01_level_alpha = if (enable_ms_01) rep(1.0, n_levels) else numeric(0),
      log_lambda_gp_01_level_rho = if (enable_ms_01) pmax(invgamma::rinvgamma(n_levels, log_lambda_gp_01_level_rho_alpha, log_lambda_gp_01_level_rho_beta), ms_gp_grid_step) else numeric(0),
      log_lambda_gp_01_level_intercept_sd = if (enable_ms_01 && any_re_level) abs(rnorm(n_levels, sd = log_lambda_gp_01_level_intercept_sd_sd)) else numeric(0),
      log_lambda_gp_01_level_eta = if (enable_ms_01 && n_gp_groups_ms_baseline_01 > 0) {
        matrix(rnorm(n_gp_groups_ms_baseline_01 * n_ms_gp_cal_knots), nrow = n_gp_groups_ms_baseline_01, ncol = n_ms_gp_cal_knots)
      },
      raw_log_lambda_gp_01_level_intercept = if (n_enabled_groups_ms_baseline_01 > 0) {
        as.array(rnorm(n_enabled_groups_ms_baseline_01))
      },

      # --- 0→2 Transition (Death without progression) ---
      log_lambda_gp_02_pop_intercept = if (enable_ms_02) {
        array(rnorm(1, log_lambda_gp_02_pop_intercept_mean, log_lambda_gp_02_pop_intercept_sd), dim = 1)
      },
      log_lambda_gp_02_pop_alpha = if (enable_ms_02 && !isTRUE(share_dead_gp_shape == 1L)) array(1.0, dim = 1),
      log_lambda_gp_02_pop_rho = if (enable_ms_02 && !isTRUE(share_dead_gp_shape == 1L)) {
        array(max(invgamma::rinvgamma(1, log_lambda_gp_02_pop_rho_alpha, log_lambda_gp_02_pop_rho_beta), ms_gp_grid_step), dim = 1)
      },
      log_lambda_gp_02_pop_eta = if (enable_ms_02 && !isTRUE(share_dead_gp_shape == 1L)) rnorm(n_ms_gp_cal_knots),

      # --- Shared "dead" GP shape (0→2 + 1→2 clock-forward, when share_dead_gp_shape=1) ---
      log_lambda_gp_dead_pop_alpha = if (isTRUE(share_dead_gp_shape == 1L)) array(1.0, dim = 1),
      log_lambda_gp_dead_pop_rho = if (isTRUE(share_dead_gp_shape == 1L)) {
        array(max(invgamma::rinvgamma(1, log_lambda_gp_dead_pop_rho_alpha,
                                        log_lambda_gp_dead_pop_rho_beta), ms_gp_grid_step), dim = 1)
      },
      log_lambda_gp_dead_pop_eta = if (isTRUE(share_dead_gp_shape == 1L)) rnorm(n_ms_gp_cal_knots),
      log_lambda_gp_02_level_alpha = if (enable_ms_02) rep(1.0, n_levels) else numeric(0),
      log_lambda_gp_02_level_rho = if (enable_ms_02) pmax(invgamma::rinvgamma(n_levels, log_lambda_gp_02_level_rho_alpha, log_lambda_gp_02_level_rho_beta), ms_gp_grid_step) else numeric(0),
      log_lambda_gp_02_level_intercept_sd = if (enable_ms_02 && any_re_level) abs(rnorm(n_levels, sd = log_lambda_gp_02_level_intercept_sd_sd)) else numeric(0),
      log_lambda_gp_02_level_eta = if (enable_ms_02 && n_gp_groups_ms_baseline_02 > 0) {
        matrix(rnorm(n_gp_groups_ms_baseline_02 * n_ms_gp_cal_knots), nrow = n_gp_groups_ms_baseline_02, ncol = n_ms_gp_cal_knots)
      },
      raw_log_lambda_gp_02_level_intercept = if (n_enabled_groups_ms_baseline_02 > 0) {
        as.array(rnorm(n_enabled_groups_ms_baseline_02))
      },

      # --- 1→2 Transition: Sojourn Time GP ---
      log_lambda_gp_12_s_pop_intercept = if (need_12_s_gp) {
        array(rnorm(1, log_lambda_gp_12_s_pop_intercept_mean, log_lambda_gp_12_s_pop_intercept_sd), dim = 1)
      },
      log_lambda_gp_12_s_pop_alpha = if (need_12_s_gp) array(1.0, dim = 1),
      log_lambda_gp_12_s_pop_rho = if (need_12_s_gp) {
        array(max(invgamma::rinvgamma(1, log_lambda_gp_12_s_pop_rho_alpha, log_lambda_gp_12_s_pop_rho_beta), ms_gp_grid_step), dim = 1)
      },
      log_lambda_gp_12_s_pop_eta = if (need_12_s_gp) rnorm(n_ms_gp_sojourn_knots),
      log_lambda_gp_12_s_level_alpha = if (need_12_s_gp) rep(1.0, n_levels) else numeric(0),
      log_lambda_gp_12_s_level_rho = if (need_12_s_gp) pmax(invgamma::rinvgamma(n_levels, log_lambda_gp_12_s_level_rho_alpha, log_lambda_gp_12_s_level_rho_beta), ms_gp_grid_step) else numeric(0),
      log_lambda_gp_12_s_level_intercept_sd = if (need_12_s_gp && any_re_level) abs(rnorm(n_levels, sd = log_lambda_gp_12_s_level_intercept_sd_sd)) else numeric(0),
      log_lambda_gp_12_s_level_eta = if (need_12_s_gp && n_gp_groups_ms_baseline_12_s > 0) {
        matrix(rnorm(n_gp_groups_ms_baseline_12_s * n_ms_gp_sojourn_knots), nrow = n_gp_groups_ms_baseline_12_s, ncol = n_ms_gp_sojourn_knots)
      },
      raw_log_lambda_gp_12_s_level_intercept = if (n_enabled_groups_ms_baseline_12_s > 0) {
        as.array(rnorm(n_enabled_groups_ms_baseline_12_s))
      },

      # --- 1→2 Transition: Clock-forward Time GP ---
      log_lambda_gp_12_t_pop_intercept = if (need_12_t_gp) {
        array(rnorm(1, log_lambda_gp_12_t_pop_intercept_mean, log_lambda_gp_12_t_pop_intercept_sd), dim = 1)
      },
      log_lambda_gp_12_t_pop_alpha = if (need_12_t_gp && !isTRUE(share_dead_gp_shape == 1L)) array(1.0, dim = 1),
      log_lambda_gp_12_t_pop_rho = if (need_12_t_gp && !isTRUE(share_dead_gp_shape == 1L)) {
        array(max(invgamma::rinvgamma(1, log_lambda_gp_12_t_pop_rho_alpha, log_lambda_gp_12_t_pop_rho_beta), ms_gp_grid_step), dim = 1)
      },
      log_lambda_gp_12_t_pop_eta = if (need_12_t_gp && !isTRUE(share_dead_gp_shape == 1L)) rnorm(n_ms_gp_cal_knots),
      log_lambda_gp_12_t_level_alpha = if (need_12_t_gp) rep(1.0, n_levels) else numeric(0),
      log_lambda_gp_12_t_level_rho = if (need_12_t_gp) pmax(invgamma::rinvgamma(n_levels, log_lambda_gp_12_t_level_rho_alpha, log_lambda_gp_12_t_level_rho_beta), ms_gp_grid_step) else numeric(0),
      log_lambda_gp_12_t_level_intercept_sd = if (need_12_t_gp && any_re_level) abs(rnorm(n_levels, sd = log_lambda_gp_12_t_level_intercept_sd_sd)) else numeric(0),
      log_lambda_gp_12_t_level_eta = if (need_12_t_gp && n_gp_groups_ms_baseline_12_t > 0) {
        matrix(rnorm(n_gp_groups_ms_baseline_12_t * n_ms_gp_cal_knots), nrow = n_gp_groups_ms_baseline_12_t, ncol = n_ms_gp_cal_knots)
      },
      raw_log_lambda_gp_12_t_level_intercept = if (n_enabled_groups_ms_baseline_12_t > 0) {
        as.array(rnorm(n_enabled_groups_ms_baseline_12_t))
      },

      # --- 0→3 Transition (Dropout: GP baseline hazard with N-level hierarchy) ---
      log_lambda_gp_03_pop_intercept = if (enable_ms_03) {
        array(rnorm(1, log_lambda_gp_03_pop_intercept_mean, log_lambda_gp_03_pop_intercept_sd), dim = 1)
      },
      log_lambda_gp_03_pop_alpha = if (enable_ms_03) array(1.0, dim = 1),
      log_lambda_gp_03_pop_rho = if (enable_ms_03) {
        array(max(invgamma::rinvgamma(1, log_lambda_gp_03_pop_rho_alpha, log_lambda_gp_03_pop_rho_beta), ms_gp_grid_step), dim = 1)
      },
      log_lambda_gp_03_pop_eta = if (enable_ms_03) rnorm(n_ms_gp_cal_knots),
      log_lambda_gp_03_level_alpha = if (enable_ms_03) rep(1.0, n_levels) else numeric(0),
      log_lambda_gp_03_level_rho = if (enable_ms_03) {
        pmax(invgamma::rinvgamma(n_levels, log_lambda_gp_03_level_rho_alpha, log_lambda_gp_03_level_rho_beta), ms_gp_grid_step)
      } else numeric(0),
      log_lambda_gp_03_level_intercept_sd = if (enable_ms_03 && any_re_level) {
        abs(rnorm(n_levels, sd = log_lambda_gp_03_level_intercept_sd_sd))
      } else numeric(0),
      log_lambda_gp_03_level_eta = if (enable_ms_03 && n_gp_groups_ms_baseline_03 > 0) {
        matrix(rnorm(n_gp_groups_ms_baseline_03 * n_ms_gp_cal_knots),
               nrow = n_gp_groups_ms_baseline_03, ncol = n_ms_gp_cal_knots)
      },
      raw_log_lambda_gp_03_level_intercept = if (n_enabled_groups_ms_baseline_03 > 0) {
        as.array(rnorm(n_enabled_groups_ms_baseline_03))
      },

      # --- 3→2 Transition (Off-trial death: sojourn time GP, semi-Markov) ---
      log_lambda_gp_32_s_pop_intercept = if (enable_ms_32) {
        array(rnorm(1, log_lambda_gp_32_s_pop_intercept_mean, log_lambda_gp_32_s_pop_intercept_sd), dim = 1)
      },
      log_lambda_gp_32_s_pop_alpha = if (enable_ms_32) array(1.0, dim = 1),
      log_lambda_gp_32_s_pop_rho = if (enable_ms_32) {
        array(max(invgamma::rinvgamma(1, log_lambda_gp_32_s_pop_rho_alpha, log_lambda_gp_32_s_pop_rho_beta), ms_gp_grid_step), dim = 1)
      },
      log_lambda_gp_32_s_pop_eta = if (enable_ms_32) rnorm(n_ms_gp_sojourn_32_knots),
      log_lambda_gp_32_s_level_alpha = if (enable_ms_32) rep(1.0, n_levels) else numeric(0),
      log_lambda_gp_32_s_level_rho = if (enable_ms_32) {
        pmax(invgamma::rinvgamma(n_levels, log_lambda_gp_32_s_level_rho_alpha, log_lambda_gp_32_s_level_rho_beta), ms_gp_grid_step)
      } else numeric(0),
      log_lambda_gp_32_s_level_intercept_sd = if (enable_ms_32 && any_re_level) {
        abs(rnorm(n_levels, sd = log_lambda_gp_32_s_level_intercept_sd_sd))
      } else numeric(0),
      log_lambda_gp_32_s_level_eta = if (enable_ms_32 && n_gp_groups_ms_baseline_32 > 0) {
        matrix(rnorm(n_gp_groups_ms_baseline_32 * n_ms_gp_sojourn_32_knots),
               nrow = n_gp_groups_ms_baseline_32, ncol = n_ms_gp_sojourn_32_knots)
      },
      raw_log_lambda_gp_32_s_level_intercept = if (n_enabled_groups_ms_baseline_32 > 0) {
        as.array(rnorm(n_enabled_groups_ms_baseline_32))
      },

      # --- Time-varying Covariate Coefficients ---
      # Use a conservative init scale: Q row norms ~ sqrt(n_patients), so even
      # prior-SD draws can produce Q*coef >> hazard cap. Scale to keep Q*coef < 10.
      qr_init_scale = min(1.0, 10.0 / (sqrt(n_patients - 1) * sqrt(max(n_time_invariant_covar, 1L)))),
      # Compute sizes locally from flags
      time_varying_coef_01 = {
        n_tv_01 <- if (isTRUE(enable_ms_visit_gated_01 == 1L)) 1L else n_time_varying_covar
        if (enable_ms_01 && enable_ms_pop_time_varying_cov && n_tv_01 > 0)
          as.array(rnorm(n_tv_01, time_varying_coef_01_mean, time_varying_coef_01_sd * qr_init_scale))
        else NULL
      },
      time_varying_coef_02 = {
        n_tv_02 <- if (isTRUE(enable_ms_02_time_varying_cov == 1L)) n_time_varying_covar else 0L
        if (enable_ms_02 && enable_ms_pop_time_varying_cov && n_tv_02 > 0)
          as.array(rnorm(n_tv_02, time_varying_coef_02_mean, time_varying_coef_02_sd * qr_init_scale))
        else NULL
      },
      time_varying_coef_12 = if (enable_ms_12 && enable_ms_pop_time_varying_cov && n_time_varying_covar > 0) {
        as.array(rnorm(n_time_varying_covar, time_varying_coef_12_mean, time_varying_coef_12_sd * qr_init_scale))
      },

      # --- Time-invariant Covariate Coefficients (QR space) ---
      time_invariant_coef_qr_01 = if (enable_ms_01 && enable_ms_pop_time_invariant_cov && n_time_invariant_covar > 0) {
        rnorm(n_time_invariant_covar, time_invariant_coef_01_mean, time_invariant_coef_01_sd * qr_init_scale)
      },
      time_invariant_coef_qr_02 = if (enable_ms_02 && enable_ms_pop_time_invariant_cov && n_time_invariant_covar > 0) {
        rnorm(n_time_invariant_covar, time_invariant_coef_02_mean, time_invariant_coef_02_sd * qr_init_scale)
      },
      time_invariant_coef_qr_12 = if (enable_ms_12 && enable_ms_pop_time_invariant_cov && n_time_invariant_covar > 0) {
        rnorm(n_time_invariant_covar, time_invariant_coef_12_mean, time_invariant_coef_12_sd * qr_init_scale)
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

      # Student-t hierarchy: multistate nu starting values
      ms_nu_baseline_level = if (isTRUE(enable_student_t_hierarchy == 1L)) {
        as.array(pmax(2.1, rgamma(n_levels, 2, 0.1)))
      },
      ms_nu_slope_level = if (isTRUE(enable_student_t_hierarchy == 1L)) {
        as.array(pmax(2.1, rgamma(n_levels, 2, 0.1)))
      },

      # PSA-at-state-entry covariate coefficients
      coef_log_psa_12 = if (enable_ms_12 && isTRUE(enable_ms_12_entry_psa_cov == 1L)) {
        array(rnorm(1, coef_log_psa_12_mean, coef_log_psa_12_sd * 0.3), dim = 1)
      },
      coef_log_psa_32 = if (enable_ms_32 && isTRUE(enable_ms_32_entry_psa_cov == 1L)) {
        array(rnorm(1, coef_log_psa_32_mean, coef_log_psa_32_sd * 0.3), dim = 1)
      },
    )
  })
}

# nolint end: object_usage_linter
