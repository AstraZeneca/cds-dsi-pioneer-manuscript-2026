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

    # --- Decomposed baseline-hazard schema (per-transition x per-level) ---
    # Reconstruct each slot's legacy 0-4 mode from the decomposed arrays so the
    # per-slot group counts mirror Stan's transformed_data.stan exactly.
    # Slot order matches flags.stan: 1=01,2=02,3=03,4=12_s,5=12_t,6=32.
    ms_slot_active <- c(
      as.integer(enable_ms_01), as.integer(enable_ms_02), as.integer(enable_ms_03),
      as.integer(need_12_s_gp), as.integer(need_12_t_gp), as.integer(enable_ms_32)
    )
    reconstruct_legacy_mode <- function(slot_idx) {
      if (!ms_slot_active[slot_idx]) return(rep(0L, n_levels))
      gp_row <- enable_ms_level_gp[slot_idx, ]
      mode_row <- ms_level_intercept_mode[slot_idx, ]
      dplyr::case_when(
        gp_row == 1L   ~ 3L,   # RE_GP
        mode_row == 0L ~ 0L,   # NONE
        mode_row == 1L ~ 1L,   # FE
        mode_row == 2L ~ 2L,   # RE_NCP
        TRUE           ~ 4L    # RE_CP (decomposed mode 3)
      )
    }
    ms_legacy_mode <- lapply(seq_len(6L), reconstruct_legacy_mode)

    n_enabled_slot <- vapply(ms_legacy_mode,
      function(m) sum(n_groups_per_level[m > 0L]), integer(1))
    n_gp_slot <- vapply(ms_legacy_mode,
      function(m) sum(n_groups_per_level[m == 3L]), integer(1))
    any_re_slot <- vapply(ms_legacy_mode,
      function(m) any(m == 2L | m == 3L | m == 4L), logical(1))

    # Enabled group counts: truthy (> 0) — both intercept-only and GP count.
    n_enabled_groups_ms_baseline_01 <- n_enabled_slot[1]
    n_enabled_groups_ms_baseline_02 <- n_enabled_slot[2]
    n_enabled_groups_ms_baseline_03 <- n_enabled_slot[3]
    n_enabled_groups_ms_baseline_12_s <- n_enabled_slot[4]
    n_enabled_groups_ms_baseline_12_t <- n_enabled_slot[5]
    n_enabled_groups_ms_baseline_32 <- n_enabled_slot[6]
    n_enabled_groups_ms_slope <- sum(n_groups_per_level[enable_ms_level_cov == 1])
    # any_re per slot (used to gate the per-slot intercept-SD inits)
    any_re_level_01   <- any_re_slot[1]
    any_re_level_02   <- any_re_slot[2]
    any_re_level_03   <- any_re_slot[3]
    any_re_level_12_s <- any_re_slot[4]
    any_re_level_12_t <- any_re_slot[5]
    any_re_level_32   <- any_re_slot[6]

    # GP-only group counts (legacy mode == 3): used for _level_eta matrix sizing
    n_gp_groups_ms_baseline_01 <- n_gp_slot[1]
    n_gp_groups_ms_baseline_02 <- n_gp_slot[2]
    n_gp_groups_ms_baseline_03 <- n_gp_slot[3]
    n_gp_groups_ms_baseline_12_s <- n_gp_slot[4]
    n_gp_groups_ms_baseline_12_t <- n_gp_slot[5]
    n_gp_groups_ms_baseline_32 <- n_gp_slot[6]

    # GP knot counts (coarse grid, matches Stan transformed_data ceiling division)
    n_ms_gp_cal_knots        <- ceiling(max_all_t / ms_gp_grid_step)
    n_ms_gp_sojourn_knots    <- ceiling(ms_max_sojourn_t / ms_gp_grid_step)
    n_ms_gp_sojourn_32_knots <- ceiling(ms_max_sojourn_t_32 / ms_gp_grid_step)

    # --- Correlated intercept block inits (Phase 2) ---
    # Patient (last) level uses the forecast-patient count, mirroring Stan's
    # n_forecast_groups_per_level. Empty when no block is configured.
    n_forecast_groups_per_level <- n_groups_per_level
    if (exists("n_forecast_patients", inherits = FALSE)) {
      n_forecast_groups_per_level[n_levels] <- n_forecast_patients
    }
    .ms_corr <- ms_corr_blocks(
      ms_level_intercept_corr_group, ms_slot_active, n_forecast_groups_per_level
    )
    .ms_corr_inits <- ms_corr_block_inits(.ms_corr, z_sd = 0.3)

    c(tibble::lst(
      # --- 0→1 Transition (Progression / PFS event) ---
      log_lambda_gp_01_pop_intercept = if (enable_ms_01) {
        array(rnorm(1, log_lambda_gp_01_pop_intercept_mean, log_lambda_gp_01_pop_intercept_sd), dim = 1)
      },
      log_lambda_gp_01_pop_alpha = if (enable_ms_01) array(1.0, dim = 1),
      log_lambda_gp_01_pop_rho = if (enable_ms_01) {
        array(max(invgamma::rinvgamma(1, log_lambda_gp_01_pop_rho_alpha, log_lambda_gp_01_pop_rho_beta), ms_gp_grid_step), dim = 1)
      },
      log_lambda_gp_01_pop_eta = if (enable_ms_01) rnorm(n_ms_gp_cal_knots),
      # 0->1 baseline log-time trend (Fix A, hierarchical — mirrors the GP):
      # population slope, per-group raw slope deviation (NCP, zeros — same length
      # as the GP intercept raw vector, n_enabled_groups_ms_baseline_01), and
      # per-level slope SD. All size-0 when the flag is off, matching the Stan
      # parameter sizing exactly.
      log_lambda_trend_01_pop_slope =
        if (isTRUE(enable_ms_baseline_trend_01 == 1L)) array(0.1, dim = 1) else numeric(0),
      raw_log_lambda_trend_01_level =
        if (isTRUE(enable_ms_baseline_trend_01 == 1L) && n_enabled_groups_ms_baseline_01 > 0) {
          as.array(rep(0, n_enabled_groups_ms_baseline_01))
        } else {
          numeric(0)
        },
      log_lambda_trend_01_level_sd =
        if (isTRUE(enable_ms_baseline_trend_01 == 1L) && any_re_level_01) {
          rep(0.1, n_levels)
        } else {
          numeric(0)
        },
      log_lambda_gp_01_level_alpha = if (enable_ms_01) rep(1.0, n_levels) else numeric(0),
      log_lambda_gp_01_level_rho = if (enable_ms_01) pmax(invgamma::rinvgamma(n_levels, log_lambda_gp_01_level_rho_alpha, log_lambda_gp_01_level_rho_beta), ms_gp_grid_step) else numeric(0),
      log_lambda_gp_01_level_intercept_sd = if (enable_ms_01 && any_re_level_01) abs(rnorm(n_levels, sd = log_lambda_gp_01_level_intercept_sd_sd)) else numeric(0),
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
      log_lambda_gp_02_level_intercept_sd = if (enable_ms_02 && any_re_level_02) abs(rnorm(n_levels, sd = log_lambda_gp_02_level_intercept_sd_sd)) else numeric(0),
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
      log_lambda_gp_12_s_level_intercept_sd = if (need_12_s_gp && any_re_level_12_s) abs(rnorm(n_levels, sd = log_lambda_gp_12_s_level_intercept_sd_sd)) else numeric(0),
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
      log_lambda_gp_12_t_level_intercept_sd = if (need_12_t_gp && any_re_level_12_t) abs(rnorm(n_levels, sd = log_lambda_gp_12_t_level_intercept_sd_sd)) else numeric(0),
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
      log_lambda_gp_03_level_intercept_sd = if (enable_ms_03 && any_re_level_03) {
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
      log_lambda_gp_32_s_level_intercept_sd = if (enable_ms_32 && any_re_level_32) {
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
      # Compute sizes locally from flags. 0->1 dimensioning matches
      # stan/modules/multistate/parameters.stan: single coefficient only in
      # observed visit-gated mode; full n_time_varying_covar in continuous
      # and latent visit-gated modes.
      time_varying_coef_01 = {
        observed_vg <- isTRUE(enable_ms_visit_gated_01 == 1L) &&
          !isTRUE(enable_ms_visit_gated_latent_01 == 1L)
        n_tv_01 <- if (observed_vg) 1L else n_time_varying_covar
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
      time_varying_coef_03 = {
        n_tv_03 <- if (isTRUE(enable_ms_03_time_varying_cov == 1L)) n_time_varying_covar else 0L
        if (enable_ms_03 && enable_ms_pop_time_varying_cov && n_tv_03 > 0)
          as.array(rnorm(n_tv_03, time_varying_coef_03_mean, time_varying_coef_03_sd * qr_init_scale))
        else NULL
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
      time_invariant_coef_qr_03 = if (enable_ms_03 && enable_ms_pop_time_invariant_cov &&
                                       isTRUE(enable_ms_03_time_invariant_cov == 1L) &&
                                       n_time_invariant_covar > 0) {
        rnorm(n_time_invariant_covar, time_invariant_coef_03_mean, time_invariant_coef_03_sd * qr_init_scale)
      },
      time_invariant_coef_qr_32 = if (enable_ms_32 && enable_ms_pop_time_invariant_cov &&
                                       isTRUE(enable_ms_32_time_invariant_cov == 1L) &&
                                       n_time_invariant_covar > 0) {
        rnorm(n_time_invariant_covar, time_invariant_coef_32_mean, time_invariant_coef_32_sd * qr_init_scale)
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
      sd_level_slope_03 = if (enable_ms_03 && isTRUE(enable_ms_03_time_invariant_cov == 1L) && n_time_invariant_covar > 0) {
        lapply(seq_len(n_levels), function(lv) abs(rnorm(n_time_invariant_covar, sd = sd_level_slope_03_sd[[lv]])))
      },
      sd_level_slope_32 = if (enable_ms_32 && isTRUE(enable_ms_32_time_invariant_cov == 1L) && n_time_invariant_covar > 0) {
        lapply(seq_len(n_levels), function(lv) abs(rnorm(n_time_invariant_covar, sd = sd_level_slope_32_sd[[lv]])))
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
      raw_level_slope_03 = if (enable_ms_03 && isTRUE(enable_ms_03_time_invariant_cov == 1L) && n_time_invariant_covar > 0) {
        matrix(rnorm(n_enabled_groups_ms_slope * n_time_invariant_covar, sd = 0.5),
               nrow = n_enabled_groups_ms_slope, ncol = n_time_invariant_covar)
      },
      raw_level_slope_32 = if (enable_ms_32 && isTRUE(enable_ms_32_time_invariant_cov == 1L) && n_time_invariant_covar > 0) {
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
      coef_log_entry_covar_12 = if (enable_ms_12 && isTRUE(enable_ms_12_entry_covar == 1L)) {
        array(rnorm(1, coef_log_entry_covar_12_mean, coef_log_entry_covar_12_sd * 0.3), dim = 1)
      },
      coef_log_entry_covar_32 = if (enable_ms_32 && isTRUE(enable_ms_32_entry_covar == 1L)) {
        array(rnorm(1, coef_log_entry_covar_32_mean, coef_log_entry_covar_32_sd * 0.3), dim = 1)
      },
    ), .ms_corr_inits)
  })
}

# nolint end: object_usage_linter
