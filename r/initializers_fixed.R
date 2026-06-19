# Fixed initializers for testing chain stability
# Uses reproducible per-chain values with realistic heterogeneity

# =========================================================================
# Shared multistate initialization helper
# =========================================================================
# Called by both the joint model (create_tumor_ssls_initializer_fixed) and
# the standalone multistate model (create_ms_standalone_initializer_fixed).
# Must be called inside a with(stan_data, { ... }) block so that all
# multistate fields (flags, hyperparameters, dimensions) are in scope.

ms_init_values_fixed <- function(env) {
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

    # Mirrors Stan's n_forecast_groups_per_level: patient level uses forecast
    # count, not the full patient count (background patients have no parameters).
    n_fgpl <- n_groups_per_level
    if (exists("n_forecast_patients", inherits = FALSE)) {
      n_fgpl[n_levels] <- n_forecast_patients
    }

    n_enabled_slot <- vapply(ms_legacy_mode,
      function(m) sum(n_fgpl[m > 0L]), integer(1))
    n_gp_slot <- vapply(ms_legacy_mode,
      function(m) sum(n_fgpl[m == 3L]), integer(1))
    # raw/cp split mirrors Stan's split_cp_ncp_pos (hierarchy.stanfunctions):
    # raw bucket = modes {FE=1, RE=2, RE_GP=3}; cp bucket = mode {RE_CP=4}.
    n_raw_slot <- vapply(ms_legacy_mode,
      function(m) sum(n_fgpl[m == 1L | m == 2L | m == 3L]), integer(1))
    n_cp_slot <- vapply(ms_legacy_mode,
      function(m) sum(n_fgpl[m == 4L]), integer(1))
    any_re_slot <- vapply(ms_legacy_mode,
      function(m) any(m == 2L | m == 3L | m == 4L), logical(1))

    # Enabled group counts - truthy (> 0) for both intercept-only and GP modes
    n_enabled_groups_ms_baseline_01 <- n_enabled_slot[1]
    n_enabled_groups_ms_baseline_02 <- n_enabled_slot[2]
    n_enabled_groups_ms_baseline_03 <- n_enabled_slot[3]
    n_enabled_groups_ms_baseline_12_s <- n_enabled_slot[4]
    n_enabled_groups_ms_baseline_12_t <- n_enabled_slot[5]
    n_enabled_groups_ms_baseline_32 <- n_enabled_slot[6]
    n_enabled_groups_ms_slope <- sum(n_fgpl[enable_ms_level_cov == 1])

    # GP-only group counts (legacy mode == 3, for eta matrix sizing)
    n_gp_groups_ms_baseline_01 <- n_gp_slot[1]
    n_gp_groups_ms_baseline_02 <- n_gp_slot[2]
    n_gp_groups_ms_baseline_03 <- n_gp_slot[3]
    n_gp_groups_ms_baseline_12_s <- n_gp_slot[4]
    n_gp_groups_ms_baseline_12_t <- n_gp_slot[5]
    n_gp_groups_ms_baseline_32 <- n_gp_slot[6]

    # raw/cp group counts: size the raw_*_level_intercept (NCP) and
    # cp_*_level_intercept (centered) vectors respectively. Mirror the per-slot
    # scalar aliases in stan/modules/multistate/transformed_data.stan.
    n_raw_groups_ms_baseline_01 <- n_raw_slot[1]
    n_raw_groups_ms_baseline_02 <- n_raw_slot[2]
    n_raw_groups_ms_baseline_03 <- n_raw_slot[3]
    n_raw_groups_ms_baseline_12_s <- n_raw_slot[4]
    n_raw_groups_ms_baseline_12_t <- n_raw_slot[5]
    n_raw_groups_ms_baseline_32 <- n_raw_slot[6]
    n_cp_groups_ms_baseline_01 <- n_cp_slot[1]
    n_cp_groups_ms_baseline_02 <- n_cp_slot[2]
    n_cp_groups_ms_baseline_03 <- n_cp_slot[3]
    n_cp_groups_ms_baseline_12_s <- n_cp_slot[4]
    n_cp_groups_ms_baseline_12_t <- n_cp_slot[5]
    n_cp_groups_ms_baseline_32 <- n_cp_slot[6]

    # any_re per slot (used to gate the per-slot intercept-SD inits)
    any_re_level_01   <- any_re_slot[1]
    any_re_level_02   <- any_re_slot[2]
    any_re_level_03   <- any_re_slot[3]
    any_re_level_12_s <- any_re_slot[4]
    any_re_level_12_t <- any_re_slot[5]
    any_re_level_32   <- any_re_slot[6]

    # GP knot counts (matches Stan transformed_data ceiling division)
    n_ms_gp_cal_knots        <- ceiling(max_all_t / ms_gp_grid_step)
    n_ms_gp_sojourn_knots    <- ceiling(ms_max_sojourn_t / ms_gp_grid_step)
    n_ms_gp_sojourn_32_knots <- ceiling(ms_max_sojourn_t_32 / ms_gp_grid_step)

    # --- Correlated intercept block inits (Phase 2) ---
    .ms_corr <- ms_corr_blocks(
      ms_level_intercept_corr_group, ms_slot_active, n_fgpl
    )
    .ms_corr_inits <- ms_corr_block_inits(.ms_corr, z_sd = 0.3)

    c(tibble::lst(
      # --- 0→1 Transition ---
      log_lambda_gp_01_pop_intercept = if (enable_ms_01) array(-4.5, dim = 1),
      log_lambda_gp_01_pop_alpha = if (enable_ms_01) array(1.0, dim = 1),
      log_lambda_gp_01_pop_rho = if (enable_ms_01) array(1.4, dim = 1),
      log_lambda_gp_01_pop_eta = if (enable_ms_01) rnorm(n_ms_gp_cal_knots, sd = 0.1),
      # 0->1 baseline log-time trend (Fix A, hierarchical — mirrors the GP):
      # population slope, per-group raw slope deviation (NCP, zeros — same length
      # as the GP intercept raw vector, n_enabled_groups_ms_baseline_01), and
      # per-level slope SD. All size-0 when the flag is off, matching the Stan
      # parameter sizing exactly.
      log_lambda_trend_01_pop_slope =
        if (isTRUE(enable_ms_baseline_trend_01 == 1L)) array(0.1, dim = 1) else numeric(0),
      raw_log_lambda_trend_01_level =
        if (isTRUE(enable_ms_baseline_trend_01 == 1L) && n_raw_groups_ms_baseline_01 > 0) {
          rep(0, n_raw_groups_ms_baseline_01)
        } else {
          numeric(0)
        },
      log_lambda_trend_01_level_sd =
        if (isTRUE(enable_ms_baseline_trend_01 == 1L) && any_re_level_01) {
          rep(0.1, n_levels)
        } else {
          numeric(0)
        },
      log_lambda_gp_01_level_alpha = rep(1.0, n_levels),
      log_lambda_gp_01_level_rho = rep(1.4, n_levels),
      log_lambda_gp_01_level_intercept_sd = if (any_re_level_01) rep(0.1, n_levels) else numeric(0),
      log_lambda_gp_01_level_eta = if (enable_ms_01 && n_gp_groups_ms_baseline_01 > 0) {
        matrix(rnorm(n_gp_groups_ms_baseline_01 * n_ms_gp_cal_knots, sd = 0.1),
               nrow = n_gp_groups_ms_baseline_01, ncol = n_ms_gp_cal_knots)
      },
      raw_log_lambda_gp_01_level_intercept = if (n_raw_groups_ms_baseline_01 > 0) rep(0, n_raw_groups_ms_baseline_01),
      # cp path (RE_CP, mode 4): centered, sampled ~ normal(0, sd[lv]); fixed start at 0.
      cp_log_lambda_gp_01_level_intercept = if (n_cp_groups_ms_baseline_01 > 0) rep(0, n_cp_groups_ms_baseline_01),

      # --- 0→2 Transition ---
      log_lambda_gp_02_pop_intercept = if (enable_ms_02) array(-4.5, dim = 1),
      log_lambda_gp_02_pop_alpha = if (enable_ms_02) array(1.0, dim = 1),
      log_lambda_gp_02_pop_rho = if (enable_ms_02) array(1.4, dim = 1),
      log_lambda_gp_02_pop_eta = if (enable_ms_02) rnorm(n_ms_gp_cal_knots, sd = 0.1),
      log_lambda_gp_02_level_alpha = if (enable_ms_02) rep(1.0, n_levels) else numeric(0),
      log_lambda_gp_02_level_rho = if (enable_ms_02) rep(1.4, n_levels) else numeric(0),
      log_lambda_gp_02_level_intercept_sd = if (enable_ms_02 && any_re_level_02) rep(0.1, n_levels) else numeric(0),
      log_lambda_gp_02_level_eta = if (enable_ms_02 && n_gp_groups_ms_baseline_02 > 0) {
        matrix(rnorm(n_gp_groups_ms_baseline_02 * n_ms_gp_cal_knots, sd = 0.1),
               nrow = n_gp_groups_ms_baseline_02, ncol = n_ms_gp_cal_knots)
      },
      raw_log_lambda_gp_02_level_intercept = if (n_raw_groups_ms_baseline_02 > 0) rep(0, n_raw_groups_ms_baseline_02),
      cp_log_lambda_gp_02_level_intercept = if (n_cp_groups_ms_baseline_02 > 0) rep(0, n_cp_groups_ms_baseline_02),

      # --- 1→2 Sojourn GP ---
      log_lambda_gp_12_s_pop_intercept = if (need_12_s_gp) array(-4.5, dim = 1),
      log_lambda_gp_12_s_pop_alpha = if (need_12_s_gp) array(1.0, dim = 1),
      log_lambda_gp_12_s_pop_rho = if (need_12_s_gp) array(1.4, dim = 1),
      log_lambda_gp_12_s_pop_eta = if (need_12_s_gp) rnorm(n_ms_gp_sojourn_knots, sd = 0.1),
      log_lambda_gp_12_s_level_alpha = if (need_12_s_gp) rep(1.0, n_levels) else numeric(0),
      log_lambda_gp_12_s_level_rho = if (need_12_s_gp) rep(1.4, n_levels) else numeric(0),
      log_lambda_gp_12_s_level_intercept_sd = if (need_12_s_gp && any_re_level_12_s) rep(0.1, n_levels) else numeric(0),
      log_lambda_gp_12_s_level_eta = if (need_12_s_gp && n_gp_groups_ms_baseline_12_s > 0) {
        matrix(rnorm(n_gp_groups_ms_baseline_12_s * n_ms_gp_sojourn_knots, sd = 0.1),
               nrow = n_gp_groups_ms_baseline_12_s, ncol = n_ms_gp_sojourn_knots)
      },
      raw_log_lambda_gp_12_s_level_intercept = if (n_raw_groups_ms_baseline_12_s > 0) rep(0, n_raw_groups_ms_baseline_12_s),
      cp_log_lambda_gp_12_s_level_intercept = if (n_cp_groups_ms_baseline_12_s > 0) rep(0, n_cp_groups_ms_baseline_12_s),

      # --- 1→2 Clock-forward GP ---
      log_lambda_gp_12_t_pop_intercept = if (need_12_t_gp) array(-4.5, dim = 1),
      log_lambda_gp_12_t_pop_alpha = if (need_12_t_gp) array(1.0, dim = 1),
      log_lambda_gp_12_t_pop_rho = if (need_12_t_gp) array(1.4, dim = 1),
      log_lambda_gp_12_t_pop_eta = if (need_12_t_gp) rnorm(n_ms_gp_cal_knots, sd = 0.1),
      log_lambda_gp_12_t_level_alpha = if (need_12_t_gp) rep(1.0, n_levels) else numeric(0),
      log_lambda_gp_12_t_level_rho = if (need_12_t_gp) rep(1.4, n_levels) else numeric(0),
      log_lambda_gp_12_t_level_intercept_sd = if (need_12_t_gp && any_re_level_12_t) rep(0.1, n_levels) else numeric(0),
      log_lambda_gp_12_t_level_eta = if (need_12_t_gp && n_gp_groups_ms_baseline_12_t > 0) {
        matrix(rnorm(n_gp_groups_ms_baseline_12_t * n_ms_gp_cal_knots, sd = 0.1),
               nrow = n_gp_groups_ms_baseline_12_t, ncol = n_ms_gp_cal_knots)
      },
      raw_log_lambda_gp_12_t_level_intercept = if (n_raw_groups_ms_baseline_12_t > 0) rep(0, n_raw_groups_ms_baseline_12_t),
      cp_log_lambda_gp_12_t_level_intercept = if (n_cp_groups_ms_baseline_12_t > 0) rep(0, n_cp_groups_ms_baseline_12_t),

      # --- 0→3 Transition (Dropout: GP baseline hazard) ---
      log_lambda_gp_03_pop_intercept = if (enable_ms_03) array(-4.5, dim = 1),
      log_lambda_gp_03_pop_alpha = if (enable_ms_03) array(1.0, dim = 1),
      log_lambda_gp_03_pop_rho = if (enable_ms_03) array(1.4, dim = 1),
      log_lambda_gp_03_pop_eta = if (enable_ms_03) rnorm(n_ms_gp_cal_knots, sd = 0.1),
      log_lambda_gp_03_level_alpha = if (enable_ms_03) rep(1.0, n_levels) else numeric(0),
      log_lambda_gp_03_level_rho = if (enable_ms_03) rep(1.4, n_levels) else numeric(0),
      log_lambda_gp_03_level_intercept_sd = if (enable_ms_03 && any_re_level_03) rep(0.1, n_levels) else numeric(0),
      log_lambda_gp_03_level_eta = if (enable_ms_03 && n_gp_groups_ms_baseline_03 > 0) {
        matrix(rnorm(n_gp_groups_ms_baseline_03 * n_ms_gp_cal_knots, sd = 0.1),
               nrow = n_gp_groups_ms_baseline_03, ncol = n_ms_gp_cal_knots)
      },
      raw_log_lambda_gp_03_level_intercept = if (n_raw_groups_ms_baseline_03 > 0) rep(0, n_raw_groups_ms_baseline_03),
      cp_log_lambda_gp_03_level_intercept = if (n_cp_groups_ms_baseline_03 > 0) rep(0, n_cp_groups_ms_baseline_03),

      # --- 3→2 Transition (Off-trial death: sojourn time GP) ---
      log_lambda_gp_32_s_pop_intercept = if (enable_ms_32) array(-4.5, dim = 1),
      log_lambda_gp_32_s_pop_alpha = if (enable_ms_32) array(1.0, dim = 1),
      log_lambda_gp_32_s_pop_rho = if (enable_ms_32) array(1.4, dim = 1),
      log_lambda_gp_32_s_pop_eta = if (enable_ms_32) rnorm(n_ms_gp_sojourn_32_knots, sd = 0.1),
      log_lambda_gp_32_s_level_alpha = if (enable_ms_32) rep(1.0, n_levels) else numeric(0),
      log_lambda_gp_32_s_level_rho = if (enable_ms_32) rep(1.4, n_levels) else numeric(0),
      log_lambda_gp_32_s_level_intercept_sd = if (enable_ms_32 && any_re_level_32) rep(0.1, n_levels) else numeric(0),
      log_lambda_gp_32_s_level_eta = if (enable_ms_32 && n_gp_groups_ms_baseline_32 > 0) {
        matrix(rnorm(n_gp_groups_ms_baseline_32 * n_ms_gp_sojourn_32_knots, sd = 0.1),
               nrow = n_gp_groups_ms_baseline_32, ncol = n_ms_gp_sojourn_32_knots)
      },
      raw_log_lambda_gp_32_s_level_intercept = if (n_raw_groups_ms_baseline_32 > 0) rep(0, n_raw_groups_ms_baseline_32),
      cp_log_lambda_gp_32_s_level_intercept = if (n_cp_groups_ms_baseline_32 > 0) rep(0, n_cp_groups_ms_baseline_32),

      # --- Time-varying covariate coefficients ---
      # 0->1 dimensioning matches stan/modules/multistate/parameters.stan:
      # single coefficient only in observed visit-gated mode; full
      # n_time_varying_covar in continuous and latent visit-gated modes.
      time_varying_coef_01 = {
        observed_vg <- isTRUE(enable_ms_visit_gated_01 == 1L) &&
          !isTRUE(enable_ms_visit_gated_latent_01 == 1L)
        n_tv_01 <- if (observed_vg) 1L else n_time_varying_covar
        if (enable_ms_01 && enable_ms_pop_time_varying_cov && n_tv_01 > 0)
          as.array(rep(0, n_tv_01))
        else NULL
      },
      time_varying_coef_02 = {
        n_tv_02 <- if (isTRUE(enable_ms_02_time_varying_cov == 1L)) n_time_varying_covar else 0L
        if (enable_ms_02 && enable_ms_pop_time_varying_cov && n_tv_02 > 0)
          as.array(rep(0, n_tv_02))
        else NULL
      },
      time_varying_coef_12 = if (enable_ms_12 && enable_ms_pop_time_varying_cov && n_time_varying_covar > 0) {
        as.array(rep(0, n_time_varying_covar))
      },
      time_varying_coef_03 = {
        n_tv_03 <- if (isTRUE(enable_ms_03_time_varying_cov == 1L)) n_time_varying_covar else 0L
        if (enable_ms_03 && enable_ms_pop_time_varying_cov && n_tv_03 > 0)
          as.array(rep(0, n_tv_03))
        else NULL
      },

      # --- Time-invariant covariate coefficients ---
      time_invariant_coef_qr_01 = if (enable_ms_01 && enable_ms_pop_time_invariant_cov && n_time_invariant_covar > 0) {
        rep(0, n_time_invariant_covar)
      },
      time_invariant_coef_qr_02 = if (enable_ms_02 && enable_ms_pop_time_invariant_cov && n_time_invariant_covar > 0) {
        rep(0, n_time_invariant_covar)
      },
      time_invariant_coef_qr_12 = if (enable_ms_12 && enable_ms_pop_time_invariant_cov && n_time_invariant_covar > 0) {
        rep(0, n_time_invariant_covar)
      },
      time_invariant_coef_qr_03 = if (enable_ms_03 && enable_ms_pop_time_invariant_cov &&
                                       isTRUE(enable_ms_03_time_invariant_cov == 1L) &&
                                       n_time_invariant_covar > 0) {
        rep(0, n_time_invariant_covar)
      },
      time_invariant_coef_qr_32 = if (enable_ms_32 && enable_ms_pop_time_invariant_cov &&
                                       isTRUE(enable_ms_32_time_invariant_cov == 1L) &&
                                       n_time_invariant_covar > 0) {
        rep(0, n_time_invariant_covar)
      },

      # --- Multi-level random slopes ---
      sd_level_slope_01 = if (enable_ms_01 && n_time_invariant_covar > 0) {
        lapply(seq_len(n_levels), function(lv) rep(0.1, n_time_invariant_covar))
      },
      sd_level_slope_02 = if (enable_ms_02 && n_time_invariant_covar > 0) {
        lapply(seq_len(n_levels), function(lv) rep(0.1, n_time_invariant_covar))
      },
      sd_level_slope_12 = if (enable_ms_12 && n_time_invariant_covar > 0) {
        lapply(seq_len(n_levels), function(lv) rep(0.1, n_time_invariant_covar))
      },
      sd_level_slope_03 = if (enable_ms_03 && isTRUE(enable_ms_03_time_invariant_cov == 1L) && n_time_invariant_covar > 0) {
        lapply(seq_len(n_levels), function(lv) rep(0.1, n_time_invariant_covar))
      },
      sd_level_slope_32 = if (enable_ms_32 && isTRUE(enable_ms_32_time_invariant_cov == 1L) && n_time_invariant_covar > 0) {
        lapply(seq_len(n_levels), function(lv) rep(0.1, n_time_invariant_covar))
      },
      raw_level_slope_01 = if (enable_ms_01 && n_time_invariant_covar > 0) {
        matrix(0, nrow = n_enabled_groups_ms_slope, ncol = n_time_invariant_covar)
      },
      raw_level_slope_02 = if (enable_ms_02 && n_time_invariant_covar > 0) {
        matrix(0, nrow = n_enabled_groups_ms_slope, ncol = n_time_invariant_covar)
      },
      raw_level_slope_12 = if (enable_ms_12 && n_time_invariant_covar > 0) {
        matrix(0, nrow = n_enabled_groups_ms_slope, ncol = n_time_invariant_covar)
      },
      raw_level_slope_03 = if (enable_ms_03 && isTRUE(enable_ms_03_time_invariant_cov == 1L) && n_time_invariant_covar > 0) {
        matrix(0, nrow = n_enabled_groups_ms_slope, ncol = n_time_invariant_covar)
      },
      raw_level_slope_32 = if (enable_ms_32 && isTRUE(enable_ms_32_time_invariant_cov == 1L) && n_time_invariant_covar > 0) {
        matrix(0, nrow = n_enabled_groups_ms_slope, ncol = n_time_invariant_covar)
      },
    ), .ms_corr_inits)
  })
}

# =========================================================================
# Joint model initializer (tumor + multistate)
# =========================================================================

create_tumor_ssls_initializer_fixed <- function(stan_data, save_dir = NULL, run_id = NULL) {
  # Capture helper in closure so it's available when called in a different process
  .ms_init <- ms_init_values_fixed
  function(chain_id) {
    # Compute max_t_width (same as in Stan's transformed data)
    min_all_t <- min(stan_data$t_patient_visits)
    max_all_t <- max(max(stan_data$t_patient_visits) + 1, stan_data$extend_max_all_t %||% 0)
    max_t_width <- max_all_t - min_all_t + 1

    with(stan_data, {
      # Mirrors Stan's n_forecast_groups_per_level: patient level uses forecast
      # count only — background patients have no explicit NCP parameters.
      n_fgpl_outer <- n_groups_per_level
      if (!is.null(n_forecast_patients)) n_fgpl_outer[n_levels] <- n_forecast_patients

      # Split raw/cp counts matching Stan's routing by level mode:
      # - raw: modes 1 (FE), 2 (RE), 3 (RE+GP) — NCP parameterisation
      # - cp:  mode 4 (RE_CP) — centred parameterisation on natural scale
      n_raw_groups_tr_intercept  <- sum(n_fgpl_outer[enable_level_intercept_tr %in% c(1L, 2L, 3L)])
      n_cp_groups_tr_intercept   <- sum(n_fgpl_outer[enable_level_intercept_tr == 4L])
      n_raw_groups_frac_intercept <- sum(n_fgpl_outer[enable_level_intercept_frac %in% c(1L, 2L, 3L)])
      n_cp_groups_frac_intercept  <- sum(n_fgpl_outer[enable_level_intercept_frac == 4L])
      n_raw_groups_init_intercept <- sum(n_fgpl_outer[enable_level_intercept_init %in% c(1L, 2L, 3L)])
      n_cp_groups_init_intercept  <- sum(n_fgpl_outer[enable_level_intercept_init == 4L])

      # Slope group counts split by mode
      n_raw_groups_tr_slope  <- sum(n_fgpl_outer[(enable_level_intercept_tr %in% c(1L, 2L, 3L)) & (enable_level_cov_tr == 1L)])
      n_cp_groups_tr_slope   <- sum(n_fgpl_outer[(enable_level_intercept_tr == 4L) & (enable_level_cov_tr == 1L)])
      n_raw_groups_frac_slope <- sum(n_fgpl_outer[(enable_level_intercept_frac %in% c(1L, 2L, 3L)) & (enable_level_cov_frac == 1L)])
      n_cp_groups_frac_slope  <- sum(n_fgpl_outer[(enable_level_intercept_frac == 4L) & (enable_level_cov_frac == 1L)])
      n_raw_groups_init_slope <- sum(n_fgpl_outer[(enable_level_intercept_init %in% c(1L, 2L, 3L)) & (enable_level_cov_init == 1L)])
      n_cp_groups_init_slope  <- sum(n_fgpl_outer[(enable_level_intercept_init == 4L) & (enable_level_cov_init == 1L)])

      # gr_decay level group counts (0 when the module/levels are off). The flag
      # arrays are only referenced inside the enable_gr_decay branch so an OFF
      # model (or one where Task 9 wiring is absent) defaults to all-NONE => 0.
      gr_decay_int_modes <- if (isTRUE(enable_gr_decay == 1L)) enable_level_intercept_gr_decay else rep(0L, n_levels)
      gr_decay_cov_modes <- if (isTRUE(enable_gr_decay == 1L)) enable_level_cov_gr_decay else rep(0L, n_levels)
      n_raw_groups_gr_decay_intercept <- sum(n_groups_per_level[gr_decay_int_modes %in% c(1L, 2L, 3L)])
      n_cp_groups_gr_decay_intercept  <- sum(n_groups_per_level[gr_decay_int_modes == 4L])
      n_raw_groups_gr_decay_slope <- sum(n_groups_per_level[(gr_decay_int_modes %in% c(1L, 2L, 3L)) & (gr_decay_cov_modes == 1L)])
      n_cp_groups_gr_decay_slope  <- sum(n_groups_per_level[(gr_decay_int_modes == 4L) & (gr_decay_cov_modes == 1L)])

      # Level positions for indexing flattened arrays
      level_pos <- c(1L, cumsum(n_fgpl_outer) + 1L)

      # SD initial values per level (n_levels elements) — used for RE SD parameter init
      tr_sd_level   <- rep(0.40, n_levels)
      frac_sd_level <- rep(0.40, n_levels)
      init_sd_level <- rep(0.50, n_levels)

      # Generate raw NCP values for raw-routed levels (modes 1/2/3, not mode 4).
      tr_raw_level   <- rnorm(n_raw_groups_tr_intercept,   sd = 0.2)
      frac_raw_level <- rnorm(n_raw_groups_frac_intercept, sd = 0.2)
      init_raw_level <- rnorm(n_raw_groups_init_intercept, sd = 0.2)

      # Centred draws for RE_CP levels (mode 4): sampled on the natural scale ~N(0, sd).
      # rnorm(0, ...) returns numeric(0) when no level uses mode 4 — Stan accepts zero-length arrays.
      tr_cp_level_intercept   <- rnorm(n_cp_groups_tr_intercept,   0, tr_sd_level[enable_level_intercept_tr   == 4L] * 0.2)
      frac_cp_level_intercept <- rnorm(n_cp_groups_frac_intercept, 0, frac_sd_level[enable_level_intercept_frac == 4L] * 0.2)
      init_cp_level_intercept <- rnorm(n_cp_groups_init_intercept, 0, init_sd_level[enable_level_intercept_init == 4L] * 0.2)

      # gr_decay level RE/RE_CP intercept inits (mirrors frac: per-level SD init,
      # NCP raw draws for modes 1/2/3, centred draws for mode 4 on the natural scale).
      gr_decay_sd_level <- rep(0.30, n_levels)
      gr_decay_raw_level <- rnorm(n_raw_groups_gr_decay_intercept, sd = 0.2)
      gr_decay_cp_level_intercept_init <- rnorm(n_cp_groups_gr_decay_intercept, 0, gr_decay_sd_level[gr_decay_int_modes == 4L] * 0.2)

      # Tumor dynamics parameters
      tumor_init <- tibble::lst(
        # Population-level parameters - initialize at prior means
        tr_loc_pop = -2.0,
        frac_logit_loc_pop = 1.5,
        init_logit_loc_pop = 0.0,
        # Gompertz decay intercept: drawn near the prior mean when enabled, length-0
        # when off. Prior defaults mirror r/priors.R (log(0.02), sd 0.75) — priors are
        # not in scope inside the initializer (with(stan_data) only).
        gr_decay_log_loc_pop = if (isTRUE(enable_gr_decay == 1L)) {
          as.array(rnorm(1, log(0.02), 0.75))
        } else {
          numeric(0)
        },
        gr_decay_coef_qr_pop = if (isTRUE(enable_gr_decay == 1L) && isTRUE(enable_pop_cov_gr_decay == 1L) && n_covar > 0) {
          as.array(rep(0, n_covar))
        } else {
          numeric(0)
        },

        # gr_decay level hierarchy (mirrors frac). Gated on enable_gr_decay so an
        # OFF model gets zero-length / zero-size values matching Stan's collapsed
        # param sizing. sd_level_slope is ALWAYS length n_levels (Stan declares it
        # array[n_levels] vector[n_covar]).
        gr_decay_sd_level_intercept_raw = if (isTRUE(enable_gr_decay == 1L)) {
          as.array(gr_decay_sd_level[gr_decay_int_modes %in% c(2L, 4L)])
        } else {
          numeric(0)
        },
        gr_decay_raw_level_intercept = if (isTRUE(enable_gr_decay == 1L)) gr_decay_raw_level else numeric(0),
        gr_decay_cp_level_intercept = if (isTRUE(enable_gr_decay == 1L)) gr_decay_cp_level_intercept_init else numeric(0),
        gr_decay_sd_level_slope = if (isTRUE(enable_gr_decay == 1L) && n_covar > 0) {
          lapply(seq_len(n_levels), function(lv) rep(if (lv == n_levels) 0.03 else 0.05, n_covar))
        } else {
          lapply(seq_len(n_levels), function(lv) numeric(0))
        },
        gr_decay_raw_level_slope = if (isTRUE(enable_gr_decay == 1L) && n_covar > 0) matrix(0, n_raw_groups_gr_decay_slope, n_covar) else matrix(0, 0, max(n_covar, 0)),
        gr_decay_cp_level_slope = if (isTRUE(enable_gr_decay == 1L) && n_covar > 0) matrix(0, n_cp_groups_gr_decay_slope, n_covar) else matrix(0, 0, max(n_covar, 0)),

        # Level-indexed intercept SDs (RE/RE_CP levels: modes 2 and 4), raw NCP values,
        # and centred values for RE_CP levels (mode 4).
        tr_sd_level_intercept_raw = as.array(tr_sd_level[enable_level_intercept_tr %in% c(2L, 4L)]),
        tr_raw_level_intercept = tr_raw_level,
        tr_cp_level_intercept = tr_cp_level_intercept,
        frac_sd_level_intercept_raw = as.array(frac_sd_level[enable_level_intercept_frac %in% c(2L, 4L)]),
        frac_raw_level_intercept = frac_raw_level,
        frac_cp_level_intercept = frac_cp_level_intercept,
        init_sd_level_intercept_raw = as.array(init_sd_level[enable_level_intercept_init %in% c(2L, 4L)]),
        init_raw_level_intercept = init_raw_level,
        init_cp_level_intercept = init_cp_level_intercept,

        # Covariate effects - all zero
        tr_coef_qr_pop = if (n_covar > 0 && enable_pop_cov_tr) array(rep(0, n_covar), dim = n_covar),
        frac_coef_qr_pop = if (n_covar > 0 && enable_pop_cov_frac) array(rep(0, n_covar), dim = n_covar),
        init_coef_qr_pop = if (n_covar > 0 && enable_pop_cov_init) array(rep(0, n_covar), dim = n_covar),

        # Level-indexed slope SDs, raw effects (NCP for modes 1/2/3), and cp effects (mode 4)
        tr_sd_level_slope = if (n_covar > 0) lapply(seq_len(n_levels), function(lv) {
          rep(if (lv == n_levels) 0.2 else 0.1, n_covar)
        }),
        tr_raw_level_slope = if (n_covar > 0) matrix(0, n_raw_groups_tr_slope, n_covar),
        tr_cp_level_slope = if (n_covar > 0) matrix(0, n_cp_groups_tr_slope, n_covar),
        frac_sd_level_slope = if (n_covar > 0) lapply(seq_len(n_levels), function(lv) {
          rep(if (lv == n_levels) 0.2 else 0.1, n_covar)
        }),
        frac_raw_level_slope = if (n_covar > 0) matrix(0, n_raw_groups_frac_slope, n_covar),
        frac_cp_level_slope = if (n_covar > 0) matrix(0, n_cp_groups_frac_slope, n_covar),
        init_sd_level_slope = if (n_covar > 0) lapply(seq_len(n_levels), function(lv) {
          rep(if (lv == n_levels) 0.2 else 0.1, n_covar)
        }),
        init_raw_level_slope = if (n_covar > 0) matrix(0, n_raw_groups_init_slope, n_covar),
        init_cp_level_slope = if (n_covar > 0) matrix(0, n_cp_groups_init_slope, n_covar),

        # Patient-level process noise (all disabled for now)
        tr_raw_patient_process_noise = if (enable_patient_process_noise_tr) {
          matrix(0, nrow = n_patients, ncol = max_t_width)
        },
        tr_log_sd_pop_process_noise = if (enable_patient_process_noise_tr) array(log(0.05)),
        tr_sd_patient_log_sd_process_noise = if (enable_patient_process_noise_tr) array(0.1),
        tr_raw_patient_log_sd_process_noise = if (enable_patient_process_noise_sd_tr) rep(0, n_patients),
        tr_logit_phi_pop_process_noise = if (enable_patient_process_noise_tr) array(2),
        tr_sd_patient_phi_process_noise = if (enable_patient_process_noise_tr) array(0.05),
        tr_raw_patient_phi_process_noise = if (enable_patient_process_noise_phi_tr) rep(0, n_patients),

        # Population-level process noise
        tr_raw_pop_process_noise = if (enable_pop_process_noise_tr) {
          rep(0, max_t_width)
        },
        tr_log_sd_pop_process_noise_pop = if (enable_pop_process_noise_tr) array(log(0.05)),
        tr_logit_phi_pop_process_noise_pop = if (enable_pop_process_noise_tr) array(2),

        # Measurement error - realistic value
        measure_sd_sld = 0.15,
      )

      # Combine tumor + multistate init values
      c(tumor_init, .ms_init(environment()))
    }) |> purrr::compact() |>
      (\(init_vals) {
        # Save init values to disk if save_dir is provided
        if (!is.null(save_dir)) {
          dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
          filename <- if (!is.null(run_id)) {
            sprintf("init-%s-chain-%d.json", run_id, chain_id)
          } else {
            sprintf("init-chain-%d.json", chain_id)
          }
          init_file <- file.path(save_dir, filename)
          jsonlite::write_json(init_vals, init_file, auto_unbox = TRUE, pretty = TRUE)
          message("Saved init values for chain ", chain_id, " to ", init_file)
        }
        init_vals
      })()
  }
}

# =========================================================================
# Pathfinder pre-solve for tumor SSLS (data-only warm-start bootstrap)
# =========================================================================

# Allow-list of POPULATION SCALAR base names seeded from a Pathfinder pre-solve.
# Shared by run_tumor_ssls_pathfinder (to select just these columns via select_draws)
# and create_tumor_ssls_pathfinder_initializer_fixed (to extract the seed values),
# so the two can never drift. GP pop intercept/alpha/rho are declared array[1] in
# Stan -> Pathfinder names them "<p>[1]"; the plain SLD locations are bare scalars.
# All transition slots are listed; disabled ones simply won't appear in the output
# (select_draws matches against stan_variables base names; absent names don't match).
tumor_ssls_pathfinder_seed_params <- function() {
  gp_slots <- c("01", "02", "03", "12_s", "12_t", "32_s")
  c(
    "tr_loc_pop", "frac_logit_loc_pop", "init_logit_loc_pop", "measure_sd_sld",
    as.vector(outer(
      sprintf("log_lambda_gp_%s_pop", gp_slots),
      c("intercept", "alpha", "rho"),
      paste, sep = "_"
    ))
  )
}

#' Run a Pathfinder pre-solve on the tumor SSLS model.
#'
#' A ~12s variational L-BFGS pass over the SAME data the posterior fit uses,
#' purely to locate the typical set. Its draws seed the population scalars of the
#' sampling init (see create_tumor_ssls_pathfinder_initializer_fixed), which
#' eliminates the cold-start inner-Laplace-solver thrash that made an unseeded
#' run intractable (~80 min/warmup-iter -> ~5 s/iter). Requires NOTHING but the
#' data — works on a first run, no prior MCMC fit needed.
#'
#' @param exe_file Path to the compiled sf-ssm-log-space executable.
#' @param stan_data Assembled Stan data list (surrogate-enabled posterior data).
#' @param seed Pathfinder RNG seed.
#' @return A draws_df of the Pathfinder draws (NOT the live fit object — the
#'   draws are extracted here so the value is qs2-serializable across targets/
#'   crew workers; the fit's tempfile CSVs would not survive serialization).
run_tumor_ssls_pathfinder <- function(exe_file, stan_data, seed = 123) {
  # Mirror sample_and_save()'s exe handling: make executable, set stan_threads
  # via the private field (cmdstanr#765 workaround — cpp_options ignored for
  # exe_file-loaded models).
  if (fs::file_exists(exe_file) && !fs::file_access(exe_file, "execute")) {
    fs::file_chmod(exe_file, "u+x")
  }
  model <- cmdstanr::cmdstan_model(exe_file = exe_file)
  model$.__enclos_env__$private$cpp_options_$stan_threads <- TRUE

  n_threads <- stan_data$n_shards %||% 1L
  fixed_init <- create_tumor_ssls_initializer_fixed(stan_data)

  pf <- model$pathfinder(
    data = stan_data,
    init = list(fixed_init(1)),
    num_paths = 1,
    single_path_draws = 40,
    draws = 40,
    max_lbfgs_iters = 30,
    history_size = 6,
    num_threads = n_threads,
    seed = seed,
    refresh = 0
  )
  # Return PLAIN draws (qs2-safe), not the tempfile-backed fit object — and only
  # the ~22 population-scalar columns the seeded initializer consumes. The full
  # all-497 draws object is ~1.8M columns (per-patient/per-visit GQ); reading just
  # the seed params keeps this cheap. Uses the project select_draws() convention,
  # which is now fit-type-agnostic (S3 cmdstan_draws_field dispatch) so it works on
  # a CmdStanPathfinder fit. any_of() takes the exact base-name vector by value and
  # silently drops names absent from the output (disabled transition slots) — unlike
  # matches() with a runtime regex variable, which select_draws' substitute(c(...))
  # NSE capture cannot resolve in a crew worker (fails with extent=0).
  posterior::as_draws_df(select_draws(pf, any_of(tumor_ssls_pathfinder_seed_params())))
}

# =========================================================================
# Pathfinder-seeded tumor SSLS initializer (population scalars only)
# =========================================================================

#' Tumor SSLS initializer seeded from a Pathfinder pre-solve.
#'
#' Reuses create_tumor_ssls_initializer_fixed() for the structurally-correct
#' skeleton (all NCP routing, conditional zero-length slots, per-level/patient/
#' frailty params), then OVERRIDES only the well-identified POPULATION SCALARS
#' with a draw from a Pathfinder fit on the same data. This puts the chains in
#' the typical set (killing cold-start trajectory blowup + the inner-Laplace
#' solver thrash) while leaving every higher-dim / conditional slot at the
#' known-valid fixed values. Works on a first run: Pathfinder is a ~12s
#' data-only bootstrap, no prior MCMC fit required.
#'
#' Scope is deliberately narrow (population locations + GP pop intercept/alpha/
#' rho scalars). Pathfinder's joint approximation is rough (high Pareto k), so
#' seeding only the low-dim params it estimates well is the low-risk choice.
#'
#' @param pathfinder_draws A draws_df (or coercible data frame) of Pathfinder
#'   draws on the same stan_data. Taking a PLAIN data frame rather than the live
#'   CmdStanPathfinder object is deliberate: a fit object holds tempfile-backed
#'   CSV paths that do not survive qs2 serialization across crew workers, so the
#'   target extracts draws once (run_tumor_ssls_pathfinder) and passes them here.
#' @param stan_data The assembled Stan data list.
#' @param save_dir,run_id Passed through to the base initializer for init JSON dumps.
#' @return An init function of chain_id, suitable for cmdstanr $sample(init=).
create_tumor_ssls_pathfinder_initializer_fixed <- function(
    pathfinder_draws, stan_data, save_dir = NULL, run_id = NULL) {
  base_init_fn <- create_tumor_ssls_initializer_fixed(stan_data, save_dir, run_id)

  draws_df <- posterior::as_draws_df(pathfinder_draws)
  pf_names <- setdiff(colnames(draws_df), c(".chain", ".iteration", ".draw"))

  # Allow-list of POPULATION SCALAR params to seed from Pathfinder (shared with
  # run_tumor_ssls_pathfinder, which subsets $draws() to exactly these). Disabled
  # transition slots are NULL in the base skeleton and skipped in the loop below.
  override_params <- tumor_ssls_pathfinder_seed_params()

  # Pull a scalar from the Pathfinder draw row by param name, trying the bare
  # name then the array-element "<p>[1]" form. Returns NULL if neither exists.
  pf_scalar <- function(draw_row, p) {
    if (p %in% pf_names) return(as.numeric(draw_row[[p]]))
    pe <- paste0(p, "[1]")
    if (pe %in% pf_names) return(as.numeric(draw_row[[pe]]))
    NULL
  }

  function(chain_id) {
    init_vals <- base_init_fn(chain_id)
    draw_row <- dplyr::slice_sample(draws_df, n = 1)

    for (p in override_params) {
      base_val <- init_vals[[p]]
      if (is.null(base_val)) next            # disabled slot -> keep skeleton's NULL
      pf_val <- pf_scalar(draw_row, p)
      if (is.null(pf_val) || !is.finite(pf_val)) next
      # Preserve array-ness: Stan rejects a plain scalar where it expects array[1].
      init_vals[[p]] <- if (is.array(base_val)) array(pf_val, dim = dim(base_val)) else pf_val
    }
    init_vals
  }
}

# =========================================================================
# Standalone multistate model initializer
# =========================================================================

create_ms_standalone_initializer_fixed <- function(stan_data, save_dir = NULL, run_id = NULL) {
  # Capture helper in closure so it's available when called in a different process
  .ms_init <- ms_init_values_fixed
  function(chain_id) {
    with(stan_data, {
      .ms_init(environment())
    }) |> purrr::compact() |>
      (\(init_vals) {
        if (!is.null(save_dir)) {
          dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
          filename <- if (!is.null(run_id)) {
            sprintf("init-%s-chain-%d.json", run_id, chain_id)
          } else {
            sprintf("init-chain-%d.json", chain_id)
          }
          init_file <- file.path(save_dir, filename)
          jsonlite::write_json(init_vals, init_file, auto_unbox = TRUE, pretty = TRUE)
          message("Saved init values for chain ", chain_id, " to ", init_file)
        }
        init_vals
      })()
  }
}
