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

    n_enabled_slot <- vapply(ms_legacy_mode,
      function(m) sum(n_groups_per_level[m > 0L]), integer(1))
    n_gp_slot <- vapply(ms_legacy_mode,
      function(m) sum(n_groups_per_level[m == 3L]), integer(1))
    any_re_slot <- vapply(ms_legacy_mode,
      function(m) any(m == 2L | m == 3L | m == 4L), logical(1))

    # Enabled group counts - truthy (> 0) for both intercept-only and GP modes
    n_enabled_groups_ms_baseline_01 <- n_enabled_slot[1]
    n_enabled_groups_ms_baseline_02 <- n_enabled_slot[2]
    n_enabled_groups_ms_baseline_03 <- n_enabled_slot[3]
    n_enabled_groups_ms_baseline_12_s <- n_enabled_slot[4]
    n_enabled_groups_ms_baseline_12_t <- n_enabled_slot[5]
    n_enabled_groups_ms_baseline_32 <- n_enabled_slot[6]
    n_enabled_groups_ms_slope <- sum(n_groups_per_level[enable_ms_level_cov == 1])

    # GP-only group counts (legacy mode == 3, for eta matrix sizing)
    n_gp_groups_ms_baseline_01 <- n_gp_slot[1]
    n_gp_groups_ms_baseline_02 <- n_gp_slot[2]
    n_gp_groups_ms_baseline_03 <- n_gp_slot[3]
    n_gp_groups_ms_baseline_12_s <- n_gp_slot[4]
    n_gp_groups_ms_baseline_12_t <- n_gp_slot[5]
    n_gp_groups_ms_baseline_32 <- n_gp_slot[6]

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
    # Mirror Stan's n_forecast_groups_per_level: the patient (last) level uses
    # the forecast-patient count, not the full patient count. Empty list when no
    # block is configured (corr_group all-zero) — keeps the scalar path intact.
    n_forecast_groups_per_level <- n_groups_per_level
    if (exists("n_forecast_patients", inherits = FALSE)) {
      n_forecast_groups_per_level[n_levels] <- n_forecast_patients
    }
    .ms_corr <- ms_corr_blocks(
      ms_level_intercept_corr_group, ms_slot_active, n_forecast_groups_per_level
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
        if (isTRUE(enable_ms_baseline_trend_01 == 1L) && n_enabled_groups_ms_baseline_01 > 0) {
          rep(0, n_enabled_groups_ms_baseline_01)
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
      raw_log_lambda_gp_01_level_intercept = if (n_enabled_groups_ms_baseline_01 > 0) rep(0, n_enabled_groups_ms_baseline_01),

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
      raw_log_lambda_gp_02_level_intercept = if (n_enabled_groups_ms_baseline_02 > 0) rep(0, n_enabled_groups_ms_baseline_02),

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
      raw_log_lambda_gp_12_s_level_intercept = if (n_enabled_groups_ms_baseline_12_s > 0) rep(0, n_enabled_groups_ms_baseline_12_s),

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
      raw_log_lambda_gp_12_t_level_intercept = if (n_enabled_groups_ms_baseline_12_t > 0) rep(0, n_enabled_groups_ms_baseline_12_t),

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
      raw_log_lambda_gp_03_level_intercept = if (n_enabled_groups_ms_baseline_03 > 0) rep(0, n_enabled_groups_ms_baseline_03),

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
      raw_log_lambda_gp_32_s_level_intercept = if (n_enabled_groups_ms_baseline_32 > 0) rep(0, n_enabled_groups_ms_baseline_32),

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
      # Multi-level hierarchy: n_levels, n_groups_per_level
      # Split raw/cp counts matching Stan's routing by level mode:
      # - raw: modes 1 (FE), 2 (RE), 3 (RE+GP) — NCP parameterisation
      # - cp:  mode 4 (RE_CP) — centred parameterisation on natural scale
      n_raw_groups_tr_intercept  <- sum(n_groups_per_level[enable_level_intercept_tr %in% c(1L, 2L, 3L)])
      n_cp_groups_tr_intercept   <- sum(n_groups_per_level[enable_level_intercept_tr == 4L])
      n_raw_groups_frac_intercept <- sum(n_groups_per_level[enable_level_intercept_frac %in% c(1L, 2L, 3L)])
      n_cp_groups_frac_intercept  <- sum(n_groups_per_level[enable_level_intercept_frac == 4L])
      n_raw_groups_init_intercept <- sum(n_groups_per_level[enable_level_intercept_init %in% c(1L, 2L, 3L)])
      n_cp_groups_init_intercept  <- sum(n_groups_per_level[enable_level_intercept_init == 4L])

      # Slope group counts split by mode
      n_raw_groups_tr_slope  <- sum(n_groups_per_level[(enable_level_intercept_tr %in% c(1L, 2L, 3L)) & (enable_level_cov_tr == 1L)])
      n_cp_groups_tr_slope   <- sum(n_groups_per_level[(enable_level_intercept_tr == 4L) & (enable_level_cov_tr == 1L)])
      n_raw_groups_frac_slope <- sum(n_groups_per_level[(enable_level_intercept_frac %in% c(1L, 2L, 3L)) & (enable_level_cov_frac == 1L)])
      n_cp_groups_frac_slope  <- sum(n_groups_per_level[(enable_level_intercept_frac == 4L) & (enable_level_cov_frac == 1L)])
      n_raw_groups_init_slope <- sum(n_groups_per_level[(enable_level_intercept_init %in% c(1L, 2L, 3L)) & (enable_level_cov_init == 1L)])
      n_cp_groups_init_slope  <- sum(n_groups_per_level[(enable_level_intercept_init == 4L) & (enable_level_cov_init == 1L)])

      # Level positions for indexing flattened arrays
      level_pos <- c(1L, cumsum(n_groups_per_level) + 1L)

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

      # Tumor dynamics parameters
      tumor_init <- tibble::lst(
        # Population-level parameters - initialize at prior means
        tr_loc_pop = -2.0,
        frac_logit_loc_pop = 1.5,
        init_logit_loc_pop = 0.0,

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
