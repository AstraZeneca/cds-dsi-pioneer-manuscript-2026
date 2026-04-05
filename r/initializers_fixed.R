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

    # Multistate enabled group counts - truthy (> 0) for both intercept-only and GP modes
    n_enabled_groups_ms_baseline <- sum(n_groups_per_level[enable_ms_level_baseline_hazard > 0])
    n_enabled_groups_ms_baseline_01 <- if (enable_ms_01) n_enabled_groups_ms_baseline else 0L
    n_enabled_groups_ms_baseline_02 <- if (enable_ms_02) n_enabled_groups_ms_baseline else 0L
    n_enabled_groups_ms_baseline_12_s <- if (need_12_s_gp) n_enabled_groups_ms_baseline else 0L
    n_enabled_groups_ms_baseline_12_t <- if (need_12_t_gp) n_enabled_groups_ms_baseline else 0L
    n_enabled_groups_ms_slope <- sum(n_groups_per_level[enable_ms_level_cov == 1])

    # GP-only group counts (mode == 3, for eta matrix sizing)
    n_gp_groups_ms_baseline <- sum(n_groups_per_level[enable_ms_level_baseline_hazard == 3L])
    n_gp_groups_ms_baseline_01 <- if (enable_ms_01) n_gp_groups_ms_baseline else 0L
    n_gp_groups_ms_baseline_02 <- if (enable_ms_02) n_gp_groups_ms_baseline else 0L
    n_gp_groups_ms_baseline_12_s <- if (need_12_s_gp) n_gp_groups_ms_baseline else 0L
    n_gp_groups_ms_baseline_12_t <- if (need_12_t_gp) n_gp_groups_ms_baseline else 0L
    n_gp_groups_ms_baseline_03 <- if (enable_ms_03) n_gp_groups_ms_baseline else 0L
    n_gp_groups_ms_baseline_32 <- if (enable_ms_32) n_gp_groups_ms_baseline else 0L
    any_re_level <- any(enable_ms_level_baseline_hazard >= 2L)

    # GP knot counts (matches Stan transformed_data ceiling division)
    n_ms_gp_cal_knots     <- ceiling(max_all_t / ms_gp_grid_step)
    n_ms_gp_sojourn_knots <- ceiling(ms_max_sojourn_t / ms_gp_grid_step)

    tibble::lst(
      # --- 0→1 Transition ---
      log_lambda_gp_01_pop_intercept = if (enable_ms_01) array(-4.5, dim = 1),
      log_lambda_gp_01_pop_alpha = if (enable_ms_01) array(1.0, dim = 1),
      log_lambda_gp_01_pop_rho = if (enable_ms_01) array(1.4, dim = 1),
      log_lambda_gp_01_pop_eta = if (enable_ms_01) rnorm(n_ms_gp_cal_knots, sd = 0.1),
      log_lambda_gp_01_level_alpha = rep(1.0, n_levels),
      log_lambda_gp_01_level_rho = rep(1.4, n_levels),
      log_lambda_gp_01_level_intercept_sd = if (any_re_level) rep(0.1, n_levels) else numeric(0),
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
      log_lambda_gp_02_level_intercept_sd = if (enable_ms_02 && any_re_level) rep(0.1, n_levels) else numeric(0),
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
      log_lambda_gp_12_s_level_intercept_sd = if (need_12_s_gp && any_re_level) rep(0.1, n_levels) else numeric(0),
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
      log_lambda_gp_12_t_level_intercept_sd = if (need_12_t_gp && any_re_level) rep(0.1, n_levels) else numeric(0),
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
      log_lambda_gp_03_level_intercept_sd = if (enable_ms_03 && any_re_level) rep(0.1, n_levels) else numeric(0),
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
      log_lambda_gp_32_s_level_intercept_sd = if (enable_ms_32 && any_re_level) rep(0.1, n_levels) else numeric(0),
      log_lambda_gp_32_s_level_eta = if (enable_ms_32 && n_gp_groups_ms_baseline_32 > 0) {
        matrix(rnorm(n_gp_groups_ms_baseline_32 * n_ms_gp_sojourn_32_knots, sd = 0.1),
               nrow = n_gp_groups_ms_baseline_32, ncol = n_ms_gp_sojourn_32_knots)
      },
      raw_log_lambda_gp_32_s_level_intercept = if (n_enabled_groups_ms_baseline_32 > 0) rep(0, n_enabled_groups_ms_baseline_32),

      # --- Time-varying covariate coefficients ---
      time_varying_coef_01 = {
        n_tv_01 <- if (isTRUE(enable_ms_visit_gated_01 == 1L)) 1L else n_time_varying_covar
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
      raw_level_slope_01 = if (enable_ms_01 && n_time_invariant_covar > 0) {
        matrix(0, nrow = n_enabled_groups_ms_slope, ncol = n_time_invariant_covar)
      },
      raw_level_slope_02 = if (enable_ms_02 && n_time_invariant_covar > 0) {
        matrix(0, nrow = n_enabled_groups_ms_slope, ncol = n_time_invariant_covar)
      },
      raw_level_slope_12 = if (enable_ms_12 && n_time_invariant_covar > 0) {
        matrix(0, nrow = n_enabled_groups_ms_slope, ncol = n_time_invariant_covar)
      },
    )
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
      # Compute enabled group counts for each module (matches Stan transformed_data)
      n_enabled_groups_tr_intercept <- sum(n_groups_per_level[enable_level_intercept_tr == 1])
      n_enabled_groups_tr_slope <- sum(n_groups_per_level[enable_level_cov_tr == 1])
      n_enabled_groups_frac_intercept <- sum(n_groups_per_level[enable_level_intercept_frac == 1])
      n_enabled_groups_frac_slope <- sum(n_groups_per_level[enable_level_cov_frac == 1])
      n_enabled_groups_init_intercept <- sum(n_groups_per_level[enable_level_intercept_init == 1])
      n_enabled_groups_init_slope <- sum(n_groups_per_level[enable_level_cov_init == 1])

      # Level positions for indexing flattened arrays
      level_pos <- c(1L, cumsum(n_groups_per_level) + 1L)

      # Hierarchical SDs per level - c(trial, patient) for 2-level
      tr_sd_level <- c(0.35, 0.40)
      frac_sd_level <- c(0.35, 0.40)
      init_sd_level <- c(0.35, 0.50)

      # Generate deviations for ENABLED levels only (flattened)
      tr_level_dev <- unlist(lapply(seq_len(n_levels), function(lv) {
        if (enable_level_intercept_tr[lv]) rnorm(n_groups_per_level[lv], sd = if (lv == n_levels) 0.3 else 0.2) else NULL
      }))
      frac_level_dev <- unlist(lapply(seq_len(n_levels), function(lv) {
        if (enable_level_intercept_frac[lv]) rnorm(n_groups_per_level[lv], sd = if (lv == n_levels) 0.3 else 0.2) else NULL
      }))
      init_level_dev <- unlist(lapply(seq_len(n_levels), function(lv) {
        if (enable_level_intercept_init[lv]) rnorm(n_groups_per_level[lv], sd = if (lv == n_levels) 0.3 else 0.2) else NULL
      }))

      # Back-calculate raw values using enabled level positions
      enabled_level_pos_tr <- c(1L, cumsum(n_groups_per_level * enable_level_intercept_tr) + 1L)
      tr_raw_level <- unlist(lapply(seq_len(n_levels), function(lv) {
        if (!enable_level_intercept_tr[lv]) return(NULL)
        idx <- enabled_level_pos_tr[lv]:(enabled_level_pos_tr[lv + 1] - 1)
        tr_level_dev[idx] / tr_sd_level[lv]
      }))

      enabled_level_pos_frac <- c(1L, cumsum(n_groups_per_level * enable_level_intercept_frac) + 1L)
      frac_raw_level <- unlist(lapply(seq_len(n_levels), function(lv) {
        if (!enable_level_intercept_frac[lv]) return(NULL)
        idx <- enabled_level_pos_frac[lv]:(enabled_level_pos_frac[lv + 1] - 1)
        frac_level_dev[idx] / frac_sd_level[lv]
      }))

      enabled_level_pos_init <- c(1L, cumsum(n_groups_per_level * enable_level_intercept_init) + 1L)
      init_raw_level <- unlist(lapply(seq_len(n_levels), function(lv) {
        if (!enable_level_intercept_init[lv]) return(NULL)
        idx <- enabled_level_pos_init[lv]:(enabled_level_pos_init[lv + 1] - 1)
        init_level_dev[idx] / init_sd_level[lv]
      }))

      # Tumor dynamics parameters
      tumor_init <- tibble::lst(
        # Population-level parameters - initialize at prior means
        tr_loc_pop = -2.0,
        frac_logit_loc_pop = 1.5,
        init_logit_loc_pop = 0.0,

        # Level-indexed intercept SDs and raw values
        tr_sd_level_intercept = tr_sd_level,
        tr_raw_level_intercept = tr_raw_level,
        frac_sd_level_intercept = frac_sd_level,
        frac_raw_level_intercept = frac_raw_level,
        init_sd_level_intercept = init_sd_level,
        init_raw_level_intercept = init_raw_level,

        # Covariate effects - all zero
        tr_coef_qr_pop = if (n_covar > 0 && enable_pop_cov_tr) array(rep(0, n_covar), dim = n_covar),
        frac_coef_qr_pop = if (n_covar > 0 && enable_pop_cov_frac) array(rep(0, n_covar), dim = n_covar),
        init_coef_qr_pop = if (n_covar > 0 && enable_pop_cov_init) array(rep(0, n_covar), dim = n_covar),

        # Level-indexed slope SDs and raw effects
        tr_sd_level_slope = if (n_covar > 0) lapply(seq_len(n_levels), function(lv) {
          rep(if (lv == n_levels) 0.2 else 0.1, n_covar)
        }),
        tr_raw_level_slope = if (n_covar > 0) matrix(0, n_enabled_groups_tr_slope, n_covar),
        frac_sd_level_slope = if (n_covar > 0) lapply(seq_len(n_levels), function(lv) {
          rep(if (lv == n_levels) 0.2 else 0.1, n_covar)
        }),
        frac_raw_level_slope = if (n_covar > 0) matrix(0, n_enabled_groups_frac_slope, n_covar),
        init_sd_level_slope = if (n_covar > 0) lapply(seq_len(n_levels), function(lv) {
          rep(if (lv == n_levels) 0.2 else 0.1, n_covar)
        }),
        init_raw_level_slope = if (n_covar > 0) matrix(0, n_enabled_groups_init_slope, n_covar),

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
