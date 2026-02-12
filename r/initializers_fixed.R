# Temporary fixed initializer for testing chain stability
# Uses reproducible per-chain values with realistic heterogeneity

create_tumor_ssls_initializer_fixed <- function(stan_data, save_dir = NULL, run_id = NULL) {
  function(chain_id) {
    # Note: targets already sets a deterministic seed per target, so we don't need set.seed() here
    # The chain_id parameter is kept for interface compatibility but not used for seeding

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
      # Derived flags for 1→2 GPs
      need_12_s_gp <- enable_ms_12 && (ms_time_scale_12 == 1 || ms_time_scale_12 == 2)
      need_12_t_gp <- enable_ms_12 && (ms_time_scale_12 == 0 || ms_time_scale_12 == 2)

      # Multistate enabled group counts - SEPARATE for each transition
      n_enabled_groups_ms_baseline_01 <- if (enable_ms_01) sum(n_groups_per_level[enable_ms_level_baseline_hazard == 1]) else 0L
      n_enabled_groups_ms_baseline_02 <- if (enable_ms_02) sum(n_groups_per_level[enable_ms_level_baseline_hazard == 1]) else 0L
      n_enabled_groups_ms_baseline_12_s <- if (need_12_s_gp) sum(n_groups_per_level[enable_ms_level_baseline_hazard == 1]) else 0L
      n_enabled_groups_ms_baseline_12_t <- if (need_12_t_gp) sum(n_groups_per_level[enable_ms_level_baseline_hazard == 1]) else 0L
      n_enabled_groups_ms_slope <- sum(n_groups_per_level[enable_ms_level_cov == 1])

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

      tibble::lst(
        # Population-level parameters - initialize at prior means to avoid catastrophic random inits
        tr_loc_pop = -2.0,  # Total rate on log scale (prior mean)
        frac_logit_loc_pop = 1.5,  # Fraction on logit scale (prior mean)
        init_logit_loc_pop = 0.0,  # Initial fraction logit (prior mean)

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
        # Slope SDs are always n_levels, but raw effects are sized by enabled groups only
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

        # =====================================================================
        # MULTISTATE HAZARD MODEL - Fixed initializers for testing
        # =====================================================================

        # Derived flags for 1→2 GPs
        need_12_s_gp = enable_ms_12 && (ms_time_scale_12 == 1 || ms_time_scale_12 == 2),
        need_12_t_gp = enable_ms_12 && (ms_time_scale_12 == 0 || ms_time_scale_12 == 2),

        # --- 0→1 Transition ---
        ms_log_lambda_gp_01_pop_intercept = if (enable_ms_01) array(-4.5, dim = 1),
        ms_log_lambda_gp_01_pop_alpha = if (enable_ms_01) array(1.0, dim = 1),
        ms_log_lambda_gp_01_pop_rho = if (enable_ms_01) array(1.4, dim = 1),
        ms_log_lambda_gp_01_pop_eta = if (enable_ms_01) rnorm(max_all_t, sd = 0.1),
        ms_log_lambda_gp_01_level_alpha = rep(1.0, n_levels),
        ms_log_lambda_gp_01_level_rho = rep(1.4, n_levels),
        ms_log_lambda_gp_01_level_intercept_sd = rep(0.1, n_levels),
        ms_log_lambda_gp_01_level_eta = if (enable_ms_01 && n_enabled_groups_ms_baseline_01 > 0) {
          matrix(rnorm(n_enabled_groups_ms_baseline_01 * max_all_t, sd = 0.1),
                 nrow = n_enabled_groups_ms_baseline_01, ncol = max_all_t)
        },
        ms_raw_log_lambda_gp_01_level_intercept = if (n_enabled_groups_ms_baseline_01 > 0) rep(0, n_enabled_groups_ms_baseline_01),

        # --- 0→2 Transition ---
        ms_log_lambda_gp_02_pop_intercept = if (enable_ms_02) array(-4.5, dim = 1),
        ms_log_lambda_gp_02_pop_alpha = if (enable_ms_02) array(1.0, dim = 1),
        ms_log_lambda_gp_02_pop_rho = if (enable_ms_02) array(1.4, dim = 1),
        ms_log_lambda_gp_02_pop_eta = if (enable_ms_02) rnorm(max_all_t, sd = 0.1),
        ms_log_lambda_gp_02_level_alpha = if (enable_ms_02) rep(1.0, n_levels) else numeric(0),
        ms_log_lambda_gp_02_level_rho = if (enable_ms_02) rep(1.4, n_levels) else numeric(0),
        ms_log_lambda_gp_02_level_intercept_sd = if (enable_ms_02) rep(0.1, n_levels) else numeric(0),
        ms_log_lambda_gp_02_level_eta = if (enable_ms_02 && n_enabled_groups_ms_baseline_02 > 0) {
          matrix(rnorm(n_enabled_groups_ms_baseline_02 * max_all_t, sd = 0.1),
                 nrow = n_enabled_groups_ms_baseline_02, ncol = max_all_t)
        },
        ms_raw_log_lambda_gp_02_level_intercept = if (n_enabled_groups_ms_baseline_02 > 0) rep(0, n_enabled_groups_ms_baseline_02),

        # --- 1→2 Sojourn GP ---
        ms_log_lambda_gp_12_s_pop_intercept = if (need_12_s_gp) array(-4.5, dim = 1),
        ms_log_lambda_gp_12_s_pop_alpha = if (need_12_s_gp) array(1.0, dim = 1),
        ms_log_lambda_gp_12_s_pop_rho = if (need_12_s_gp) array(1.4, dim = 1),
        ms_log_lambda_gp_12_s_pop_eta = if (need_12_s_gp) rnorm(ms_max_sojourn_t, sd = 0.1),
        ms_log_lambda_gp_12_s_level_alpha = if (need_12_s_gp) rep(1.0, n_levels) else numeric(0),
        ms_log_lambda_gp_12_s_level_rho = if (need_12_s_gp) rep(1.4, n_levels) else numeric(0),
        ms_log_lambda_gp_12_s_level_intercept_sd = if (need_12_s_gp) rep(0.1, n_levels) else numeric(0),
        ms_log_lambda_gp_12_s_level_eta = if (need_12_s_gp && n_enabled_groups_ms_baseline_12_s > 0) {
          matrix(rnorm(n_enabled_groups_ms_baseline_12_s * ms_max_sojourn_t, sd = 0.1),
                 nrow = n_enabled_groups_ms_baseline_12_s, ncol = ms_max_sojourn_t)
        },
        ms_raw_log_lambda_gp_12_s_level_intercept = if (n_enabled_groups_ms_baseline_12_s > 0) rep(0, n_enabled_groups_ms_baseline_12_s),

        # --- 1→2 Clock-forward GP ---
        ms_log_lambda_gp_12_t_pop_intercept = if (need_12_t_gp) array(-4.5, dim = 1),
        ms_log_lambda_gp_12_t_pop_alpha = if (need_12_t_gp) array(1.0, dim = 1),
        ms_log_lambda_gp_12_t_pop_rho = if (need_12_t_gp) array(1.4, dim = 1),
        ms_log_lambda_gp_12_t_pop_eta = if (need_12_t_gp) rnorm(max_all_t, sd = 0.1),
        ms_log_lambda_gp_12_t_level_alpha = if (need_12_t_gp) rep(1.0, n_levels) else numeric(0),
        ms_log_lambda_gp_12_t_level_rho = if (need_12_t_gp) rep(1.4, n_levels) else numeric(0),
        ms_log_lambda_gp_12_t_level_intercept_sd = if (need_12_t_gp) rep(0.1, n_levels) else numeric(0),
        ms_log_lambda_gp_12_t_level_eta = if (need_12_t_gp && n_enabled_groups_ms_baseline_12_t > 0) {
          matrix(rnorm(n_enabled_groups_ms_baseline_12_t * max_all_t, sd = 0.1),
                 nrow = n_enabled_groups_ms_baseline_12_t, ncol = max_all_t)
        },
        ms_raw_log_lambda_gp_12_t_level_intercept = if (n_enabled_groups_ms_baseline_12_t > 0) rep(0, n_enabled_groups_ms_baseline_12_t),

        # --- Time-varying covariate coefficients ---
        ms_time_varying_coef_01 = if (enable_ms_01 && enable_ms_pop_time_varying_cov && n_time_varying_covar > 0) {
          rep(0, n_time_varying_covar)
        },
        ms_time_varying_coef_02 = if (enable_ms_02 && enable_ms_pop_time_varying_cov && n_time_varying_covar > 0) {
          rep(0, n_time_varying_covar)
        },
        ms_time_varying_coef_12 = if (enable_ms_12 && enable_ms_pop_time_varying_cov && n_time_varying_covar > 0) {
          rep(0, n_time_varying_covar)
        },

        # --- Time-invariant covariate coefficients ---
        ms_time_invariant_coef_qr_01 = if (enable_ms_01 && enable_ms_pop_time_invariant_cov && n_time_invariant_covar > 0) {
          rep(0, n_time_invariant_covar)
        },
        ms_time_invariant_coef_qr_02 = if (enable_ms_02 && enable_ms_pop_time_invariant_cov && n_time_invariant_covar > 0) {
          rep(0, n_time_invariant_covar)
        },
        ms_time_invariant_coef_qr_12 = if (enable_ms_12 && enable_ms_pop_time_invariant_cov && n_time_invariant_covar > 0) {
          rep(0, n_time_invariant_covar)
        },

        # --- Multi-level random slopes ---
        ms_sd_level_slope_01 = if (enable_ms_01 && n_time_invariant_covar > 0) {
          lapply(seq_len(n_levels), function(lv) rep(0.1, n_time_invariant_covar))
        },
        ms_sd_level_slope_02 = if (enable_ms_02 && n_time_invariant_covar > 0) {
          lapply(seq_len(n_levels), function(lv) rep(0.1, n_time_invariant_covar))
        },
        ms_sd_level_slope_12 = if (enable_ms_12 && n_time_invariant_covar > 0) {
          lapply(seq_len(n_levels), function(lv) rep(0.1, n_time_invariant_covar))
        },
        ms_raw_level_slope_01 = if (enable_ms_01 && n_time_invariant_covar > 0) {
          matrix(0, nrow = n_enabled_groups_ms_slope, ncol = n_time_invariant_covar)
        },
        ms_raw_level_slope_02 = if (enable_ms_02 && n_time_invariant_covar > 0) {
          matrix(0, nrow = n_enabled_groups_ms_slope, ncol = n_time_invariant_covar)
        },
        ms_raw_level_slope_12 = if (enable_ms_12 && n_time_invariant_covar > 0) {
          matrix(0, nrow = n_enabled_groups_ms_slope, ncol = n_time_invariant_covar)
        },
      )
    }) |> purrr::compact() |>
      (\(init_vals) {
        # Save init values to disk if save_dir is provided
        if (!is.null(save_dir)) {
          dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
          # Include run_id in filename to avoid collisions between runs
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
