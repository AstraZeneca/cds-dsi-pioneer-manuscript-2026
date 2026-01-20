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
      # Hierarchical SDs - realistic values based on warmup analysis
      tr_sd_trial <- 0.35
      frac_sd_trial <- 0.35
      init_sd_trial <- 0.35
      tr_sd_patient <- 0.40
      frac_sd_patient <- 0.40
      init_sd_patient <- 0.50

      # Trial-level deviations - small spread
      tr_trial_dev <- rnorm(n_trials, sd = 0.2)
      frac_trial_dev <- rnorm(n_trials, sd = 0.2)
      init_trial_dev <- rnorm(n_trials, sd = 0.2)

      # Patient-level deviations - need heterogeneity for model to fit data
      tr_patient_dev <- rnorm(n_patients, sd = 0.3)
      frac_patient_dev <- rnorm(n_patients, sd = 0.3)
      init_patient_dev <- rnorm(n_patients, sd = 0.3)

      tibble::lst(
        # Population-level parameters - initialize at prior means to avoid catastrophic random inits
        tr_intercept_pop = -2.0,  # Total rate on log scale (prior mean)
        frac_log_decrease_pop = 1.5,  # Fraction log decrease (prior mean on logit scale)
        init_logit_frac_pop = 0.0,  # Initial fraction logit (prior mean)

        # Trial-level SDs and raw values (with heterogeneity)
        tr_sd_trial_intercept = tr_sd_trial,
        tr_raw_trial_intercept = if (enable_trial_intercept_tr) tr_trial_dev / tr_sd_trial,
        frac_sd_trial_intercept = frac_sd_trial,
        frac_raw_trial_intercept = if (enable_trial_intercept_frac) frac_trial_dev / frac_sd_trial,
        init_sd_trial_intercept = init_sd_trial,
        init_raw_trial_intercept = if (enable_trial_intercept_init) init_trial_dev / init_sd_trial,

        # Patient-level SDs and raw values (with heterogeneity - critical for model fit)
        tr_sd_patient_intercept = tr_sd_patient,
        tr_raw_patient_intercept = if (enable_patient_intercept_tr) tr_patient_dev / tr_sd_patient,
        frac_sd_patient_intercept = frac_sd_patient,
        frac_raw_patient_intercept = if (enable_patient_intercept_frac) frac_patient_dev / frac_sd_patient,
        init_sd_patient_intercept = init_sd_patient,
        init_raw_patient_intercept = if (enable_patient_intercept_init) init_patient_dev / init_sd_patient,

        # Covariate effects - all zero
        # Note: use enable_pop_cov_* flags (not enable_trial_cov_*) for population-level coefficients
        tr_coef_qr_pop = if (n_covar > 0 && enable_pop_cov_tr) array(rep(0, n_covar), dim = n_covar),
        frac_coef_qr_pop = if (n_covar > 0 && enable_pop_cov_frac) array(rep(0, n_covar), dim = n_covar),
        init_coef_qr_pop = if (n_covar > 0 && enable_pop_cov_init) array(rep(0, n_covar), dim = n_covar),

        tr_sd_trial_slope = if (n_covar > 0 && enable_trial_cov_tr) array(rep(0.1, n_covar), dim = n_covar),
        tr_raw_trial_slope = if (n_covar > 0 && enable_trial_cov_tr) matrix(0, n_trials, n_covar),
        frac_sd_trial_slope = if (n_covar > 0 && enable_trial_cov_frac) array(rep(0.1, n_covar), dim = n_covar),
        frac_raw_trial_slope = if (n_covar > 0 && enable_trial_cov_frac) matrix(0, n_trials, n_covar),
        init_sd_trial_slope = if (n_covar > 0 && enable_trial_cov_init) array(rep(0.1, n_covar), dim = n_covar),
        init_raw_trial_slope = if (n_covar > 0 && enable_trial_cov_init) matrix(0, n_trials, n_covar),

        # Patient-level covariate SDs and raw effects (with small heterogeneity)
        tr_sd_patient_slope = if (n_covar > 0 && enable_patient_cov_tr) array(rep(0.2, n_covar), dim = n_covar),
        tr_raw_patient_slope = if (n_covar > 0 && enable_patient_cov_tr) matrix(rnorm(n_patients * n_covar, sd = 0.3), n_patients, n_covar),
        frac_sd_patient_slope = if (n_covar > 0 && enable_patient_cov_frac) array(rep(0.2, n_covar), dim = n_covar),
        frac_raw_patient_slope = if (n_covar > 0 && enable_patient_cov_frac) matrix(rnorm(n_patients * n_covar, sd = 0.3), n_patients, n_covar),
        init_sd_patient_slope = if (n_covar > 0 && enable_patient_cov_init) array(rep(0.2, n_covar), dim = n_covar),
        init_raw_patient_slope = if (n_covar > 0 && enable_patient_cov_init) matrix(rnorm(n_patients * n_covar, sd = 0.3), n_patients, n_covar),

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
        measure_sd = 0.15,

        # Other events baseline hazard - realistic values
        log_lambda_gp_pop_intercept = array(rep(-4.5, n_causes), dim = n_causes),
        log_lambda_gp_pop_alpha = array(rep(1.0, n_causes), dim = n_causes),
        log_lambda_gp_pop_rho = array(rep(1.4, n_causes), dim = n_causes),
        # GP eta (time effects) - small heterogeneity for smoother start
        log_lambda_gp_pop_eta = replicate(n_causes, rnorm(max_all_t, sd = 0.1), simplify = FALSE),

        # Other events trial-level baseline hazard
        log_lambda_gp_trial_alpha = if (oe_enable_trial_baseline_hazard) {
          rep(1.0, n_trials)
        },
        log_lambda_gp_trial_rho = if (oe_enable_trial_baseline_hazard) {
          rep(1.4, n_trials)
        },
        log_lambda_gp_trial_intercept_sd = if (oe_enable_trial_baseline_hazard) rep(0.1, n_causes),
        raw_log_lambda_gp_trial_intercept = if (oe_enable_trial_baseline_hazard) {
          array(replicate(n_causes, rep(0, n_trials), simplify = FALSE), dim = c(n_causes, n_trials))
        },
        log_lambda_gp_trial_eta = if (oe_enable_trial_baseline_hazard) {
          replicate(n_causes, matrix(rnorm(n_trials * max_all_t, sd = 0.1), n_trials, max_all_t), simplify = FALSE)
        },

        # Other events covariate effects - all zero
        # Note: tumor coefficients are NOT QR-transformed (unlike oe_covar_coef_qr_pop)
        oe_tumor_coef_pop = if (n_tumor_covar > 0 && oe_enable_pop_tumor_cov) {
          array(replicate(n_causes, rep(0, n_tumor_covar)), dim = c(n_causes, n_tumor_covar))
        },
        oe_sd_trial_tumor_slope = if (n_tumor_covar > 0 && oe_enable_trial_tumor_cov) {
          array(replicate(n_causes, rep(0.1, n_tumor_covar)), dim = c(n_causes, n_tumor_covar))
        },
        oe_raw_trial_tumor_slope = if (n_tumor_covar > 0 && oe_enable_trial_tumor_cov) {
          array(replicate(n_causes, matrix(0, n_trials, n_tumor_covar), simplify = FALSE), dim = c(n_causes, n_trials, n_tumor_covar))
        },

        oe_covar_coef_qr_pop = if (n_covar > 0 && oe_enable_pop_cov) {
          array(replicate(n_causes, rep(0, n_covar)), dim = c(n_causes, n_covar))
        },
        oe_sd_trial_slope = if (n_covar > 0 && oe_enable_trial_cov) {
          array(replicate(n_causes, rep(0.1, n_covar)), dim = c(n_causes, n_covar))
        },
        oe_raw_trial_slope = if (n_covar > 0 && oe_enable_trial_cov) {
          array(replicate(n_causes, matrix(0, n_trials, n_covar), simplify = FALSE), dim = c(n_causes, n_trials, n_covar))
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
