# nolint start: object_usage_linter

# Dependencies: Requires r/multi_level_hierarchy.R to be loaded
# (sourced by init_project() in .Rprofile)

#' Transform elicited priors from original covariate space to QR space
#'
#' When covariates are correlated, the QR decomposition's R matrix introduces
#' rotation.  Elicited priors encode directional knowledge about individual
#' covariate effects in the **original** (standardized) space.  This function
#' maps those priors into QR space so the Stan model sees correctly rotated
#' prior hyperparameters.
#'
#' @param coef_mean Numeric vector of prior means in original space.
#' @param coef_sd   Numeric vector of prior SDs in original space.
#' @param design_matrix The covariate design matrix (n_patients x n_covar).
#' @return A list with `coef_mean_qr`, `coef_sd_qr`, and `R_stan`.
transform_priors_to_qr_space <- function(coef_mean, coef_sd, design_matrix) {
  n <- nrow(design_matrix)
  R_stan <- qr.R(qr(design_matrix)) / sqrt(n - 1)

  mu_qr <- as.vector(R_stan %*% coef_mean)
  sigma_sq_qr <- diag(R_stan %*% diag(coef_sd^2) %*% t(R_stan))
  sd_qr <- sqrt(sigma_sq_qr)

  list(
    coef_mean_qr = mu_qr,
    coef_sd_qr = sd_qr,
    R_stan = R_stan
  )
}

#' Get shared multistate GP hyperpriors
#'
#' Returns hyperparameters for all multistate transitions (0→1, 0→2, 1→2,
#' 0→3, 3→2) declared unconditionally in multistate/hyperparams.stan.
#' Used by both get_tumor_priors() and get_pioneer_priors().
#'
#' @param n_levels Number of hierarchy levels
#' @param n_time_varying_covar Number of time-varying covariates
#' @param n_time_invariant_covar Number of time-invariant covariates
#' @return Named list of multistate hyperparameter values
get_multistate_priors <- function(n_levels, n_time_varying_covar, n_time_invariant_covar,
                                   enable_ms_visit_gated_01 = 0L,
                                   enable_ms_visit_gated_latent_01 = 0L,
                                   enable_ms_02_time_varying_cov = 1L,
                                   enable_ms_03_time_varying_cov = 0L,
                                   enable_ms_03_time_invariant_cov = 0L,
                                   enable_ms_32_time_invariant_cov = 0L) {
  # 0->1 dimensioning matches stan/modules/multistate/parameters.stan:
  #   single coefficient only in observed visit-gated mode; full
  #   n_time_varying_covar in continuous mode and latent visit-gated mode.
  n_tv_01 <- if (enable_ms_visit_gated_01 && !enable_ms_visit_gated_latent_01) {
    1L
  } else {
    n_time_varying_covar
  }
  n_tv_02 <- if (enable_ms_02_time_varying_cov) n_time_varying_covar else 0L
  n_tv_03 <- if (enable_ms_03_time_varying_cov) n_time_varying_covar else 0L
  n_ti_03 <- if (enable_ms_03_time_invariant_cov) n_time_invariant_covar else 0L
  n_ti_32 <- if (enable_ms_32_time_invariant_cov) n_time_invariant_covar else 0L
  lst(
    # 0→1 transition
    log_lambda_gp_01_pop_intercept_mean = -4.5,
    log_lambda_gp_01_pop_intercept_sd = 1.0,
    log_lambda_gp_01_pop_alpha_alpha = 3.0,
    log_lambda_gp_01_pop_alpha_beta = 1.0,
    # Pop-level rho: invgamma(8, 135) → mode=15w, median≈18w, P(rho>50)=0.7%
    # Raised alpha from 3→8 to tame the heavy right tail of invgamma.
    # With alpha=3 and ~16 rho params, P(any rho>50)≈87% → GP kernel overflow
    # at init. alpha=8 preserves the mode and boundary-avoiding behavior while
    # keeping P(any rho>50)≈3%.
    log_lambda_gp_01_pop_rho_alpha = 8.0,
    log_lambda_gp_01_pop_rho_beta = 135.0,

    # 0→2 transition
    log_lambda_gp_02_pop_intercept_mean = -4.5,
    log_lambda_gp_02_pop_intercept_sd = 1.0,
    log_lambda_gp_02_pop_alpha_alpha = 3.0,
    log_lambda_gp_02_pop_alpha_beta = 1.0,
    log_lambda_gp_02_pop_rho_alpha = 8.0,
    log_lambda_gp_02_pop_rho_beta = 135.0,

    # 1→2 transition (sojourn time GP)
    log_lambda_gp_12_s_pop_intercept_mean = -4.5,
    log_lambda_gp_12_s_pop_intercept_sd = 1.0,
    log_lambda_gp_12_s_pop_alpha_alpha = 3.0,
    log_lambda_gp_12_s_pop_alpha_beta = 1.0,
    log_lambda_gp_12_s_pop_rho_alpha = 8.0,
    log_lambda_gp_12_s_pop_rho_beta = 135.0,

    # 1→2 transition (clock-forward time GP)
    log_lambda_gp_12_t_pop_intercept_mean = -4.5,
    log_lambda_gp_12_t_pop_intercept_sd = 1.0,
    log_lambda_gp_12_t_pop_alpha_alpha = 3.0,
    log_lambda_gp_12_t_pop_alpha_beta = 1.0,
    log_lambda_gp_12_t_pop_rho_alpha = 8.0,
    log_lambda_gp_12_t_pop_rho_beta = 135.0,

    # Shared dead GP shape hyperparameters (used when share_dead_gp_shape=1)
    log_lambda_gp_dead_pop_alpha_alpha = 3.0,
    log_lambda_gp_dead_pop_alpha_beta  = 1.0,
    log_lambda_gp_dead_pop_rho_alpha   = 8.0,
    log_lambda_gp_dead_pop_rho_beta    = 135.0,

    # Level-level rho: invgamma(8, 90) → mode=10w, median≈12w, P(rho>50)=0.05%
    log_lambda_gp_01_level_intercept_sd_sd = rep(0.5, n_levels),
    log_lambda_gp_01_level_alpha_alpha = rep(3.0, n_levels),
    log_lambda_gp_01_level_alpha_beta = rep(1.0, n_levels),
    log_lambda_gp_01_level_rho_alpha = rep(8.0, n_levels),
    log_lambda_gp_01_level_rho_beta = rep(90.0, n_levels),

    log_lambda_gp_02_level_intercept_sd_sd = rep(0.5, n_levels),
    log_lambda_gp_02_level_alpha_alpha = rep(3.0, n_levels),
    log_lambda_gp_02_level_alpha_beta = rep(1.0, n_levels),
    log_lambda_gp_02_level_rho_alpha = rep(8.0, n_levels),
    log_lambda_gp_02_level_rho_beta = rep(90.0, n_levels),

    log_lambda_gp_12_s_level_intercept_sd_sd = rep(0.5, n_levels),
    log_lambda_gp_12_s_level_alpha_alpha = rep(3.0, n_levels),
    log_lambda_gp_12_s_level_alpha_beta = rep(1.0, n_levels),
    log_lambda_gp_12_s_level_rho_alpha = rep(8.0, n_levels),
    log_lambda_gp_12_s_level_rho_beta = rep(90.0, n_levels),

    log_lambda_gp_12_t_level_intercept_sd_sd = rep(0.5, n_levels),
    log_lambda_gp_12_t_level_alpha_alpha = rep(3.0, n_levels),
    log_lambda_gp_12_t_level_alpha_beta = rep(1.0, n_levels),
    log_lambda_gp_12_t_level_rho_alpha = rep(8.0, n_levels),
    log_lambda_gp_12_t_level_rho_beta = rep(90.0, n_levels),

    # 0→3 Dropout GP
    log_lambda_gp_03_pop_intercept_mean = -4.5,
    log_lambda_gp_03_pop_intercept_sd = 1.0,
    log_lambda_gp_03_pop_alpha_alpha = 3.0,
    log_lambda_gp_03_pop_alpha_beta = 1.0,
    log_lambda_gp_03_pop_rho_alpha = 8.0,
    log_lambda_gp_03_pop_rho_beta = 135.0,
    log_lambda_gp_03_level_intercept_sd_sd = rep(2.0, n_levels),
    log_lambda_gp_03_level_alpha_alpha = rep(3.0, n_levels),
    log_lambda_gp_03_level_alpha_beta = rep(1.0, n_levels),
    log_lambda_gp_03_level_rho_alpha = rep(8.0, n_levels),
    log_lambda_gp_03_level_rho_beta = rep(90.0, n_levels),

    # 3→2 Off-trial death GP (sojourn time, semi-Markov)
    log_lambda_gp_32_s_pop_intercept_mean = -4.5,
    log_lambda_gp_32_s_pop_intercept_sd = 1.0,
    log_lambda_gp_32_s_pop_alpha_alpha = 3.0,
    log_lambda_gp_32_s_pop_alpha_beta = 1.0,
    log_lambda_gp_32_s_pop_rho_alpha = 8.0,
    log_lambda_gp_32_s_pop_rho_beta = 135.0,
    log_lambda_gp_32_s_level_intercept_sd_sd = rep(0.5, n_levels),
    log_lambda_gp_32_s_level_alpha_alpha = rep(3.0, n_levels),
    log_lambda_gp_32_s_level_alpha_beta = rep(1.0, n_levels),
    log_lambda_gp_32_s_level_rho_alpha = rep(8.0, n_levels),
    log_lambda_gp_32_s_level_rho_beta = rep(90.0, n_levels),

    # Fixed-effect prior SD for level intercepts (used when flag == 1)
    fe_log_lambda_gp_01_level_intercept_sd = rep(1.0, n_levels),
    fe_log_lambda_gp_02_level_intercept_sd = rep(1.0, n_levels),
    fe_log_lambda_gp_12_s_level_intercept_sd = rep(1.0, n_levels),
    fe_log_lambda_gp_12_t_level_intercept_sd = rep(1.0, n_levels),
    fe_log_lambda_gp_03_level_intercept_sd = rep(1.0, n_levels),
    fe_log_lambda_gp_32_s_level_intercept_sd = rep(1.0, n_levels),

    # Time-varying covariate coefficient hyperparameters
    # as.array() ensures length-1 vectors are passed as 1-element arrays to CmdStan,
    # not as scalars — required when visit-gated mode sets vector size to 1.
    time_varying_coef_01_mean = as.array(rep(0, n_tv_01)),
    time_varying_coef_01_sd = as.array(rep(0.5, n_tv_01)),
    time_varying_coef_02_mean = as.array(rep(0, n_tv_02)),
    time_varying_coef_02_sd = as.array(rep(0.5, n_tv_02)),
    time_varying_coef_12_mean = as.array(rep(0, n_time_varying_covar)),
    time_varying_coef_12_sd = as.array(rep(0.5, n_time_varying_covar)),
    time_varying_coef_03_mean = as.array(rep(0, n_tv_03)),
    time_varying_coef_03_sd = as.array(rep(0.5, n_tv_03)),

    # Time-invariant covariate coefficient hyperparameters (QR space)
    time_invariant_coef_01_mean = rep(0, n_time_invariant_covar),
    time_invariant_coef_01_sd = rep(1, n_time_invariant_covar),
    time_invariant_coef_02_mean = rep(0, n_time_invariant_covar),
    time_invariant_coef_02_sd = rep(1, n_time_invariant_covar),
    time_invariant_coef_12_mean = rep(0, n_time_invariant_covar),
    time_invariant_coef_12_sd = rep(1, n_time_invariant_covar),
    time_invariant_coef_03_mean = rep(0, n_ti_03),
    time_invariant_coef_03_sd = rep(1, n_ti_03),
    time_invariant_coef_32_mean = rep(0, n_ti_32),
    time_invariant_coef_32_sd = rep(1, n_ti_32),

    # Multi-level random slope SD hyperpriors
    sd_level_slope_01_sd = lapply(seq_len(n_levels), function(lv) rep(0.15, n_time_invariant_covar)),
    sd_level_slope_02_sd = lapply(seq_len(n_levels), function(lv) rep(0.15, n_time_invariant_covar)),
    sd_level_slope_12_sd = lapply(seq_len(n_levels), function(lv) rep(0.15, n_time_invariant_covar)),
    sd_level_slope_03_sd = lapply(seq_len(n_levels), function(lv) rep(0.15, n_time_invariant_covar)),
    sd_level_slope_32_sd = lapply(seq_len(n_levels), function(lv) rep(0.15, n_time_invariant_covar)),

    # Burden-at-state-entry covariate coefficient hyperparameters
    # Normal(0, 0.5): matches the scale of standardized log-burden (~1 IQR unit
    # = 1 SD after standardization), so a unit change in standardized log-burden
    # gives a hazard ratio of exp(±0.5) ≈ 1.65, allowing meaningful but not
    # extreme effects. PSA (pioneer) is the only model wired up to use these.
    coef_log_entry_covar_12_mean = 0.0,
    coef_log_entry_covar_12_sd   = 0.5,
    coef_log_entry_covar_32_mean = 0.0,
    coef_log_entry_covar_32_sd   = 0.5,

    # Student-t hierarchy nu hyperparameters for multistate
    ms_nu_baseline_level_prior_alpha = rep(2, n_levels),
    ms_nu_baseline_level_prior_beta  = rep(0.1, n_levels),
    ms_nu_slope_level_prior_alpha    = rep(2, n_levels),
    ms_nu_slope_level_prior_beta     = rep(0.1, n_levels),

    # Correlated intercept block (cross-transition frailty) LKJ shape.
    # eta = 2 gently concentrates toward the identity (rho ~ 0), symmetric about
    # zero — does not bake in the expected negative sign; the data reveal it.
    # Shared by every configured (level, group) block. Harmless when no block is
    # configured (no L_ms_intercept_corr parameter exists in that case).
    ms_intercept_corr_eta = 2
  )
}

get_tumor_priors <- function(stan_data, coef_elicited_priors,
                             covar_design_matrix = NULL) {
  # Directly specified priors (simplified)
  # Choose log-total rate prior similar to historical center; adjust if needed.
  # Updated: shifted mean from -1.0 to -2.0 to reduce prior-posterior conflict
  tr_loc_pop_mean <- -2.0 # centered closer to typical posterior (~-2.6)
  tr_loc_pop_sd <- 0.8 # tightened from 1.2 to reduce total variance with hierarchies
  # Fraction (logit) prior - updated to reduce severe prior-posterior conflict
  # Old prior at -1.0 (~27% fraction) conflicted with posterior at +2.4 (~92%)
  frac_logit_loc_pop_mean <- 1.5 # shifted from -1.0 to reduce 4+ SD conflict
  frac_logit_loc_pop_sd <- 1.0 # slightly wider to allow data to inform
  # Initial proportion logit
  init_logit_loc_pop_mean <- 0.0 # formerly pop_decrease_prop_logis_mean
  init_logit_loc_pop_sd <- 0.8 # tightened from 1.5 (was extremely wide!)

  # Get n_levels from stan_data (default 2 for backward compat)
  n_levels <- stan_data$n_levels %||% 2L
  n_covar <- stan_data$n_covar

  # Multistate covariate dimensions (default 0 if not specified)
  n_time_varying_covar <- stan_data$n_time_varying_covar %||% 0L
  n_time_invariant_covar <- stan_data$n_time_invariant_covar %||% n_covar

  # Transform elicited priors to QR space (frac and init modules)
  # TR and MS use isotropic N(0,1) priors which are rotation-invariant
  if (n_covar > 0 && !is.null(covar_design_matrix)) {
    qr_frac <- transform_priors_to_qr_space(
      coef_elicited_priors$coef_mean,
      coef_elicited_priors$coef_sd,
      covar_design_matrix
    )
    qr_init <- transform_priors_to_qr_space(
      coef_elicited_priors$coef_mean,
      coef_elicited_priors$coef_sd,
      covar_design_matrix
    )
  } else {
    # Fallback: use original-space priors directly (no design matrix available)
    qr_frac <- list(
      coef_mean_qr = coef_elicited_priors$coef_mean,
      coef_sd_qr = coef_elicited_priors$coef_sd
    )
    qr_init <- list(
      coef_mean_qr = coef_elicited_priors$coef_mean,
      coef_sd_qr = coef_elicited_priors$coef_sd
    )
  }

  lst(
    # GP hyperparameters
    pop_tumor_gp_rho_meanlog = 3,
    pop_tumor_gp_rho_sdlog = 0.6,
    log_patient_tumor_gp_rho_sd_sd = 1.75,
    pop_decrease_process_sd_sd = 0.1,
    pop_growth_process_sd_sd = 0.1,
    process_corr_param = 2.0,
    # inv_gamma prior for measure_sd_sld keeps mass away from zero
    # mode = beta/(alpha+1) = 0.75/6 = 0.125 (at typical posterior)
    measure_sd_sld_alpha = 5,
    measure_sd_sld_beta = 0.75,
    # Process parameters
    decrease_process_alpha = 9.7,
    decrease_process_beta = 38.4,
    growth_process_alpha = 9.7,
    growth_process_beta = 38.4,

    # Total rate module hyperparams (multi-level): c(trial, patient)
    tr_loc_pop_mean = tr_loc_pop_mean,
    tr_loc_pop_sd = tr_loc_pop_sd,
    tr_sd_level_intercept_sd = rep(0.35, n_levels),
    tr_fe_sd_level_intercept = rep(0, n_levels),
    # SD sub-hierarchy hyperparameters (issue #110)
    tr_log_sd_level_intercept_pop_mean  = rep(log(1.0), n_levels),
    tr_log_sd_level_intercept_pop_sd    = rep(0.5, n_levels),
    tr_sd_hyperscale_level_intercept_sd = matrix(0.3, n_levels, n_levels),
    tr_fe_sd_hyperscale_level_intercept = matrix(0.3, n_levels, n_levels),
    tr_nu_level_prior_alpha = rep(2, n_levels),
    tr_nu_level_prior_beta = rep(0.1, n_levels),
    tr_coef_qr_pop_mean = as.array(rep(0, n_covar)),
    tr_coef_qr_pop_sd = as.array(rep(1, n_covar)),
    tr_sd_level_slope_sd = list(
      trial = rep(0.15, n_covar),
      patient = rep(0.10, n_covar)
    ),

    # Patient-level process noise (AR(1) time-varying rates per patient)
    # Note: These are deviations in log-rates (decrease/growth), which integrate over time
    # Even small rate deviations accumulate into substantial tumor trajectory effects
    tr_log_sd_pop_process_noise_mean      = -3,    # log(0.05) ≈ -3, median σ ≈ 0.05 (5% deviations)
    tr_log_sd_pop_process_noise_sd        = 0.5,   # allows σ ~[0.02, 0.13] (95% CI)
    tr_logit_phi_pop_process_noise_mean   = 1.4,   # logit(0.8) ≈ 1.39, favors high correlation
    tr_logit_phi_pop_process_noise_sd     = 0.5,   # allows φ ~[0.6, 0.9] (95% CI)
    tr_log_sd_patient_process_noise_sd    = 0.3,   # patient-level variation in log(σ)
    tr_phi_patient_process_noise_sd       = 0.3,   # patient-level variation on logit(φ) scale

    # Population-level process noise (shared AR(1) temporal trend across all patients)
    tr_log_sd_pop_process_noise_pop_mean  = -3,    # log(0.05), conservative
    tr_log_sd_pop_process_noise_pop_sd    = 0.5,   # allows σ ~[0.02, 0.13]
    tr_logit_phi_pop_process_noise_pop_mean = 1.4, # favors high correlation
    tr_logit_phi_pop_process_noise_pop_sd = 0.5,   # allows φ ~[0.6, 0.9]

    # Fraction module hyperparams (multi-level): c(trial, patient)
    frac_logit_loc_pop_mean = frac_logit_loc_pop_mean,
    frac_logit_loc_pop_sd = frac_logit_loc_pop_sd,
    frac_sd_level_intercept_sd = rep(0.25, n_levels),
    frac_fe_sd_level_intercept = rep(0, n_levels),
    frac_nu_level_prior_alpha = rep(2, n_levels),
    frac_nu_level_prior_beta = rep(0.1, n_levels),
    frac_coef_qr_pop_mean = as.array(qr_frac$coef_mean_qr),
    frac_coef_qr_pop_sd = as.array(qr_frac$coef_sd_qr),
    frac_sd_level_slope_sd = list(
      trial = rep(0.05, n_covar),
      patient = rep(0.03, n_covar)
    ),

    # Initial state proportion module hyperparams (multi-level): c(trial, patient)
    init_logit_loc_pop_mean = init_logit_loc_pop_mean,
    init_logit_loc_pop_sd = init_logit_loc_pop_sd,
    init_sd_level_intercept_sd = c(0.6, 0.5),
    init_fe_sd_level_intercept = rep(0, n_levels),
    init_nu_level_prior_alpha = rep(2, n_levels),
    init_nu_level_prior_beta = rep(0.1, n_levels),
    init_coef_qr_pop_mean = as.array(qr_init$coef_mean_qr),
    init_coef_qr_pop_sd = as.array(qr_init$coef_sd_qr),
    init_sd_level_slope_sd = list(
      trial = rep(0.10, n_covar),
      patient = rep(0.08, n_covar)
    ),

    # Growth lag
    growth_lag_mean = 2.7,
    growth_lag_sd = 0.5,
    patient_log_growth_lag_sd_sd = 1,
    log_growth_transition_rate_sd = 1,

    rate_corr_param = 2.0,

    # =========================================================================
    # Multistate Hazard Model Hyperparameters
    # =========================================================================
    # Supports configurable transitions:
    #   - 0→1: Progression / PFS event
    #   - 0→2: Death without progression
    #   - 1→2: Post-progression death (sojourn _s and clock-forward _t GPs)

    # --- Population-level baseline hazard GP hyperparameters ---
    # Relaxed priors to avoid divergences
    # inv_gamma(3, 1.0) - wider prior on alpha
    # inv_gamma(5, 5.0) - wider prior on rho

    # 0→1 transition
    log_lambda_gp_01_pop_intercept_mean = -4.5,
    log_lambda_gp_01_pop_intercept_sd = 1.0,  # Wider intercept prior
    log_lambda_gp_01_pop_alpha_alpha = 3.0,   # Relaxed
    log_lambda_gp_01_pop_alpha_beta = 1.0,
    log_lambda_gp_01_pop_rho_alpha = 5.0,     # Relaxed
    log_lambda_gp_01_pop_rho_beta = 5.0,

    # 0→2 transition
    log_lambda_gp_02_pop_intercept_mean = -4.5,
    log_lambda_gp_02_pop_intercept_sd = 1.0,
    log_lambda_gp_02_pop_alpha_alpha = 3.0,
    log_lambda_gp_02_pop_alpha_beta = 1.0,
    log_lambda_gp_02_pop_rho_alpha = 5.0,
    log_lambda_gp_02_pop_rho_beta = 5.0,

    # 1→2 transition (sojourn time GP)
    log_lambda_gp_12_s_pop_intercept_mean = -4.5,
    log_lambda_gp_12_s_pop_intercept_sd = 1.0,
    log_lambda_gp_12_s_pop_alpha_alpha = 3.0,
    log_lambda_gp_12_s_pop_alpha_beta = 1.0,
    log_lambda_gp_12_s_pop_rho_alpha = 5.0,
    log_lambda_gp_12_s_pop_rho_beta = 5.0,

    # 1→2 transition (clock-forward time GP)
    log_lambda_gp_12_t_pop_intercept_mean = -4.5,
    log_lambda_gp_12_t_pop_intercept_sd = 1.0,
    log_lambda_gp_12_t_pop_alpha_alpha = 3.0,
    log_lambda_gp_12_t_pop_alpha_beta = 1.0,
    log_lambda_gp_12_t_pop_rho_alpha = 5.0,
    log_lambda_gp_12_t_pop_rho_beta = 5.0,

    # --- Level-level baseline hazard GP hyperparameters (N-level hierarchy) ---
    # Relaxed priors to avoid divergences

    # 0→1 transition (per level)
    log_lambda_gp_01_level_intercept_sd_sd = rep(0.5, n_levels),  # Wider
    log_lambda_gp_01_level_alpha_alpha = rep(3.0, n_levels),      # Relaxed
    log_lambda_gp_01_level_alpha_beta = rep(1.0, n_levels),
    log_lambda_gp_01_level_rho_alpha = rep(4.0, n_levels),        # Relaxed
    log_lambda_gp_01_level_rho_beta = rep(4.0, n_levels),

    # 0→2 transition (per level)
    log_lambda_gp_02_level_intercept_sd_sd = rep(0.5, n_levels),
    log_lambda_gp_02_level_alpha_alpha = rep(3.0, n_levels),
    log_lambda_gp_02_level_alpha_beta = rep(1.0, n_levels),
    log_lambda_gp_02_level_rho_alpha = rep(4.0, n_levels),
    log_lambda_gp_02_level_rho_beta = rep(4.0, n_levels),

    # 1→2 sojourn time GP (per level)
    log_lambda_gp_12_s_level_intercept_sd_sd = rep(0.5, n_levels),
    log_lambda_gp_12_s_level_alpha_alpha = rep(3.0, n_levels),
    log_lambda_gp_12_s_level_alpha_beta = rep(1.0, n_levels),
    log_lambda_gp_12_s_level_rho_alpha = rep(4.0, n_levels),
    log_lambda_gp_12_s_level_rho_beta = rep(4.0, n_levels),

    # 1→2 clock-forward time GP (per level)
    log_lambda_gp_12_t_level_intercept_sd_sd = rep(0.5, n_levels),
    log_lambda_gp_12_t_level_alpha_alpha = rep(3.0, n_levels),
    log_lambda_gp_12_t_level_alpha_beta = rep(1.0, n_levels),
    log_lambda_gp_12_t_level_rho_alpha = rep(4.0, n_levels),
    log_lambda_gp_12_t_level_rho_beta = rep(4.0, n_levels),

    # --- 0→3 Dropout GP (N-level hierarchy, clock-forward time) ---
    log_lambda_gp_03_pop_intercept_mean = -4.5,
    log_lambda_gp_03_pop_intercept_sd = 1.0,
    log_lambda_gp_03_pop_alpha_alpha = 3.0,
    log_lambda_gp_03_pop_alpha_beta = 1.0,
    log_lambda_gp_03_pop_rho_alpha = 5.0,
    log_lambda_gp_03_pop_rho_beta = 5.0,
    log_lambda_gp_03_level_intercept_sd_sd = rep(2.0, n_levels),
    log_lambda_gp_03_level_alpha_alpha = rep(3.0, n_levels),
    log_lambda_gp_03_level_alpha_beta = rep(1.0, n_levels),
    log_lambda_gp_03_level_rho_alpha = rep(5.0, n_levels),
    log_lambda_gp_03_level_rho_beta = rep(5.0, n_levels),

    # --- 3→2 Baseline Hazard GP: Population-level (sojourn time, semi-Markov) ---
    log_lambda_gp_32_s_pop_intercept_mean = -4.5,
    log_lambda_gp_32_s_pop_intercept_sd = 1.0,
    log_lambda_gp_32_s_pop_alpha_alpha = 3.0,
    log_lambda_gp_32_s_pop_alpha_beta = 1.0,
    log_lambda_gp_32_s_pop_rho_alpha = 5.0,
    log_lambda_gp_32_s_pop_rho_beta = 5.0,

    # --- 3→2 Baseline Hazard GP: Level-level ---
    log_lambda_gp_32_s_level_intercept_sd_sd = rep(0.5, n_levels),
    log_lambda_gp_32_s_level_alpha_alpha = rep(3.0, n_levels),
    log_lambda_gp_32_s_level_alpha_beta = rep(1.0, n_levels),
    log_lambda_gp_32_s_level_rho_alpha = rep(4.0, n_levels),
    log_lambda_gp_32_s_level_rho_beta = rep(4.0, n_levels),

    # --- Time-varying covariate coefficient hyperparameters ---
    time_varying_coef_01_mean = as.array(rep(0, n_time_varying_covar)),
    time_varying_coef_01_sd = as.array(rep(0.5, n_time_varying_covar)),
    time_varying_coef_02_mean = as.array(rep(0, n_time_varying_covar)),
    time_varying_coef_02_sd = as.array(rep(0.5, n_time_varying_covar)),
    time_varying_coef_12_mean = as.array(rep(0, n_time_varying_covar)),
    time_varying_coef_12_sd = as.array(rep(0.5, n_time_varying_covar)),

    # --- Time-invariant covariate coefficient hyperparameters (QR space) ---
    time_invariant_coef_01_mean = rep(0, n_time_invariant_covar),
    time_invariant_coef_01_sd = rep(1, n_time_invariant_covar),
    time_invariant_coef_02_mean = rep(0, n_time_invariant_covar),
    time_invariant_coef_02_sd = rep(1, n_time_invariant_covar),
    time_invariant_coef_12_mean = rep(0, n_time_invariant_covar),
    time_invariant_coef_12_sd = rep(1, n_time_invariant_covar),

    # --- Multi-level random slope SD hyperpriors (per level, per transition) ---
    sd_level_slope_01_sd = lapply(seq_len(n_levels), function(lv) {
      rep(0.15, n_time_invariant_covar)
    }),
    sd_level_slope_02_sd = lapply(seq_len(n_levels), function(lv) {
      rep(0.15, n_time_invariant_covar)
    }),
    sd_level_slope_12_sd = lapply(seq_len(n_levels), function(lv) {
      rep(0.15, n_time_invariant_covar)
    }),

    log_lod_sd = 0.2
  ) |>
    list_assign(!!!get_multistate_priors(n_levels, n_time_varying_covar, n_time_invariant_covar,
                                        enable_ms_visit_gated_01        = stan_data$enable_ms_visit_gated_01 %||% 0L,
                                        enable_ms_visit_gated_latent_01 = stan_data$enable_ms_visit_gated_latent_01 %||% 0L,
                                        enable_ms_02_time_varying_cov   = stan_data$enable_ms_02_time_varying_cov %||% 0L,
                                        enable_ms_03_time_varying_cov   = stan_data$enable_ms_03_time_varying_cov %||% 0L,
                                        enable_ms_03_time_invariant_cov = stan_data$enable_ms_03_time_invariant_cov %||% 0L,
                                        enable_ms_32_time_invariant_cov = stan_data$enable_ms_32_time_invariant_cov %||% 0L))
}

get_pfs_priors <- function() {
  lst(
    log_lambda_gp_intercept_mean = -4.5,
    log_lambda_gp_intercept_sd = 0.25,
    log_lambda_gp_alpha_sd = 0.5,
    log_lambda_gp_rho_alpha = 7.3,
    log_lambda_gp_rho_beta = 7.5,
    log_lambda_gp_trial_alpha_sd = 0.25,
    log_lambda_gp_trial_intercept_sd_sd = 0.25,

    tumor_stim_pop_intercept_sd = 0.35,
    tumor_stim_pop_coef_sd = c(0.15, 0.15, 0.15, 0.05, 0.05),
    tumor_stim_trial_coef_sd_sd = c(0.3, 0.125, 0.125, 0.125, 0.125, 0.125),
    tumor_stim_location_coef_sd_sd = tumor_stim_trial_coef_sd_sd,

    orr_pop_coef_sd = 0.125
  )
}

get_pfs_conf_resp_priors <- function(stan_data) {
  get_pfs_priors() |>
    list_assign(
      log_lambda_gp_intercept_mean = -2.5,
      log_lambda_gp_intercept_sd = 0.5,
      log_lambda_gp_alpha_sd = 1,
      log_lambda_gp_rho_alpha = 3,
      log_lambda_gp_rho_beta = 10,

      covar_effect_mean = rep(0, stan_data$n_covar),
      covar_effect_sd = rep(0.15, stan_data$n_covar),
      tumor_stim_pop_coef_sd = c(0.2, 0.2, 0.2, 0.15, 0.15),
      conf_resp_effect_mean = 0,
      conf_resp_effect_sd = 0.3,

      covar_trial_sd_sd = 0.1,
      covar_trial_corr_eta = 2
    )
}

get_confirmed_resp_priors <- function() {
  lst(
    crcr_covar_effect_mean = 0,
    crcr_covar_effect_sd = 0.2,
    crcr_tumor_stim_pop_coef_sd = c(0.5, 0.5, 0.2, 0.15, 0.15),

    log_crcr_lambda_gp_intercept_mean = rep(-2.5, 2),
    log_crcr_lambda_gp_intercept_sd = 0.75,

    log_crcr_lambda_gp_alpha_sd = 0.6,
    log_crcr_lambda_gp_rho_alpha = 3,
    log_crcr_lambda_gp_rho_beta = 10,
    log_crcr_lambda_gp_trial_alpha_sd = 0.25,
    log_crcr_lambda_gp_trial_intercept_sd_sd = 0.25,

    crcr_covar_trial_sd_sd = 0.1,
    crcr_covar_trial_corr_eta = 2
  )
}

# nolint end: object_usage_linter
