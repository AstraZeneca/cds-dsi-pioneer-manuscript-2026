library(testthat)
library(cmdstanr)
library(here)
library(posterior)
library(stringr)

# Backward-compat anchor for the gr_decay (Gompertz growth-rate decay) hierarchy.
#
# The module was extended from a pooled population-scalar kappa to a frac-style
# level RE hierarchy (population -> trial-arm -> patient). This test proves the
# hierarchy is a STRICT GENERALIZATION of the old pooled scalar:
#
#   Config A (pop-only): every level intercept mode = NONE
#     -> kappa_i == exp(gr_decay_log_loc_pop) for every patient.
#   Config B (trial-arm + patient RE, raw draws = 0): RE levels enabled but the
#     raw NCP draws are fixed to 0 -> sd * 0 = 0 contribution on the log scale
#     -> kappa_i must EQUAL Config A's kappa.
#
# Strategy 1: the test Stan model #includes the REAL module fragments
# (modules/gr_decay/{transformed_data,parameters,transformed_parameters}.stan)
# plus minimal hierarchy scaffolding, so this exercises the production code path.
#
# kappa is made deterministic via fixed_param sampling with init = 0 for every
# gr_decay parameter and the pop intercept pinned to log_loc_pop_fixed.

test_that("gr_decay kappa hierarchy reduces to pop scalar at zero RE effect", {
  stan_file <- here::here("tests/testthat/stan/test_gr_decay_hierarchy.stan")

  # Stale-exe guard: drop any previously compiled binary so include changes
  # are always picked up (mirrors sibling gr_decay tests).
  exe <- tools::file_path_sans_ext(stan_file)
  if (file.exists(exe)) file.remove(exe)

  mod <- cmdstan_model(
    stan_file,
    include_paths = c(here::here("stan"), here::here("stan/tumor")),
    quiet = FALSE,
    force_recompile = TRUE
  )

  # ----- synthetic hierarchy: 6 patients across 2 trial-arms, 2 levels -----
  # NO real subject IDs — pure synthetic integer indices.
  n_patients <- 6L
  n_levels <- 2L                      # level 1 = trial-arm, level 2 = patient
  arm_of_patient <- c(1L, 1L, 1L, 2L, 2L, 2L)   # patient -> arm map (synthetic)
  n_arms <- 2L
  n_groups_per_level <- c(n_arms, n_patients)
  # patient_level_groups[, 1] = arm, [, 2] = identity (patient level)
  patient_level_groups <- cbind(arm_of_patient, seq_len(n_patients))
  # patient_trial must match patient_level_groups[, 1] (validated in Stan).
  patient_trial <- arm_of_patient

  log_loc_pop_fixed <- -3.9           # kappa = exp(-3.9) ~ 0.02024
  kappa_expected <- exp(log_loc_pop_fixed)

  base_data <- list(
    n_trials = n_arms,
    n_patients = n_patients,
    patient_trial = patient_trial,
    n_levels = n_levels,
    n_groups_per_level = n_groups_per_level,
    patient_level_groups = patient_level_groups,
    n_covar = 0L,
    covar_design_matrix = matrix(0, n_patients, 0),
    n_forecast_patients = n_patients,
    forecast_patient_idx = seq_len(n_patients),
    enable_gr_decay = 1L,
    enable_pop_cov_gr_decay = 0L,
    enable_level_cov_gr_decay = rep(0L, n_levels),
    # hyperparams (priors-only; FE SDs unused since no FE levels here)
    gr_decay_log_loc_pop_mean = log_loc_pop_fixed,
    gr_decay_log_loc_pop_sd = 1.0,
    gr_decay_coef_qr_pop_mean = numeric(0),
    gr_decay_coef_qr_pop_sd = numeric(0),
    gr_decay_sd_level_intercept_sd = rep(1.0, n_levels),
    gr_decay_fe_sd_level_intercept = rep(1.0, n_levels),
    gr_decay_sd_level_slope_sd = replicate(n_levels, numeric(0), simplify = FALSE)
  )

  run_kappa <- function(level_modes) {
    data <- base_data
    data$enable_level_intercept_gr_decay <- level_modes

    # Number of RE/RE_CP intercept levels -> sizes gr_decay_sd_level_intercept_raw.
    n_re <- sum(level_modes %in% c(2L, 4L))
    # Number of raw (FE/RE/RE_GP) intercept groups -> gr_decay_raw_level_intercept.
    # Modes 1/2/3 contribute their level's forecast group count; 0/4 contribute 0.
    raw_levels <- which(level_modes %in% c(1L, 2L, 3L))
    n_forecast_groups <- n_groups_per_level
    n_forecast_groups[n_levels] <- data$n_forecast_patients
    n_raw_intercept <- sum(n_forecast_groups[raw_levels])
    # CP groups (mode 4) -> gr_decay_cp_level_intercept.
    cp_levels <- which(level_modes == 4L)
    n_cp_intercept <- sum(n_forecast_groups[cp_levels])

    # Pin the pop intercept and force every RE raw NCP draw to 0, so each level's
    # contribution is sd * raw = sd * 0 = 0 on the log scale regardless of sd.
    # The SD free params carry a <lower=0> constraint, so they must be initialised
    # strictly positive (a 0 init lands on the boundary -> unconstrained log(0) =
    # -inf and the chain rejects). 1.0 is arbitrary: it never enters kappa because
    # the raw draw it multiplies is 0. The (n_covar = 0) slope params are 0-sized;
    # cmdstanr auto-inits them — omit to avoid empty matrix-dim round-trip mismatches.
    init_list <- list(
      gr_decay_log_loc_pop = array(log_loc_pop_fixed, dim = 1),
      gr_decay_sd_level_intercept_raw = rep(1, n_re),
      gr_decay_raw_level_intercept = rep(0, n_raw_intercept),
      gr_decay_cp_level_intercept = rep(0, n_cp_intercept)
    )

    fit <- mod$sample(
      data = data,
      init = list(init_list),
      chains = 1,
      iter_sampling = 1,
      iter_warmup = 0,
      fixed_param = TRUE,
      show_messages = FALSE,
      refresh = 0
    )
    draws <- posterior::as_draws_df(fit$draws("kappa_out"))
    vapply(
      seq_len(n_patients),
      function(i) as.numeric(draws[[sprintf("kappa_out[%d]", i)]][1]),
      numeric(1)
    )
  }

  # Config A: pop-only (all level modes NONE).
  kappa_A <- run_kappa(c(0L, 0L))
  expect_equal(kappa_A, rep(kappa_expected, n_patients), tolerance = 1e-6)

  # Config B: trial-arm + patient RE, but raw draws = 0.
  kappa_B <- run_kappa(c(2L, 2L))
  expect_equal(kappa_B, rep(kappa_expected, n_patients), tolerance = 1e-6)

  # The key reduction claim: RE-with-zero-effect == pop-only, exactly.
  expect_equal(kappa_B, kappa_A, tolerance = 1e-6)
})
