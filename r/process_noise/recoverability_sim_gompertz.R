# =============================================================================
# Recoverability simulation for Gompertz growth-rate decay (kappa)
# (Task 12 of docs/superpowers/plans/2026-06-12-gompertz-growth-decay.md)
# =============================================================================
#
# DIAGNOSTIC, NOT a unit test. Characterises the kappa-detection threshold under
# realistic follow-up. Simulates N patients from the model's OWN generative form
# with a KNOWN population decay rate kappa_true, warping ONLY the growth arm:
#
#   state_g(t) = init_g + growth_rate * phi(t),   phi(t) = (1 - exp(-kappa*t)) / kappa
#   state_d(t) = init_d - decrease_rate * t            (decrease arm NOT warped)
#   SLD(t)     = B0 * ( pi_decrease * exp(state_d - init_d... )  -- two-component log-space:
#   SLD(t)     = exp(state_d(t)) + exp(state_g(t)) , scaled so SLD(0) = B0.
#
# then fits stan/tumor/sf-ssm-log-space.stan with enable_gr_decay = 1L and checks
# whether the population decay rate (gr_decay_kappa_pop) is recovered and whether
# the posterior CONTRACTS off the weakly-informative prior.
#
# kappa is WEAKLY IDENTIFIED, not aliased (cf. the static compartment in
# recoverability_sim.R which collapsed): d state_g / d kappa ~ -growth_rate*t^2/2
# is nonzero, so signal CAN move it -- but long, post-nadir follow-up is what
# carries that signal. This sim quantifies how much follow-up is needed.
# Cross-ref [[lfo-forecast-regrowth-bias]]: most SCLC/SCLC patients are
# censored while still declining or <=~7wk post-nadir, so the regrowth limb that
# informs kappa is short -> expect weak contraction at empirical follow-up.
#
# SYNTHETIC data only -- no data-prep pipeline.
#
# STATUS (2026-06-12): data generator validated numerically sane (worst-case
# kappa=0 +2sd patient reaches ~12x baseline over the 72-wk window, no
# log_sum_exp overflow). Warmup with the minimal initializer + patient-level RE
# on tr/frac is slow to leave the initial region for this model size (N=300);
# the known-good static sim uses the same minimal init, so a fuller init of the
# patient raw-effect vectors (tr_raw_patient_*, frac_raw_patient_*) and/or a
# longer adapt phase is the likely tuning needed before the table is harvested.
# This is a diagnostic, NOT a code gate -- the warp correctness is proven by the
# committed Stan tests (cross-branch invariance + long-horizon plateau).
#
# Run:  Rscript r/process_noise/recoverability_sim_gompertz.R
# =============================================================================

suppressWarnings(suppressMessages({
  if (!exists("init_project")) source(here::here(".Rprofile"))
  init_project()
  # Explicitly load the sclc tumor helpers (see recoverability_sim.R note):
  # init_project() only sources them when TAR_PROJECT == "sclc".
  source(here::here("r", "sclc", "priors.R"))
  source(here::here("r", "sclc", "prepare_analysis_data.R"))
  source(here::here("r", "sclc", "initializers.R"))
  source(here::here("r", "sclc", "initializers_fixed.R"))
  library(tidyverse)
  library(cmdstanr)
  library(posterior)
}))

set.seed(20260612)
options(cmdstanr_warn_inits = FALSE)

# ---------------------------------------------------------------------------
# 1. Configuration
# ---------------------------------------------------------------------------
N_PATIENTS <- 300L

# Population truth (two-component log-space, baseline-normalised proportions).
# A clear regrowth limb requires a non-trivial growing fraction.
PI_DECREASE_TRUE <- 0.65   # fraction of baseline SLD that responds and shrinks
PI_GROWTH_TRUE   <- 0.35   # fraction that grows back (drives the regrowth limb)
stopifnot(abs(PI_DECREASE_TRUE + PI_GROWTH_TRUE - 1) < 1e-9)

DECREASE_RATE_MEAN <- 0.06   # ~6%/wk shrinkage of the decreasing fraction
# Underlying (un-decayed) growth of the growing fraction. Kept modest so the
# kappa_true = 0 baseline (pure exponential) does not overflow log_sum_exp over
# the follow-up window: at the 72-wk horizon below, exp(0.03*72) ~ 9x, and even
# the +2sd patient (rate ~ 0.055/wk) stays finite. This is the explosive case
# the feature exists to tame, so the data generator must stay numerically sane.
GROWTH_RATE_MEAN   <- 0.03   # ~3%/wk
RATE_PATIENT_SD    <- 0.25   # log-normal spread of per-patient rates

MEASURE_SD_TRUE   <- 0.12    # log(SLD) measurement noise (matches measure_sd_sld mode)
BASELINE_LOG_MEAN <- log(7)  # ~7 cm baseline SLD
BASELINE_LOG_SD   <- 0.5

# kappa sweep: 0 (no decay / pure exponential) -> 0.1 (strong deceleration).
# 0.02 is the prior center (log(0.02)); spans prior center to clearly-bending.
KAPPA_SWEEP <- c(0, 0.01, 0.02, 0.05, 0.1)

# Follow-up: a LONG-but-realistic schedule so the regrowth limb is observed (the
# best case for detecting kappa). Weeks 0,6,...,72 (~1.4 yr) -> 13 visits. Long
# enough to see deceleration, short enough that the kappa=0 baseline stays finite.
VISIT_WEEKS <- as.integer(seq(0L, 72L, by = 6L))

# Gompertz growth-time warp (matches the Stan helper growth_warp()).
phi <- function(t, kappa) if (abs(kappa) < 1e-10) t else (1 - exp(-kappa * t)) / kappa

# ---------------------------------------------------------------------------
# 2. Generative form: warp the GROWTH arm only by phi(t; kappa_true)
# ---------------------------------------------------------------------------
simulate_patient <- function(i, kappa_true) {
  b0  <- exp(rnorm(1, BASELINE_LOG_MEAN, BASELINE_LOG_SD))
  r_d <- DECREASE_RATE_MEAN * exp(rnorm(1, 0, RATE_PATIENT_SD))
  r_g <- GROWTH_RATE_MEAN   * exp(rnorm(1, 0, RATE_PATIENT_SD))

  # Two-component log-space states, normalised so SLD(0) = b0.
  # decrease arm: linear in t (NOT warped). growth arm: warped by phi.
  state_d <- log(PI_DECREASE_TRUE) - r_d * VISIT_WEEKS
  state_g <- log(PI_GROWTH_TRUE)   + r_g * phi(VISIT_WEEKS, kappa_true)
  mean_sld_cm <- b0 * (exp(state_d) + exp(state_g))

  obs_sld_cm <- exp(log(mean_sld_cm) + rnorm(length(VISIT_WEEKS), 0, MEASURE_SD_TRUE))
  obs_sld_mm <- obs_sld_cm * 10  # prepare_tumor_stan_data divides mmsumdiam by 10

  pct <- obs_sld_cm / obs_sld_cm[1] - 1
  resp <- case_when(
    obs_sld_cm < 0.1 * obs_sld_cm[1] ~ "CR",
    pct <= -0.30 ~ "PR",
    pct >= 0.20 ~ "PD",
    TRUE ~ "SD"
  )

  visit_data <- tibble(
    week = VISIT_WEEKS,
    ady = VISIT_WEEKS * 7L + 1L,
    mmsumdiam = obs_sld_mm,
    response = resp,
    det_response = resp
  )

  patient_max_t <- max(VISIT_WEEKS)
  tibble(
    usubjid = sprintf("SIM-%04d", i),
    trial = "sim", group = "sim",
    visit_data = list(visit_data),
    patient_t_width = patient_max_t + 1L,
    patient_max_t = patient_max_t,
    calendar_day = 1L, calendar_week = 1L,
    pfs = patient_max_t, right_censored = TRUE,
    det_pfs = patient_max_t, det_right_censored = TRUE,
    det_interval_censored = 0L, interval_censored = 0L,
    death_week = 0L, ms_prog_deterministic = 0L,
    ms_pattern = factor(
      "admin_censored",
      levels = c("admin_censored", "true_dropout", "progressed_alive",
                 "progressed_died", "died_on_trial", "died_off_trial")
    )
  )
}

# ---------------------------------------------------------------------------
# 3. Assemble stan-data + fit for one kappa_true, return the kappa posterior
# ---------------------------------------------------------------------------
none_mode <- level_intercept_mode["none"]
re_mode   <- level_intercept_mode["re"]

fit_one_kappa <- function(kappa_true, mod) {
  analysis_data <- map(seq_len(N_PATIENTS), simulate_patient, kappa_true = kappa_true) |>
    list_rbind()
  analysis_data$trial <- factor(analysis_data$trial)
  analysis_data$group <- factor(analysis_data$group)

  covar_design_matrix <- array(numeric(0), dim = c(nrow(analysis_data), 0))
  colnames(covar_design_matrix) <- character(0)

  base_stan_data <- prepare_tumor_stan_data(
    analysis_data, covar_design_matrix,
    cond_group = list(), pfs_quantiles = c(0.25, 0.5, 0.75),
    extend_max_all_t = max(VISIT_WEEKS) + 100L,
    forecast_observation_interval = 6L, group_col = "group"
  )
  base_stan_data$n_time_varying_covar <- 0L
  base_stan_data$n_time_invariant_covar <- 0L

  stan_settings <- list(
    fit_tumor_data = TRUE, fit_multistate_data = FALSE,
    enable_states_full_grid = FALSE, sf_rep_T = 20, debug = FALSE,
    enable_ms_01 = FALSE, enable_ms_02 = FALSE, enable_ms_12 = FALSE,
    enable_ms_03 = FALSE, enable_ms_32 = FALSE, ms_time_scale_12 = 1L,
    enable_ms_baseline_trend_01 = 0L,
    enable_ms_pop_time_varying_cov = FALSE, enable_ms_pop_time_invariant_cov = FALSE,
    enable_ms_level_cov = c(trial = FALSE, patient = FALSE),
    enable_ms_visit_gated_01 = 0L, enable_ms_visit_gated_latent_01 = 0L,
    share_dead_gp_shape = 0L,
    enable_ms_02_time_varying_cov = 0L, enable_ms_03_time_invariant_cov = 0L,
    enable_ms_03_time_varying_cov = 0L, enable_ms_32_time_invariant_cov = 0L,
    enable_ms_12_entry_covar = 0L, enable_ms_32_entry_covar = 0L,
    entry_covar_12 = numeric(0), entry_covar_32 = numeric(0),

    enable_level_intercept_tr = c(trial = none_mode, patient = re_mode),
    enable_level_cov_tr = c(trial = FALSE, patient = FALSE),
    enable_pop_cov_tr = FALSE,
    enable_pop_process_noise_tr = FALSE, enable_patient_process_noise_tr = FALSE,
    enable_patient_process_noise_sd_tr = FALSE, enable_patient_process_noise_phi_tr = FALSE,

    enable_level_intercept_frac = c(trial = none_mode, patient = re_mode),
    enable_level_cov_frac = c(trial = FALSE, patient = FALSE),
    enable_pop_cov_frac = FALSE,

    enable_level_intercept_init = c(trial = none_mode, patient = none_mode),
    enable_level_cov_init = c(trial = FALSE, patient = FALSE),
    enable_pop_cov_init = FALSE,
    enable_static_init = 0L,

    # The feature under test: Gompertz growth-rate decay, pop-intercept only.
    enable_gr_decay = 1L,
    enable_pop_cov_gr_decay = 0L,

    pfs_timepoints = c(24L, 48L), n_pfs_timepoints = 2L, n_shards = 1L
  )

  ms_decomp <- decompose_ms_level_baseline_hazard(c(trial = 0L, patient = 0L))
  elicited_priors <- prepare_elicited_priors(covar_design_matrix)
  tumor_priors <- get_tumor_priors(
    c(stan_settings, base_stan_data), elicited_priors, covar_design_matrix
  )

  stan_data <- base_stan_data |>
    add_tumor_priors(tumor_priors) |>
    (\(d) c(d, stan_settings))() |>
    (\(d) c(d, ms_decomp))() |>
    (\(d) c(d, derive_ms_fields(analysis_data, "none")))() |>
    (\(d) c(d, list(lfo_eval_trial = 1L)))()

  # Minimal initializer (tumor-side pop params + the gr_decay intercept).
  init_fn <- function(chain_id) {
    list(
      tr_loc_pop = log(DECREASE_RATE_MEAN) + rnorm(1, 0, 0.1),
      frac_logit_loc_pop = rnorm(1, 0, 0.3),
      init_logit_loc_pop = qlogis(PI_DECREASE_TRUE) + rnorm(1, 0, 0.2),
      gr_decay_log_loc_pop = as.array(log(0.02) + rnorm(1, 0, 0.3)),
      measure_sd_sld = MEASURE_SD_TRUE,
      log_lod = log(0.1)
    )
  }

  out_dir <- file.path(tempdir(), sprintf("gompertz_sim_k%g", kappa_true))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  fit <- mod$sample(
    data = stan_data, init = init_fn,
    chains = 2, parallel_chains = 2, threads_per_chain = 1,
    iter_warmup = 500, iter_sampling = 500,
    adapt_delta = 0.95, max_treedepth = 12,
    refresh = 100, seed = 20260612, output_dir = out_dir
  )

  kd <- fit$draws(variables = "gr_decay_kappa_pop", format = "draws_matrix")
  diag <- fit$diagnostic_summary()
  tibble(
    kappa_true = kappa_true,
    kappa_post_mean = mean(kd[, "gr_decay_kappa_pop"]),
    q05 = unname(quantile(kd[, "gr_decay_kappa_pop"], 0.05)),
    q95 = unname(quantile(kd[, "gr_decay_kappa_pop"], 0.95)),
    n_divergent = sum(diag$num_divergent)
  )
}

# ---------------------------------------------------------------------------
# 4. Compile once, sweep kappa_true
# ---------------------------------------------------------------------------
mod <- cmdstan_model(
  here::here("stan/tumor/sf-ssm-log-space.stan"),
  include_paths = c(here::here("stan"), here::here("stan/tumor")),
  cpp_options = list(stan_threads = TRUE)
)

# Prior 90% interval for kappa (log-normal: log(0.02), sd 0.75) for contraction ref.
prior_q05 <- exp(log(0.02) - qnorm(0.95) * 0.75)
prior_q95 <- exp(log(0.02) + qnorm(0.95) * 0.75)
prior_width <- prior_q95 - prior_q05

results <- map(KAPPA_SWEEP, fit_one_kappa, mod = mod) |> list_rbind() |>
  mutate(
    post_width = q95 - q05,
    ci_covers_truth = kappa_true >= q05 & kappa_true <= q95,
    # Contraction: posterior interval materially narrower than the prior.
    contracted = post_width < 0.6 * prior_width
  )

cat("\n===================== GOMPERTZ kappa RECOVERABILITY =====================\n")
cat(sprintf("Prior 90%% interval for kappa: [%.4f, %.4f] (width %.4f)\n",
            prior_q05, prior_q95, prior_width))
cat(sprintf("Follow-up: %d visits, weeks %d-%d; N = %d/arm.\n\n",
            length(VISIT_WEEKS), min(VISIT_WEEKS), max(VISIT_WEEKS), N_PATIENTS))
print(as.data.frame(results), digits = 3)
cat("\nReading: 'contracted' = posterior interval < 60% of prior width (data moved kappa).\n")
cat("'ci_covers_truth' = 90% posterior CI contains kappa_true.\n")
cat("=========================================================================\n")

invisible(results)
