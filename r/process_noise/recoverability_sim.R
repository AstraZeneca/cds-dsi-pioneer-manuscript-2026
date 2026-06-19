# =============================================================================
# Recoverability simulation for the static initial-fraction compartment
# (Task 10 of docs/superpowers/plans/2026-06-12-static-init-fraction.md)
# =============================================================================
#
# GATING VALIDATION. Simulate N patients from a KNOWN population baseline split
#   pi_decrease ~ 0.4, pi_static ~ 0.3, pi_growth ~ 0.3
# using the model's own generative form
#   SLD(t) = B0 * ( pi_decrease * exp(-r_d * t)      # exp(state_decrease)
#                 + pi_static                         # static compartment (rate = 0)
#                 + pi_growth   * exp( r_g * t) )     # exp(state_growth)
# fit stan/tumor/sf-ssm-log-space.stan with enable_static_init = 1L, and check
# whether the population static fraction (init_pi_static_pop) is recovered.
#
# This is exploratory analysis, NOT a unit test: it produces a short table and a
# verdict, not a green/red assertion. SYNTHETIC data only -- no data-prep pipeline.
#
# Run:  Rscript r/process_noise/recoverability_sim.R
#
# -----------------------------------------------------------------------------
# RESULT (recorded 2026-06-12, N = 300, 9 visits to week 48, 2x500/500, 0 div):
#
#                 variable  true   post.mean   90% CI            covers?
#   init_pi_decrease_pop    0.40   0.615       [0.576, 0.655]    NO  (inflated)
#   init_pi_static_pop      0.30   0.011       [0.005, 0.017]    NO  (COLLAPSED ~0)
#   init_pi_growth_pop      0.30   0.374       [0.334, 0.414]    NO  (slightly high)
#
#   VERDICT: static fraction is NOT recoverable at this configuration. The
#   posterior collapses the static mass to ~1%, and the truly-flat patients are
#   absorbed into the DECREASE compartment (pi_decrease inflated 0.40 -> 0.62 ~
#   0.40 + most of the 0.30 static mass). Zero divergences and tight CIs => this
#   is a genuine identifiability result, NOT a sampling artefact.
#
#   ROOT CAUSE: with free per-patient regression/growth rates (tr/frac patient
#   RE), a flat trajectory is observationally near-equivalent to a decrease
#   compartment whose rate ~ 0. The tumor-burden likelihood alone cannot
#   separate "static compartment" from "decrease compartment with rate -> 0",
#   so the extra static degree of freedom is not informed and shrinks to the
#   prior's lower tail. Longer follow-up does not help: a flat patient stays
#   flat, which the decay-rate->0 decrease compartment reproduces exactly.
#
#   GATE OUTCOME: FAIL. The static compartment is not separately identifiable
#   from longitudinal SLD under this parameterization. Viable approaches before
#   any real fit would need EITHER (a) a hard rate-sign / lower-bound constraint
#   that forbids the decrease compartment from going flat, OR (b) external
#   information that pins the static fraction (e.g. tying it to RECIST stable-
#   disease status), OR (c) a tighter, more informative static prior justified
#   by clinical priors rather than the data. As-is, flipping enable_static_init
#   on a real fit would just relabel decrease mass, not recover stable disease.
# -----------------------------------------------------------------------------

suppressWarnings(suppressMessages({
  if (!exists("init_project")) source(here::here(".Rprofile"))
  init_project()
  # Explicitly load the sclc tumor helpers. init_project() only sources them
  # when TAR_PROJECT == "sclc"; for other projects the generic
  # r/prepare_analysis_data.R (tumor-history based) shadows the version this sim
  # needs (visit_data based prepare_tumor_stan_data + elicited-prior helpers).
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
# 1. Known population truth
# ---------------------------------------------------------------------------
N_PATIENTS <- 300L

PI_DECREASE_TRUE <- 0.40
PI_STATIC_TRUE   <- 0.30
PI_GROWTH_TRUE   <- 0.30
stopifnot(abs(PI_DECREASE_TRUE + PI_STATIC_TRUE + PI_GROWTH_TRUE - 1) < 1e-9)

# Realistic per-week log-rates (decrease shrinks, growth expands). These are the
# population tr/frac means; per-patient rates get a small log-normal jitter so
# the trajectories are not identical (mirrors the hierarchical model).
DECREASE_RATE_MEAN <- 0.06   # ~6% shrinkage per week of the decreasing fraction
GROWTH_RATE_MEAN   <- 0.03   # ~3% growth per week of the growing fraction
RATE_PATIENT_SD    <- 0.30   # log-normal spread of per-patient rates

# Measurement noise on log(SLD): matches measure_sd_sld prior mode (~0.125).
MEASURE_SD_TRUE <- 0.12

# Baseline SLD (cm). Log-normal around ~7 cm (= 70 mm), realistic for solid tumours.
BASELINE_LOG_MEAN <- log(7)
BASELINE_LOG_SD   <- 0.5

# ---------------------------------------------------------------------------
# 2. Visit schedule (realistic: baseline + 6-weekly assessments to ~week 48)
# ---------------------------------------------------------------------------
VISIT_WEEKS <- as.integer(seq(0L, 48L, by = 6L))  # weeks 0,6,...,48 -> 9 visits

# ---------------------------------------------------------------------------
# 3. Generate synthetic patients from the generative form
# ---------------------------------------------------------------------------
simulate_patient <- function(i) {
  b0 <- exp(rnorm(1, BASELINE_LOG_MEAN, BASELINE_LOG_SD))   # baseline SLD (cm)
  # Per-patient rates (log-normal jitter around the population means).
  r_d <- DECREASE_RATE_MEAN * exp(rnorm(1, 0, RATE_PATIENT_SD))
  r_g <- GROWTH_RATE_MEAN   * exp(rnorm(1, 0, RATE_PATIENT_SD))

  # Mean SLD at each visit from the three-compartment generative form.
  comp <- PI_DECREASE_TRUE * exp(-r_d * VISIT_WEEKS) +
    PI_STATIC_TRUE +
    PI_GROWTH_TRUE * exp(r_g * VISIT_WEEKS)
  mean_sld_cm <- b0 * comp

  # Measurement noise on the log scale.
  obs_sld_cm <- exp(log(mean_sld_cm) + rnorm(length(VISIT_WEEKS), 0, MEASURE_SD_TRUE))
  obs_sld_mm <- obs_sld_cm * 10  # prepare_tumor_stan_data divides mmsumdiam by 10

  # RECIST category from % change vs baseline (CR if ~0, PR <= -30%, PD >= +20%).
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
    trial = "sim",
    group = "sim",
    visit_data = list(visit_data),
    patient_t_width = patient_max_t + 1L,
    patient_max_t = patient_max_t,
    calendar_day = 1L,
    calendar_week = 1L,
    # All patients administratively right-censored, no events -> ms_mode "none".
    pfs = patient_max_t,
    right_censored = TRUE,
    det_pfs = patient_max_t,
    det_right_censored = TRUE,
    det_interval_censored = 0L,
    interval_censored = 0L,
    death_week = 0L,
    ms_prog_deterministic = 0L,
    # ms_pattern set directly: every patient admin-censored (no events).
    ms_pattern = factor(
      "admin_censored",
      levels = c(
        "admin_censored", "true_dropout", "progressed_alive",
        "progressed_died", "died_on_trial", "died_off_trial"
      )
    )
  )
}

analysis_data <- map(seq_len(N_PATIENTS), simulate_patient) |> list_rbind()
analysis_data$trial <- factor(analysis_data$trial)
analysis_data$group <- factor(analysis_data$group)

cat(sprintf("Simulated %d patients, %d visits each (weeks %s).\n",
            nrow(analysis_data), length(VISIT_WEEKS),
            paste(range(VISIT_WEEKS), collapse = "-")))

# ---------------------------------------------------------------------------
# 4. Assemble stan-data (no covariates, multistate OFF) via the real helpers
# ---------------------------------------------------------------------------
covar_design_matrix <- array(numeric(0), dim = c(nrow(analysis_data), 0))
colnames(covar_design_matrix) <- character(0)

base_stan_data <- prepare_tumor_stan_data(
  analysis_data,
  covar_design_matrix,
  cond_group = list(),
  pfs_quantiles = c(0.25, 0.5, 0.75),
  extend_max_all_t = max(VISIT_WEEKS) + 50L,
  forecast_observation_interval = 6L,
  group_col = "group"
)

# Multistate is fully disabled in this sim, so zero out the covariate-dimension
# fields. With MS off, the model sizes the time-varying / time-invariant coef
# hyperparam arrays to 0 (or to n_time_varying_covar for the always-declared
# 1->2 block); matching n_*_covar = 0 keeps get_tumor_priors() and the model in
# agreement (otherwise time_varying_coef_01_mean is dims-0 in the model but
# length-3 from the priors).
base_stan_data$n_time_varying_covar <- 0L
base_stan_data$n_time_invariant_covar <- 0L

# Init intercepts at "none" so the baseline split is PURELY population-level and
# init_pi_static_pop is directly comparable to PI_STATIC_TRUE (no patient RE on
# the simplex). tr/frac keep patient-level RE so trajectories vary.
none_mode <- level_intercept_mode["none"]
re_mode   <- level_intercept_mode["re"]

stan_settings <- list(
  fit_tumor_data = TRUE,
  fit_multistate_data = FALSE,
  enable_states_full_grid = FALSE,
  sf_rep_T = 20,
  debug = FALSE,

  # Multistate fully disabled.
  enable_ms_01 = FALSE, enable_ms_02 = FALSE, enable_ms_12 = FALSE,
  enable_ms_03 = FALSE, enable_ms_32 = FALSE,
  ms_time_scale_12 = 1L,
  enable_ms_baseline_trend_01 = 0L,
  enable_ms_pop_time_varying_cov = FALSE,
  enable_ms_pop_time_invariant_cov = FALSE,
  enable_ms_level_cov = c(trial = FALSE, patient = FALSE),
  enable_ms_visit_gated_01 = 0L,
  enable_ms_visit_gated_latent_01 = 0L,
  share_dead_gp_shape = 0L,
  enable_ms_02_time_varying_cov = 0L,
  enable_ms_03_time_invariant_cov = 0L,
  enable_ms_03_time_varying_cov = 0L,
  enable_ms_32_time_invariant_cov = 0L,
  enable_ms_12_entry_covar = 0L,
  enable_ms_32_entry_covar = 0L,
  entry_covar_12 = numeric(0),
  entry_covar_32 = numeric(0),

  # TR / frac: patient-level RE (trajectory variation), no covariates.
  enable_level_intercept_tr = c(trial = none_mode, patient = re_mode),
  enable_level_cov_tr = c(trial = FALSE, patient = FALSE),
  enable_pop_cov_tr = FALSE,
  enable_pop_process_noise_tr = FALSE,
  enable_patient_process_noise_tr = FALSE,
  enable_patient_process_noise_sd_tr = FALSE,
  enable_patient_process_noise_phi_tr = FALSE,

  enable_level_intercept_frac = c(trial = none_mode, patient = re_mode),
  enable_level_cov_frac = c(trial = FALSE, patient = FALSE),
  enable_pop_cov_frac = FALSE,

  # Init: NO level intercepts -> population-only simplex; static compartment ON.
  enable_level_intercept_init = c(trial = none_mode, patient = none_mode),
  enable_level_cov_init = c(trial = FALSE, patient = FALSE),
  enable_pop_cov_init = FALSE,
  enable_static_init = 1L,

  # Gompertz growth-rate decay OFF for this sim, but the master gate and the two
  # per-level flag arrays are now unconditionally required in the Stan data block.
  enable_gr_decay = 0L,
  enable_pop_cov_gr_decay = 0L,
  enable_level_intercept_gr_decay = c(trial = none_mode, patient = none_mode),
  enable_level_cov_gr_decay = c(trial = FALSE, patient = FALSE),

  pfs_timepoints = c(24L, 48L),
  n_pfs_timepoints = 2L,

  n_shards = 1L
)

# Decomposed MS baseline-hazard arrays (all-off -> matches disabled MS).
ms_decomp <- decompose_ms_level_baseline_hazard(c(trial = 0L, patient = 0L))

# Priors (no covariates -> empty QR blocks).
elicited_priors <- prepare_elicited_priors(covar_design_matrix)
tumor_priors <- get_tumor_priors(
  c(stan_settings, base_stan_data),
  elicited_priors,
  covar_design_matrix
)

stan_data <- base_stan_data |>
  add_tumor_priors(tumor_priors) |>
  (\(d) c(d, stan_settings))() |>
  (\(d) c(d, ms_decomp))() |>
  (\(d) c(d, derive_ms_fields(analysis_data, "none")))() |>
  (\(d) c(d, list(lfo_eval_trial = 1L)))()

# ---------------------------------------------------------------------------
# 5. Compile and fit (short run: 2 chains, 500 warmup / 500 sampling)
# ---------------------------------------------------------------------------
mod <- cmdstan_model(
  here::here("stan/tumor/sf-ssm-log-space.stan"),
  include_paths = c(here::here("stan"), here::here("stan/tumor")),
  cpp_options = list(stan_threads = TRUE)
)

# Minimal initializer. The full MS-aware fixed initializer sizes its MS GP /
# slope arrays from the slot config and clashes with this MS-disabled model
# (e.g. log_lambda_gp_01_level_alpha would be length 2 but is size 0 here). We
# only need sane starts for the tumor-side population parameters; cmdstanr
# random-inits the remaining (size-0 MS) parameters, which is a no-op for them.
static_logit_true <- qlogis(PI_STATIC_TRUE / (1 - PI_DECREASE_TRUE))
init_fn <- function(chain_id) {
  list(
    tr_loc_pop = log(DECREASE_RATE_MEAN) + rnorm(1, 0, 0.1),
    frac_logit_loc_pop = rnorm(1, 0, 0.3),
    init_logit_loc_pop = qlogis(PI_DECREASE_TRUE) + rnorm(1, 0, 0.2),
    init_logit_static_loc_pop = as.array(static_logit_true + rnorm(1, 0, 0.2)),
    measure_sd_sld = MEASURE_SD_TRUE,
    log_lod = log(0.1)
  )
}

out_dir <- file.path(tempdir(), "recoverability_sim_fit")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

fit <- mod$sample(
  data = stan_data,
  init = init_fn,
  chains = 2,
  parallel_chains = 2,
  threads_per_chain = 1,
  iter_warmup = 500,
  iter_sampling = 500,
  adapt_delta = 0.95,
  max_treedepth = 12,
  refresh = 100,
  seed = 20260612,
  output_dir = out_dir
)

# ---------------------------------------------------------------------------
# 6. Extract population fractions and compare to truth
# ---------------------------------------------------------------------------
pi_vars <- c("init_pi_decrease_pop", "init_pi_static_pop", "init_pi_growth_pop")
draws_mat <- fit$draws(variables = pi_vars, format = "draws_matrix")

truth_vec <- c(PI_DECREASE_TRUE, PI_STATIC_TRUE, PI_GROWTH_TRUE)

comparison <- tibble(
  variable = pi_vars,
  true_value = truth_vec,
  posterior_mean = map_dbl(pi_vars, ~ mean(draws_mat[, .x])),
  q05 = map_dbl(pi_vars, ~ unname(quantile(draws_mat[, .x], 0.05))),
  q95 = map_dbl(pi_vars, ~ unname(quantile(draws_mat[, .x], 0.95)))
) |>
  mutate(
    ci_covers_truth = true_value >= q05 & true_value <= q95,
    collapsed_to_zero = q95 < 0.02
  )

# Divergences
diag <- fit$diagnostic_summary()
n_divergent <- sum(diag$num_divergent)

cat("\n========================= RECOVERABILITY RESULT =========================\n")
print(as.data.frame(comparison), digits = 3)
cat(sprintf("\nDivergent transitions (post-warmup, all chains): %d\n", n_divergent))

static_row <- comparison |> filter(variable == "init_pi_static_pop")
static_recovered <- static_row$ci_covers_truth && !static_row$collapsed_to_zero

cat(sprintf(
  "\nVERDICT: static fraction %s (90%% CI [%.3f, %.3f] %s true %.3f; collapsed_to_zero = %s)\n",
  if (static_recovered) "RECOVERED" else "NOT recoverable",
  static_row$q05, static_row$q95,
  if (static_row$ci_covers_truth) "covers" else "MISSES",
  static_row$true_value,
  static_row$collapsed_to_zero
))
cat("=========================================================================\n")

invisible(list(comparison = comparison, n_divergent = n_divergent,
               static_recovered = static_recovered))
