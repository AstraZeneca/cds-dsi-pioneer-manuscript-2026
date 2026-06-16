# Phase 1: 3D bi-exponential TUMOR (SLD) Laplace experiment — GO/NO-GO gate
# ============================================================================
# Validates the built-in laplace_marginal_tol (Stan >= 2.39) on the REAL tumor
# observation model (log-Gaussian SLD + below-LOD left-censoring) before we
# wire it into stan/tumor/sf-ssls-lfo.stan to marginalize the patient-level
# latents of a backgrounded historical trial.
#
# Compares Full HMC (mode 0) vs built-in Laplace (mode 1) on synthetic SLD data
# that INCLUDES censored (below-LOD) visits, since the normal_lcdf tail is the
# part most likely to stress the Gaussian (log-concavity) approximation.
#
# PASS: all population parameters agree within 5 MCSE.
# ============================================================================

library(cmdstanr)
library(posterior)
library(dplyr)
library(tibble)
library(purrr)

# Released CmdStan 2.39 ships laplace_marginal_tol (the old experiments used the
# unreleased cmdstan-develop — no longer needed).
set_cmdstan_path("~/.cmdstan/cmdstan-2.39.0")

# --- Simulate data ---------------------------------------------------------
set.seed(42)
n_patients <- 12
n_visits_per_patient <- 8  # weeks 1, 5, 9, ..., (1-indexed; dt = idx - 1)

true_tr_loc_pop     <- -3.0
true_tr_sd          <- 0.3
true_frac_logit_pop <- 0.5
true_frac_sd        <- 0.3
true_init_logit_pop <- 0.0
true_init_sd        <- 0.3
true_measure_sd     <- 0.15

# Normalized LOD: observations whose simulated value falls below this are
# recorded as censored (encoded as <= 0 in normalized_obs, matching Stan).
lod <- 0.08
log_lod <- log(lod)

z_tr_true   <- rnorm(n_patients)
z_frac_true <- rnorm(n_patients)
z_init_true <- rnorm(n_patients)

visit_weeks <- seq(1, by = 4, length.out = n_visits_per_patient)

all_obs <- c()
all_visit_idx <- c()
patient_visit_pos <- integer(n_patients + 1)
patient_visit_pos[1] <- 1L
n_censored <- 0L

for (i in seq_len(n_patients)) {
  tr_loc     <- true_tr_loc_pop     + true_tr_sd   * z_tr_true[i]
  frac_logit <- true_frac_logit_pop + true_frac_sd * z_frac_true[i]
  init_logit <- true_init_logit_pop + true_init_sd * z_init_true[i]

  log_dec_frac <- plogis(frac_logit, log.p = TRUE)
  log_gro_frac <- plogis(frac_logit, lower.tail = FALSE, log.p = TRUE)
  dec_rate <- exp(tr_loc + log_dec_frac)
  gro_rate <- exp(tr_loc + log_gro_frac)
  init_log_dec <- plogis(init_logit, log.p = TRUE)
  init_log_gro <- plogis(init_logit, lower.tail = FALSE, log.p = TRUE)

  for (v in seq_along(visit_weeks)) {
    dt <- visit_weeks[v] - 1
    state_dec <- init_log_dec - dec_rate * dt
    state_gro <- init_log_gro + gro_rate * dt
    log_pred <- matrixStats::logSumExp(c(state_dec, state_gro))
    log_obs <- rnorm(1, mean = log_pred, sd = true_measure_sd)
    obs <- exp(log_obs)
    if (obs < lod) {
      obs <- 0  # encode as censored
      n_censored <- n_censored + 1L
    }
    all_obs <- c(all_obs, obs)
    all_visit_idx <- c(all_visit_idx, visit_weeks[v])
  }

  patient_visit_pos[i + 1] <- patient_visit_pos[i] + n_visits_per_patient
}

stan_data <- list(
  n_patients = n_patients,
  n_total_visits = length(all_obs),
  normalized_obs = all_obs,
  patient_visit_pos = patient_visit_pos,
  visit_time_idx = all_visit_idx,
  measure_sd = true_measure_sd,
  log_lod = log_lod
)

cat(sprintf("Data: %d patients, %d visits, %d censored (%.0f%%)\n",
            n_patients, length(all_obs), n_censored,
            100 * n_censored / length(all_obs)))

# --- Compile ---------------------------------------------------------------
model <- cmdstan_model(
  "stan/experiments/laplace_tumor_test.stan",
  include_paths = c("stan", "stan/tumor"),
  quiet = FALSE
)

pop_params <- c("tr_loc_pop", "tr_sd", "frac_logit_pop", "frac_sd",
                "init_logit_pop", "init_sd")

# Tight inits near plausible values so the Laplace marginal (inner Newton solve)
# is never evaluated at extreme hyperparameters where it returns -inf. SDs are
# floored away from 0 to keep the NCP scaling well-conditioned.
init_fn <- function() {
  list(
    tr_loc_pop = rnorm(1, -3, 0.1),
    tr_sd = 0.3 + abs(rnorm(1, 0, 0.05)),
    frac_logit_pop = rnorm(1, 0.5, 0.1),
    frac_sd = 0.3 + abs(rnorm(1, 0, 0.05)),
    init_logit_pop = rnorm(1, 0, 0.1),
    init_sd = 0.3 + abs(rnorm(1, 0, 0.05))
  )
}

fit_modes <- list()
mode_names <- c("full_hmc", "laplace")
timings <- list()

for (mode in 0:1) {
  name <- mode_names[mode + 1]
  cat(sprintf("\n=== Fitting mode %d (%s) ===\n", mode, name))

  init_fn_mode <- if (mode == 0) {
    function() c(init_fn(), list(
      z_tr = rnorm(n_patients, 0, 0.5),
      z_frac = rnorm(n_patients, 0, 0.5),
      z_init = rnorm(n_patients, 0, 0.5)
    ))
  } else {
    init_fn
  }

  # Laplace mode: tighter inits (avoid -inf marginal-density starts that killed
  # a chain in the first run) and higher adapt_delta (the marginal density is
  # stiffer because each eval runs an inner Newton solve).
  adapt_delta <- if (mode == 1) 0.99 else 0.9

  t0 <- proc.time()
  fit <- model$sample(
    data = c(stan_data, list(laplace_mode = mode)),
    init = init_fn_mode,
    seed = 123,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 1000,
    iter_sampling = 2000,
    refresh = 500,
    adapt_delta = adapt_delta
  )
  timings[[name]] <- (proc.time() - t0)["elapsed"]
  fit_modes[[name]] <- fit
}

# --- Compare posteriors ----------------------------------------------------
cat("\n\n=== Posterior comparison (population parameters) ===\n\n")
summaries <- map(mode_names, \(name) {
  fit_modes[[name]]$summary(variables = pop_params) |>
    mutate(mode = name, .before = 1)
}) |>
  list_rbind()
print(summaries |> select(mode, variable, mean, sd, q5, q95, rhat, ess_bulk))

cat("\n\n=== HMC vs Laplace comparison ===\n")
hmc_summary <- fit_modes$full_hmc$summary(variables = pop_params)
lap_summary <- fit_modes$laplace$summary(variables = pop_params)

comparison <- tibble(
  parameter = hmc_summary$variable,
  hmc_mean = hmc_summary$mean,
  laplace_mean = lap_summary$mean,
  diff = abs(hmc_summary$mean - lap_summary$mean),
  hmc_sd = hmc_summary$sd,
  laplace_sd = lap_summary$sd,
  mcse = hmc_summary$sd / sqrt(pmin(hmc_summary$ess_bulk, lap_summary$ess_bulk)),
  diff_over_mcse = abs(hmc_summary$mean - lap_summary$mean) /
    (hmc_summary$sd / sqrt(pmin(hmc_summary$ess_bulk, lap_summary$ess_bulk)))
)
print(comparison)

cat(sprintf("\n\nTimings: HMC = %.1fs, Laplace = %.1fs\n",
            timings$full_hmc, timings$laplace))

cat("\n=== True values ===\n")
cat(sprintf("tr_loc_pop:     %.2f\n", true_tr_loc_pop))
cat(sprintf("tr_sd:          %.2f\n", true_tr_sd))
cat(sprintf("frac_logit_pop: %.2f\n", true_frac_logit_pop))
cat(sprintf("frac_sd:        %.2f\n", true_frac_sd))
cat(sprintf("init_logit_pop: %.2f\n", true_init_logit_pop))
cat(sprintf("init_sd:        %.2f\n", true_init_sd))

# --- Convergence gate ------------------------------------------------------
# The diff/MCSE check is only meaningful if BOTH fits actually converged. A
# divergence-riddled, low-ESS Laplace fit inflates MCSE and can fake a "PASS",
# so we gate on sampler health FIRST.
diag_summary <- function(fit, name) {
  d <- fit$diagnostic_summary(quiet = TRUE)
  s <- fit$summary(variables = pop_params)
  n_div <- sum(d$num_divergent)
  n_draws <- fit$metadata()$iter_sampling * length(d$num_divergent)
  tibble(
    mode = name,
    n_divergent = n_div,
    pct_divergent = 100 * n_div / n_draws,
    max_rhat = max(s$rhat, na.rm = TRUE),
    min_ess_bulk = min(s$ess_bulk, na.rm = TRUE)
  )
}
diag_tbl <- bind_rows(
  diag_summary(fit_modes$full_hmc, "full_hmc"),
  diag_summary(fit_modes$laplace, "laplace")
)
cat("\n=== Convergence diagnostics ===\n")
print(diag_tbl)

lap_diag <- diag_tbl |> filter(mode == "laplace")
converged <- lap_diag$pct_divergent < 1 &&
  lap_diag$max_rhat < 1.01 &&
  lap_diag$min_ess_bulk > 400

max_ratio <- max(comparison$diff_over_mcse)
cat(sprintf("\nMax diff/MCSE ratio: %.2f\n", max_ratio))

if (!converged) {
  cat("NO-GO: Laplace fit did not converge (divergences / high R-hat / low ESS).\n")
  cat("       diff/MCSE is NOT trustworthy until this is fixed.\n")
} else if (max_ratio < 5) {
  cat("PASS: Laplace converged AND all population parameters agree within 5 MCSE.\n")
} else {
  cat("FAIL: Laplace converged but parameters differ by more than 5 MCSE",
      "(genuine approximation bias).\n")
}

cat("\nPhase 1 complete.\n")
