# Phase 3: 3D bi-exponential PSA Laplace experiment
# ============================================================================
# Validates manual 3D Laplace approximation with bi-exponential dynamics.
# Compares Full HMC (mode 0) vs Laplace (mode 1) on synthetic PSA data.
# ============================================================================

library(cmdstanr)
library(posterior)
library(dplyr)
library(tibble)
library(purrr)

# Use CmdStan develop (has laplace_marginal_tol support)
set_cmdstan_path("~/.cmdstan/cmdstan-develop")

# --- Simulate data ---
set.seed(42)
n_patients <- 10
n_visits_per_patient <- 6  # visits at weeks 0, 4, 8, ..., 20

# True population parameters
true_tr_loc_pop <- -3.0       # log total rate
true_tr_sd <- 0.3
true_frac_logit_pop <- 0.5    # logit(decrease fraction)
true_frac_sd <- 0.3
true_init_logit_pop <- 0.0    # logit(initial decrease proportion)
true_init_sd <- 0.3
true_measure_sd <- 0.15

# Simulate patient-level NCP params
z_tr_true <- rnorm(n_patients)
z_frac_true <- rnorm(n_patients)
z_init_true <- rnorm(n_patients)

# Visit times: weeks 1, 5, 9, ..., (1-indexed to match model convention)
visit_weeks <- seq(1, by = 4, length.out = n_visits_per_patient)

# Generate observations
all_obs <- c()
all_visit_idx <- c()
patient_visit_pos <- integer(n_patients + 1)
patient_visit_pos[1] <- 1L

for (i in seq_len(n_patients)) {
  tr_loc <- true_tr_loc_pop + true_tr_sd * z_tr_true[i]
  frac_logit <- true_frac_logit_pop + true_frac_sd * z_frac_true[i]
  init_logit <- true_init_logit_pop + true_init_sd * z_init_true[i]

  log_dec_frac <- plogis(frac_logit, log.p = TRUE)
  log_gro_frac <- plogis(frac_logit, lower.tail = FALSE, log.p = TRUE)
  dec_rate <- exp(tr_loc + log_dec_frac)
  gro_rate <- exp(tr_loc + log_gro_frac)
  init_log_dec <- plogis(init_logit, log.p = TRUE)
  init_log_gro <- plogis(init_logit, lower.tail = FALSE, log.p = TRUE)

  for (v in seq_along(visit_weeks)) {
    dt <- visit_weeks[v] - 1  # weeks since baseline
    state_dec <- init_log_dec - dec_rate * dt
    state_gro <- init_log_gro + gro_rate * dt
    log_pred <- matrixStats::logSumExp(c(state_dec, state_gro))
    log_obs <- rnorm(1, mean = log_pred, sd = true_measure_sd)
    all_obs <- c(all_obs, exp(log_obs))  # normalized observation
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
  measure_sd = true_measure_sd
)

cat(sprintf("Data: %d patients, %d total visits\n",
            n_patients, length(all_obs)))

# --- Compile model ---
model <- cmdstan_model(
  "stan/experiments/laplace_psa_test.stan",
  quiet = FALSE
)

# --- Fit both modes ---
pop_params <- c("tr_loc_pop", "tr_sd", "frac_logit_pop", "frac_sd",
                "init_logit_pop", "init_sd")

# Good initial values near truth to avoid Newton NaN during warmup
init_fn <- function() {
  list(
    tr_loc_pop = rnorm(1, -3, 0.2),
    tr_sd = abs(rnorm(1, 0.3, 0.1)),
    frac_logit_pop = rnorm(1, 0.5, 0.2),
    frac_sd = abs(rnorm(1, 0.3, 0.1)),
    init_logit_pop = rnorm(1, 0, 0.2),
    init_sd = abs(rnorm(1, 0.3, 0.1))
  )
}

fit_modes <- list()
mode_names <- c("full_hmc", "laplace")
timings <- list()

for (mode in 0:1) {
  name <- mode_names[mode + 1]
  cat(sprintf("\n=== Fitting mode %d (%s) ===\n", mode, name))

  # For HMC mode, also init z params
  init_fn_mode <- if (mode == 0) {
    function() c(init_fn(), list(
      z_tr = rnorm(n_patients, 0, 0.5),
      z_frac = rnorm(n_patients, 0, 0.5),
      z_init = rnorm(n_patients, 0, 0.5)
    ))
  } else {
    init_fn
  }

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
    adapt_delta = 0.9
  )
  timings[[name]] <- (proc.time() - t0)["elapsed"]
  fit_modes[[name]] <- fit
}

# --- Compare posteriors ---
cat("\n\n=== Posterior comparison (population parameters) ===\n\n")

summaries <- map(mode_names, \(name) {
  fit_modes[[name]]$summary(variables = pop_params) |>
    mutate(mode = name, .before = 1)
}) |>
  list_rbind()

print(summaries |> select(mode, variable, mean, sd, q5, q95, rhat, ess_bulk))

# --- Check agreement ---
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
  # MCSE based on minimum ESS
  mcse = hmc_summary$sd / sqrt(pmin(hmc_summary$ess_bulk, lap_summary$ess_bulk)),
  diff_over_mcse = abs(hmc_summary$mean - lap_summary$mean) /
    (hmc_summary$sd / sqrt(pmin(hmc_summary$ess_bulk, lap_summary$ess_bulk)))
)

print(comparison)

cat(sprintf("\n\nTimings: HMC = %.1fs, Laplace = %.1fs\n",
            timings$full_hmc, timings$laplace))

# Add true values for reference
cat("\n=== True values ===\n")
cat(sprintf("tr_loc_pop:     %.2f\n", true_tr_loc_pop))
cat(sprintf("tr_sd:          %.2f\n", true_tr_sd))
cat(sprintf("frac_logit_pop: %.2f\n", true_frac_logit_pop))
cat(sprintf("frac_sd:        %.2f\n", true_frac_sd))
cat(sprintf("init_logit_pop: %.2f\n", true_init_logit_pop))
cat(sprintf("init_sd:        %.2f\n", true_init_sd))

max_ratio <- max(comparison$diff_over_mcse)
cat(sprintf("\nMax diff/MCSE ratio: %.2f\n", max_ratio))
if (max_ratio < 5) {
  cat("PASS: All population parameters agree within 5 MCSE.\n")
} else {
  cat("WARNING: Some parameters differ by more than 5 MCSE.\n")
}

cat("\nPhase 3 complete.\n")
