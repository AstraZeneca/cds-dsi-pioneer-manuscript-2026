# Phase 1 gate: does the quadratic surrogate (Laplace, solver 1) reproduce the
# population posterior of the full bi-exponential (full HMC) for backgrounded
# patients? Data is simulated from the TRUE bi-exponential, so the surrogate is
# tested against the real generative process, not its own assumptions.
# ============================================================================
library(cmdstanr)
library(posterior)
library(dplyr)
library(tibble)
library(purrr)

set_cmdstan_path("~/.cmdstan/cmdstan-2.39.0")
set.seed(42)

n_patients <- 12
n_visits_per_patient <- 8

true <- list(tr_loc_pop = -3.0, tr_sd = 0.3, frac_logit_pop = 0.5,
             frac_sd = 0.3, init_logit_pop = 0.0, init_sd = 0.3,
             measure_sd = 0.15)
lod <- 0.85               # ~15% censored at the normalized-to-baseline scale
log_lod <- log(lod)
visit_weeks <- seq(1, by = 4, length.out = n_visits_per_patient)

# Fixed calendar anchors (weeks): baseline, mid, late — span the visit window.
anchor_times <- c(0, 12, 28)

z_tr   <- rnorm(n_patients)
z_frac <- rnorm(n_patients)
z_init <- rnorm(n_patients)

obs <- c(); idx <- c(); pid <- c()
pos <- integer(n_patients + 1); pos[1] <- 1L
for (i in seq_len(n_patients)) {
  tr_loc     <- true$tr_loc_pop     + true$tr_sd   * z_tr[i]
  frac_logit <- true$frac_logit_pop + true$frac_sd * z_frac[i]
  init_logit <- true$init_logit_pop + true$init_sd * z_init[i]
  ldf <- plogis(frac_logit, log.p = TRUE)
  lgf <- plogis(frac_logit, lower.tail = FALSE, log.p = TRUE)
  dr <- exp(tr_loc + ldf); gr <- exp(tr_loc + lgf)
  ild <- plogis(init_logit, log.p = TRUE)
  ilg <- plogis(init_logit, lower.tail = FALSE, log.p = TRUE)
  for (v in seq_along(visit_weeks)) {
    t <- visit_weeks[v] - 1
    lp <- matrixStats::logSumExp(c(ild - dr * t, ilg + gr * t))
    o <- exp(rnorm(1, lp, true$measure_sd))
    if (o < lod) o <- 0
    obs <- c(obs, o); idx <- c(idx, visit_weeks[v]); pid <- c(pid, i)
  }
  pos[i + 1] <- pos[i] + n_visits_per_patient
}

cat(sprintf("Data: %d patients, %d visits, %d censored (%.0f%%)\n",
            n_patients, length(obs), sum(obs == 0),
            100 * sum(obs == 0) / length(obs)))

stan_data <- list(
  n_patients = n_patients, n_total_visits = length(obs),
  normalized_obs = obs, patient_visit_pos = pos, visit_time = idx,
  patient_of_visit = pid, measure_sd = true$measure_sd, log_lod = log_lod,
  anchor_times = anchor_times)

model <- cmdstan_model("stan/experiments/laplace_surrogate_test.stan",
                       include_paths = c("stan", "stan/tumor"), quiet = FALSE)

pop <- c("tr_loc_pop", "tr_sd", "frac_logit_pop", "frac_sd",
         "init_logit_pop", "init_sd")
init_fn <- function() list(
  tr_loc_pop = rnorm(1, -3, 0.1), tr_sd = 0.3 + abs(rnorm(1, 0, 0.05)),
  frac_logit_pop = rnorm(1, 0.5, 0.1), frac_sd = 0.3 + abs(rnorm(1, 0, 0.05)),
  init_logit_pop = rnorm(1, 0, 0.1), init_sd = 0.3 + abs(rnorm(1, 0, 0.05)))

fits <- list(); timings <- list()
for (mode in 0:1) {
  name <- c("full_hmc", "surrogate")[mode + 1]
  cat(sprintf("\n=== mode %d (%s) ===\n", mode, name))
  init_mode <- if (mode == 0) {
    function() c(init_fn(), list(z_tr = rnorm(n_patients, 0, 0.5),
                                 z_frac = rnorm(n_patients, 0, 0.5),
                                 z_init = rnorm(n_patients, 0, 0.5)))
  } else init_fn
  t0 <- proc.time()
  fits[[name]] <- model$sample(
    data = c(stan_data, list(laplace_mode = mode)), init = init_mode,
    seed = 123, chains = 4, parallel_chains = 4,
    iter_warmup = 1000, iter_sampling = 2000, refresh = 500,
    adapt_delta = if (mode == 1) 0.95 else 0.9)
  timings[[name]] <- (proc.time() - t0)["elapsed"]
}

# --- Population comparison -------------------------------------------------
hmc <- fits$full_hmc$summary(variables = pop)
sur <- fits$surrogate$summary(variables = pop)
comparison <- tibble(
  parameter = hmc$variable, hmc_mean = hmc$mean, surrogate_mean = sur$mean,
  diff = abs(hmc$mean - sur$mean), hmc_sd = hmc$sd, surrogate_sd = sur$sd,
  mcse = hmc$sd / sqrt(pmin(hmc$ess_bulk, sur$ess_bulk)),
  diff_over_mcse = abs(hmc$mean - sur$mean) /
    (hmc$sd / sqrt(pmin(hmc$ess_bulk, sur$ess_bulk))),
  sd_ratio = sur$sd / hmc$sd)
cat("\n=== population comparison ===\n"); print(comparison)
cat(sprintf("\nTimings: HMC=%.1fs  surrogate=%.1fs\n",
            timings$full_hmc, timings$surrogate))

# --- Joint posterior: correlation-structure check --------------------------
# Subsumes target-forecast invariance (forecast is a nonlinear fn of the joint).
cor_hmc <- fits$full_hmc$draws(variables = pop, format = "draws_matrix") |>
  cor()
cor_sur <- fits$surrogate$draws(variables = pop, format = "draws_matrix") |>
  cor()
max_cor_diff <- max(abs(cor_hmc - cor_sur))
cat(sprintf("\nMax |corr difference| (joint structure): %.3f\n", max_cor_diff))

# --- Convergence gate (FIRST) ----------------------------------------------
diag_one <- function(fit, name) {
  d <- fit$diagnostic_summary(quiet = TRUE)
  s <- fit$summary(variables = pop)
  n_draws <- fit$metadata()$iter_sampling * length(d$num_divergent)
  tibble(mode = name, pct_divergent = 100 * sum(d$num_divergent) / n_draws,
         max_rhat = max(s$rhat, na.rm = TRUE),
         min_ess_bulk = min(s$ess_bulk, na.rm = TRUE))
}
diag_tbl <- bind_rows(diag_one(fits$full_hmc, "full_hmc"),
                      diag_one(fits$surrogate, "surrogate"))
cat("\n=== convergence diagnostics ===\n"); print(diag_tbl)

sur_d <- diag_tbl |> filter(mode == "surrogate")
converged <- sur_d$pct_divergent < 1 && sur_d$max_rhat < 1.01 &&
  sur_d$min_ess_bulk > 400
max_ratio <- max(comparison$diff_over_mcse)
min_sd_ratio <- min(comparison$sd_ratio)

cat(sprintf("\nMax diff/MCSE: %.2f | min SD ratio (sur/hmc): %.2f | max corr diff: %.3f\n",
            max_ratio, min_sd_ratio, max_cor_diff))
if (!converged) {
  cat("NO-GO: surrogate fit did not converge (solver 1 should be clean if",
      "the surrogate is log-concave as designed).\n")
} else if (max_ratio < 5 && min_sd_ratio > 0.8 && max_cor_diff < 0.1) {
  cat("PASS: surrogate converged AND population posterior (means, SDs, joint",
      "correlations) agrees with full HMC.\n")
} else {
  cat("FAIL: surrogate converged but population posterior differs",
      "(approximation bias) — see which of diff/MCSE, SD ratio, corr diff broke.\n")
}
cat("\nPhase 1 gate complete.\n")
