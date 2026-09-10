# Phase 1c gate: JOINT SLD + multistate(0->1) surrogate.
# Validates that marginalizing background patients' (b1,b2) from the JOINT
# likelihood (SLD + burden-coupled 0->1 hazard) recovers the population params
# AND the coupling coefficient coef_01 — the quantity Option A would lose.
#
# reference (mode 0): full HMC over per-patient (b1,b2)
# joint     (mode 1): joint-Laplace marginalization
# PASS: pop params + coef_01 agree, surrogate converges clean (solver 1).
# ============================================================================
library(cmdstanr); library(posterior); library(dplyr); library(tibble); library(purrr)
set_cmdstan_path("~/.cmdstan/cmdstan-2.39.0")
set.seed(7)

n_patients <- 30
visit_t <- seq(0, 28, by = 4)
n_wk <- 28
true <- list(b1_pop = -0.04, b2_pop = 0.001, b1_sd = 0.03, b2_sd = 0.0015,
             base_01 = -4.0, coef_01 = 0.8, measure_sd = 0.15)
median_lb <- 0.0; iqr_lb <- 1.0

z_b1 <- rnorm(n_patients); z_b2 <- rnorm(n_patients)
sld_obs <- c(); sld_time <- c(); sld_pos <- integer(n_patients + 1); sld_pos[1] <- 1L
ms_event_wk <- integer(n_patients); ms_censored <- integer(n_patients)
for (i in seq_len(n_patients)) {
  b1 <- true$b1_pop + true$b1_sd * z_b1[i]
  b2 <- true$b2_pop + true$b2_sd * z_b2[i]
  for (t in visit_t) {
    sld_obs <- c(sld_obs, rnorm(1, b1 * t + b2 * t^2, true$measure_sd))
    sld_time <- c(sld_time, t)
  }
  sld_pos[i + 1] <- sld_pos[i] + length(visit_t)
  # 0->1 event from burden-coupled hazard
  ev <- 0L
  for (w in 1:n_wk) {
    muw <- b1 * w + b2 * w^2
    haz <- exp(true$base_01 + true$coef_01 * ((muw - median_lb) / iqr_lb))
    if (runif(1) < 1 - exp(-haz)) { ev <- w; break }
  }
  if (ev == 0L) { ms_event_wk[i] <- n_wk; ms_censored[i] <- 1L }
  else { ms_event_wk[i] <- ev; ms_censored[i] <- 0L }
}
cat(sprintf("Data: %d patients, %d SLD obs, %d 0->1 events (%.0f%% censored)\n",
            n_patients, length(sld_obs), sum(ms_censored == 0),
            100 * mean(ms_censored)))

stan_data <- list(
  n_patients = n_patients, n_total_sld = length(sld_obs), n_wk = n_wk,
  sld_obs = sld_obs, sld_pos = sld_pos, sld_time = sld_time,
  ms_event_wk = ms_event_wk, ms_censored = ms_censored,
  median_lb = median_lb, iqr_lb = iqr_lb, measure_sd = true$measure_sd)

model <- cmdstan_model("stan/experiments/laplace_joint_test.stan",
                       include_paths = c("stan", "stan/tumor"), quiet = FALSE)

pop <- c("b1_pop", "b2_pop", "b1_sd", "b2_sd", "base_01", "coef_01")
init_fn <- function() list(
  b1_pop = rnorm(1, -0.04, 0.01), b2_pop = rnorm(1, 0.001, 0.0005),
  b1_sd = 0.03 + abs(rnorm(1, 0, 0.005)), b2_sd = 0.0015 + abs(rnorm(1, 0, 0.0003)),
  base_01 = rnorm(1, -4, 0.2), coef_01 = rnorm(1, 0.8, 0.1))

fits <- list()
for (mode in 0:1) {
  name <- c("reference", "joint")[mode + 1]
  cat(sprintf("\n=== mode %d (%s) ===\n", mode, name))
  init_mode <- if (mode == 0)
    function() c(init_fn(), list(z_b1 = rnorm(n_patients, 0, 0.5), z_b2 = rnorm(n_patients, 0, 0.5)))
    else init_fn
  fits[[name]] <- model$sample(
    data = c(stan_data, list(laplace_mode = mode)), init = init_mode,
    seed = 123, chains = 4, parallel_chains = 4,
    iter_warmup = 1000, iter_sampling = 2000, refresh = 500,
    adapt_delta = if (mode == 1) 0.95 else 0.9)
}

ref <- fits$reference$summary(variables = pop)
jnt <- fits$joint$summary(variables = pop)
comparison <- tibble(
  parameter = ref$variable, ref_mean = ref$mean, joint_mean = jnt$mean,
  diff = abs(ref$mean - jnt$mean), ref_sd = ref$sd, joint_sd = jnt$sd,
  mcse = ref$sd / sqrt(pmin(ref$ess_bulk, jnt$ess_bulk)),
  diff_over_mcse = abs(ref$mean - jnt$mean) / (ref$sd / sqrt(pmin(ref$ess_bulk, jnt$ess_bulk))),
  sd_ratio = jnt$sd / ref$sd)
cat("\n=== population comparison (reference HMC vs joint-Laplace) ===\n")
print(comparison)
cat(sprintf("\ncoef_01 (the coupling): true %.3f | ref %.3f | joint %.3f\n",
            true$coef_01,
            comparison$ref_mean[comparison$parameter == "coef_01"],
            comparison$joint_mean[comparison$parameter == "coef_01"]))

diag_one <- function(fit, name) {
  d <- fit$diagnostic_summary(quiet = TRUE); s <- fit$summary(variables = pop)
  n_draws <- fit$metadata()$iter_sampling * length(d$num_divergent)
  tibble(mode = name, pct_divergent = 100 * sum(d$num_divergent) / n_draws,
         max_rhat = max(s$rhat, na.rm = TRUE), min_ess_bulk = min(s$ess_bulk, na.rm = TRUE))
}
diag_tbl <- bind_rows(diag_one(fits$reference, "reference"), diag_one(fits$joint, "joint"))
cat("\n=== convergence diagnostics ===\n"); print(diag_tbl)

jd <- diag_tbl |> filter(mode == "joint")
converged <- jd$pct_divergent < 1 && jd$max_rhat < 1.01 && jd$min_ess_bulk > 400
max_ratio <- max(comparison$diff_over_mcse)
min_sd_ratio <- min(comparison$sd_ratio)
coef_ratio <- comparison$diff_over_mcse[comparison$parameter == "coef_01"]
cat(sprintf("\nMax diff/MCSE: %.2f | min SD ratio: %.2f | coef_01 diff/MCSE: %.2f\n",
            max_ratio, min_sd_ratio, coef_ratio))
if (!converged) {
  cat("NO-GO: joint-Laplace did not converge (solver 1 should be clean if log-concave).\n")
} else if (max_ratio < 5 && min_sd_ratio > 0.8) {
  cat("PASS: joint SLD+MS surrogate recovers population params AND coef_01 vs full HMC.\n")
} else {
  cat("FAIL: joint surrogate converged but a population param (check coef_01) differs.\n")
}
cat("\nPhase 1c joint-surrogate gate complete.\n")
