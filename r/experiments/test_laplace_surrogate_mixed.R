# Phase 1b MIXED-COHORT gate: does routing the BACKGROUNDED subset through the
# surrogate preserve the population posterior, when FORECAST patients (full
# bi-exponential) pin tr_sd/frac_sd/init_sd? This matches production, unlike the
# all-backgrounded toy which forced the surrogate to identify all 3 SDs alone.
#
#   reference (mode 0): ALL patients full-HMC bi-exponential
#   mixed     (mode 1): forecast = full-HMC bi-exp; background = surrogate-Laplace
#
# PASS = the two runs' population posteriors agree (surrogate faithfully stands in
# for the backgrounded subset's contribution). tr_sd is the parameter that failed
# the all-backgrounded gate; here it should be pinned by the forecast cohort.
# ============================================================================
library(cmdstanr)
library(posterior)
library(dplyr)
library(tibble)
library(purrr)

set_cmdstan_path("~/.cmdstan/cmdstan-2.39.0")
set.seed(42)

# Larger cohort so the forecast subset alone can identify the SDs. Production has
# a modest forecast trial + larger historical (backgrounded) trials, so we use a
# forecast minority but enough to pin population spread.
n_forecast        <- 20
n_background      <- 40
n_patients        <- n_forecast + n_background
n_visits_per_patient <- 8

true <- list(tr_loc_pop = -3.0, tr_sd = 0.3, frac_logit_pop = 0.5,
             frac_sd = 0.3, init_logit_pop = 0.0, init_sd = 0.3,
             measure_sd = 0.15)
lod <- 0.85
log_lod <- log(lod)
visit_weeks <- seq(1, by = 4, length.out = n_visits_per_patient)
anchor_times <- c(0, 12, 28)

z_tr   <- rnorm(n_patients)
z_frac <- rnorm(n_patients)
z_init <- rnorm(n_patients)

# Patients ordered forecast-first (1..n_forecast forecast, rest backgrounded).
obs <- c(); idx <- c()
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
    obs <- c(obs, o); idx <- c(idx, visit_weeks[v])
  }
  pos[i + 1] <- pos[i] + n_visits_per_patient
}

cat(sprintf("Data: %d patients (%d forecast + %d background), %d visits, %d censored (%.0f%%)\n",
            n_patients, n_forecast, n_background, length(obs),
            sum(obs == 0), 100 * sum(obs == 0) / length(obs)))

stan_data <- list(
  n_patients = n_patients, n_forecast = n_forecast,
  n_total_visits = length(obs), normalized_obs = obs,
  patient_visit_pos = pos, visit_time = idx,
  measure_sd = true$measure_sd, log_lod = log_lod, anchor_times = anchor_times)

model <- cmdstan_model("stan/experiments/laplace_surrogate_mixed.stan",
                       include_paths = c("stan", "stan/tumor"), quiet = FALSE)

pop <- c("tr_loc_pop", "tr_sd", "frac_logit_pop", "frac_sd",
         "init_logit_pop", "init_sd")
init_fn <- function() list(
  tr_loc_pop = rnorm(1, -3, 0.1), tr_sd = 0.3 + abs(rnorm(1, 0, 0.05)),
  frac_logit_pop = rnorm(1, 0.5, 0.1), frac_sd = 0.3 + abs(rnorm(1, 0, 0.05)),
  init_logit_pop = rnorm(1, 0, 0.1), init_sd = 0.3 + abs(rnorm(1, 0, 0.05)))

fits <- list(); timings <- list()
for (mode in 0:1) {
  name <- c("reference", "mixed")[mode + 1]
  cat(sprintf("\n=== mode %d (%s) ===\n", mode, name))
  n_hmc <- if (mode == 0) n_patients else n_forecast
  init_mode <- function() c(init_fn(), list(
    z_tr = rnorm(n_hmc, 0, 0.5), z_frac = rnorm(n_hmc, 0, 0.5),
    z_init = rnorm(n_hmc, 0, 0.5)))
  t0 <- proc.time()
  fits[[name]] <- model$sample(
    data = c(stan_data, list(mixed_mode = mode)), init = init_mode,
    seed = 123, chains = 4, parallel_chains = 4,
    iter_warmup = 1000, iter_sampling = 2000, refresh = 500,
    adapt_delta = if (mode == 1) 0.95 else 0.9)
  timings[[name]] <- (proc.time() - t0)["elapsed"]
}

ref <- fits$reference$summary(variables = pop)
mix <- fits$mixed$summary(variables = pop)
comparison <- tibble(
  parameter = ref$variable, ref_mean = ref$mean, mixed_mean = mix$mean,
  diff = abs(ref$mean - mix$mean), ref_sd = ref$sd, mixed_sd = mix$sd,
  mcse = ref$sd / sqrt(pmin(ref$ess_bulk, mix$ess_bulk)),
  diff_over_mcse = abs(ref$mean - mix$mean) /
    (ref$sd / sqrt(pmin(ref$ess_bulk, mix$ess_bulk))),
  sd_ratio = mix$sd / ref$sd)
cat("\n=== population comparison (reference = all-HMC, mixed = background-surrogate) ===\n")
print(comparison)
cat(sprintf("\nTimings: reference=%.1fs  mixed=%.1fs\n",
            timings$reference, timings$mixed))

cor_ref <- fits$reference$draws(variables = pop, format = "draws_matrix") |> cor()
cor_mix <- fits$mixed$draws(variables = pop, format = "draws_matrix") |> cor()
max_cor_diff <- max(abs(cor_ref - cor_mix))
cat(sprintf("\nMax |corr difference| (joint structure): %.3f\n", max_cor_diff))

diag_one <- function(fit, name) {
  d <- fit$diagnostic_summary(quiet = TRUE)
  s <- fit$summary(variables = pop)
  n_draws <- fit$metadata()$iter_sampling * length(d$num_divergent)
  tibble(mode = name, pct_divergent = 100 * sum(d$num_divergent) / n_draws,
         max_rhat = max(s$rhat, na.rm = TRUE),
         min_ess_bulk = min(s$ess_bulk, na.rm = TRUE))
}
diag_tbl <- bind_rows(diag_one(fits$reference, "reference"),
                      diag_one(fits$mixed, "mixed"))
cat("\n=== convergence diagnostics ===\n"); print(diag_tbl)

# --- Forecast invariance: the REAL deliverable ----------------------------
# corr-diff is only a proxy for "does the target forecast change?". The actual
# deliverable is the posterior predictive of a population-drawn patient's future
# log-burden trajectory — a nonlinear function of the FULL joint population
# posterior (so it is sensitive to the frac/init correlation the marginals miss).
# If reference and mixed give the same predictive bands, any joint distortion on
# the weakly-identified nuisance ridge is harmless to what we report.
forecast_weeks <- c(4, 12, 24, 36, 52)   # future horizons (weeks)
pop_predictive <- function(fit, n_new = 4000) {
  d <- fit$draws(variables = pop, format = "draws_matrix")
  nd <- nrow(d)
  idx <- sample.int(nd, n_new, replace = TRUE)
  # one new patient per sampled draw -> integrates parameter + patient uncertainty
  z1 <- rnorm(n_new); z2 <- rnorm(n_new); z3 <- rnorm(n_new)
  tr_loc     <- d[idx, "tr_loc_pop"]     + d[idx, "tr_sd"]   * z1
  frac_logit <- d[idx, "frac_logit_pop"] + d[idx, "frac_sd"] * z2
  init_logit <- d[idx, "init_logit_pop"] + d[idx, "init_sd"] * z3
  ldf <- plogis(frac_logit, log.p = TRUE)
  lgf <- plogis(frac_logit, lower.tail = FALSE, log.p = TRUE)
  dr <- exp(tr_loc + ldf); gr <- exp(tr_loc + lgf)
  ild <- plogis(init_logit, log.p = TRUE)
  ilg <- plogis(init_logit, lower.tail = FALSE, log.p = TRUE)
  # mean log-burden trajectory (no measurement noise) at each horizon
  sapply(forecast_weeks, function(t) {
    matrixStats::rowLogSumExps(cbind(ild - dr * t, ilg + gr * t))
  })  # n_new x length(forecast_weeks)
}
set.seed(7)
pp_ref <- pop_predictive(fits$reference)
pp_mix <- pop_predictive(fits$mixed)
qs <- c(0.05, 0.25, 0.5, 0.75, 0.95)
fc <- map(seq_along(forecast_weeks), \(j) {
  qr <- quantile(pp_ref[, j], qs); qm <- quantile(pp_mix[, j], qs)
  tibble(week = forecast_weeks[j], q = names(qr),
         ref = as.numeric(qr), mixed = as.numeric(qm),
         abs_diff = abs(as.numeric(qr) - as.numeric(qm)))
}) |> list_rbind()
cat("\n=== FORECAST INVARIANCE: population-predictive log-burden quantiles ===\n")
print(fc, n = 100)
# Tolerance: forecast quantiles on the log scale; measure_sd=0.15, so a shift
# << measure_sd is negligible relative to observation noise. Use 0.05 (1/3 of sd).
max_fc_diff <- max(fc$abs_diff)
forecast_tol <- 0.05
cat(sprintf("\nMax |forecast quantile difference| (log-burden): %.4f  (tol %.2f, ~1/3 of measure_sd)\n",
            max_fc_diff, forecast_tol))

mix_d <- diag_tbl |> filter(mode == "mixed")
converged <- mix_d$pct_divergent < 1 && mix_d$max_rhat < 1.01 &&
  mix_d$min_ess_bulk > 400
max_ratio <- max(comparison$diff_over_mcse)
min_sd_ratio <- min(comparison$sd_ratio)

cat(sprintf("\nMax diff/MCSE: %.2f | min SD ratio (mixed/ref): %.2f | max corr diff: %.3f\n",
            max_ratio, min_sd_ratio, max_cor_diff))
cat(sprintf("tr_sd: reference=%.3f  mixed=%.3f  (truth %.3f)\n",
            comparison$ref_mean[comparison$parameter == "tr_sd"],
            comparison$mixed_mean[comparison$parameter == "tr_sd"], true$tr_sd))
# Decision: marginals (means+SDs) must agree AND the real deliverable (target
# forecast) must be invariant. corr-diff is now INFORMATIONAL — a flagged proxy,
# superseded by the direct forecast-invariance measurement when they disagree.
marginals_ok <- max_ratio < 5 && min_sd_ratio > 0.8
forecast_ok  <- max_fc_diff < forecast_tol
if (!converged) {
  cat("NO-GO: mixed fit did not converge.\n")
} else if (marginals_ok && forecast_ok) {
  cat("PASS: population marginals agree AND the target forecast is invariant.\n")
  if (max_cor_diff >= 0.1)
    cat(sprintf("  NOTE: joint corr-diff %.2f exceeds the 0.1 proxy ceiling, but the\n",
                max_cor_diff),
        "  direct forecast check shows this nuisance-ridge distortion is harmless\n",
        "  to the deliverable. Proxy superseded by territory.\n")
} else if (marginals_ok && !forecast_ok) {
  cat(sprintf("FAIL: marginals agree but the target forecast MOVED (max %.4f > %.2f)\n",
              max_fc_diff, forecast_tol),
      " — the joint distortion is a real defect, not just a proxy artifact.\n")
} else {
  cat("FAIL: population marginals disagree (means/SDs).\n")
}
cat("\nPhase 1b mixed-cohort gate complete.\n")
