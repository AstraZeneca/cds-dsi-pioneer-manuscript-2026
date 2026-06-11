# Phase 1d gate: JOINT SLD + multistate(0->1, 0->3) surrogate WITH CORRELATED
# PATIENT-LEVEL FRAILTY (contract §9 step 2).
# Validates that marginalizing background patients' (b1,b2,u01,u03) recovers the
# population params, the (level, velocity) coupling coefficients on BOTH hazards,
# the frailty scales sigma_01/03, AND the 01<->03 frailty correlation rho — i.e.
# the full publication design, not just burden coupling.
#
# reference (mode 0): full HMC over per-patient (b1,b2) + correlated (u01,u03)
# joint     (mode 1): joint-Laplace marginalization (d=4 per patient)
# PASS: pop params + coefs + sigmas + rho agree, surrogate converges clean (solver 1).
# ============================================================================
library(cmdstanr); library(posterior); library(dplyr); library(tibble); library(purrr)
set_cmdstan_path("~/.cmdstan/cmdstan-2.39.0")
set.seed(11)

# Gate-design note: the frailty variances (s01,s03) and their correlation (rho)
# are only identifiable with ENOUGH events on both transitions. The first run
# (base ~ -4, ~18-32% event rates, 7+13 events) left s01/s03/rho unidentified —
# HMC itself missed truth, so the Laplace-vs-HMC comparison was meaningless.
# Fix: raise the baseline hazards so most patients have events on BOTH
# transitions (well-identified frailty), and keep n_patients modest + visits
# weekly so runtime stays reasonable (the d=4 HMC reference is the cost driver).
n_patients <- 45                   # well-identified (82%/64% event rates) yet ~2.4h not 4.3h
visit_t <- 0:28                    # WEEKLY 0->1 gating visits (was every 4w) -> more 0->1 info
n_wk <- 28
true <- list(
  b1_pop = -0.04, b2_pop = 0.001, b1_sd = 0.03, b2_sd = 0.0015,
  base_01 = -2.6, cf_lvl_01 = 0.8, cf_vel_01 = 0.5,   # raised base -> ~70% event rate
  base_03 = -2.9, cf_lvl_03 = 0.4, cf_vel_03 = 0.3,   # raised base -> ~60% event rate
  s01 = 0.35, s03 = 0.30, rho = -0.5,            # NEGATIVE corr (production motivation)
  measure_sd = 0.15)

# Correlated frailty draws: u_i = diag(s) %*% L %*% z_i
Lcorr <- t(chol(matrix(c(1, true$rho, true$rho, 1), 2)))   # lower Cholesky of corr
S <- diag(c(true$s01, true$s03))

# Per-patient latent draws (burden REs + correlated frailty), one clean pass.
z_b1 <- rnorm(n_patients); z_b2 <- rnorm(n_patients)
U <- matrix(0, n_patients, 2)               # [i, ] = (u01, u03)
for (i in seq_len(n_patients)) U[i, ] <- as.vector(S %*% Lcorr %*% rnorm(2))
b1_i <- true$b1_pop + true$b1_sd * z_b1
b2_i <- true$b2_pop + true$b2_sd * z_b2

# SLD observations + 0->1 visit-gating mask.
sld_obs <- c(); sld_time <- c(); sld_is_visit <- c()
sld_pos <- integer(n_patients + 1); sld_pos[1] <- 1L
visit_wk <- matrix(0L, n_patients, n_wk)
for (i in seq_len(n_patients)) {
  for (t in visit_t) {
    sld_obs <- c(sld_obs, rnorm(1, b1_i[i] * t + b2_i[i] * t^2, true$measure_sd))
    sld_time <- c(sld_time, t)
    sld_is_visit <- c(sld_is_visit, 1L)
    if (t >= 1 && t <= n_wk) visit_wk[i, t] <- 1L
  }
  sld_pos[i + 1] <- sld_pos[i] + length(visit_t)
}

# Empirical standardization constants (data-derived, shared by both modes —
# mirrors production transformed-data constants). Computed over all patient-weeks.
lvl_all <- as.vector(outer(b1_i, 1:n_wk) + outer(b2_i, (1:n_wk)^2))
vel_all <- as.vector(outer(b1_i, rep(1, n_wk)) + outer(2 * b2_i, 1:n_wk))
median_lvl <- median(lvl_all); iqr_lvl <- IQR(lvl_all); if (iqr_lvl == 0) iqr_lvl <- 1
median_vel <- median(vel_all); iqr_vel <- IQR(vel_all); if (iqr_vel == 0) iqr_vel <- 1
std_lvl <- function(x) (x - median_lvl) / iqr_lvl
std_vel <- function(x) (x - median_vel) / iqr_vel

# Events (0->1 visit-gated, 0->3 continuous) using the SAME standardization.
ms_event_wk_01 <- integer(n_patients); ms_censored_01 <- integer(n_patients)
ms_event_wk_03 <- integer(n_patients); ms_censored_03 <- integer(n_patients)
for (i in seq_len(n_patients)) {
  b1 <- b1_i[i]; b2 <- b2_i[i]
  u01 <- U[i, 1]; u03 <- U[i, 2]

  # 0->1 visit-gated event
  ev1 <- 0L
  for (w in 1:n_wk) {
    if (visit_wk[i, w] == 1L) {
      lvl <- b1 * w + b2 * w^2; vel <- b1 + 2 * b2 * w
      haz <- exp(true$base_01 + u01 + true$cf_lvl_01 * std_lvl(lvl) + true$cf_vel_01 * std_vel(vel))
      if (runif(1) < 1 - exp(-haz)) { ev1 <- w; break }
    }
  }
  if (ev1 == 0L) { ms_event_wk_01[i] <- n_wk; ms_censored_01[i] <- 1L }
  else { ms_event_wk_01[i] <- ev1; ms_censored_01[i] <- 0L }

  # 0->3 continuous event
  ev3 <- 0L
  for (w in 1:n_wk) {
    lvl <- b1 * w + b2 * w^2; vel <- b1 + 2 * b2 * w
    haz <- exp(true$base_03 + u03 + true$cf_lvl_03 * std_lvl(lvl) + true$cf_vel_03 * std_vel(vel))
    if (runif(1) < 1 - exp(-haz)) { ev3 <- w; break }
  }
  if (ev3 == 0L) { ms_event_wk_03[i] <- n_wk; ms_censored_03[i] <- 1L }
  else { ms_event_wk_03[i] <- ev3; ms_censored_03[i] <- 0L }
}

cat(sprintf("Data: %d patients, %d SLD obs | 0->1 events %d (%.0f%% cens) | 0->3 events %d (%.0f%% cens)\n",
            n_patients, length(sld_obs),
            sum(ms_censored_01 == 0), 100 * mean(ms_censored_01),
            sum(ms_censored_03 == 0), 100 * mean(ms_censored_03)))

stan_data <- list(
  n_patients = n_patients, n_total_sld = length(sld_obs), n_wk = n_wk,
  sld_obs = sld_obs, sld_pos = sld_pos, sld_time = sld_time, sld_is_visit = sld_is_visit,
  ms_event_wk_01 = ms_event_wk_01, ms_censored_01 = ms_censored_01,
  ms_event_wk_03 = ms_event_wk_03, ms_censored_03 = ms_censored_03,
  visit_wk = visit_wk,
  median_lvl = median_lvl, iqr_lvl = iqr_lvl, median_vel = median_vel, iqr_vel = iqr_vel,
  measure_sd = true$measure_sd)

model <- cmdstan_model("stan/experiments/laplace_joint_frailty_test.stan",
                       include_paths = c("stan", "stan/tumor"), quiet = FALSE)

pop <- c("b1_pop", "b2_pop", "b1_sd", "b2_sd",
         "base_01", "cf_lvl_01", "cf_vel_01",
         "base_03", "cf_lvl_03", "cf_vel_03",
         "s01", "s03", "rho_frailty")

init_fn <- function() list(
  b1_pop = rnorm(1, -0.04, 0.01), b2_pop = rnorm(1, 0.001, 0.0005),
  b1_sd = 0.03 + abs(rnorm(1, 0, 0.005)), b2_sd = 0.0015 + abs(rnorm(1, 0, 0.0003)),
  base_01 = rnorm(1, -4, 0.2), cf_lvl_01 = rnorm(1, 0.8, 0.1), cf_vel_01 = rnorm(1, 0.5, 0.1),
  base_03 = rnorm(1, -4.3, 0.2), cf_lvl_03 = rnorm(1, 0.4, 0.1), cf_vel_03 = rnorm(1, 0.3, 0.1),
  s01 = 0.35 + abs(rnorm(1, 0, 0.05)), s03 = 0.30 + abs(rnorm(1, 0, 0.05)),
  L_frailty = t(chol(matrix(c(1, -0.3, -0.3, 1), 2))))

fits <- list()
for (mode in 0:1) {
  name <- c("reference", "joint")[mode + 1]
  cat(sprintf("\n=== mode %d (%s) ===\n", mode, name))
  init_mode <- if (mode == 0)
    function() c(init_fn(), list(z_b1 = rnorm(n_patients, 0, 0.5),
                                 z_b2 = rnorm(n_patients, 0, 0.5),
                                 z_frailty = matrix(rnorm(2 * n_patients, 0, 0.5), 2, n_patients)))
    else init_fn
  fits[[name]] <- model$sample(
    data = c(stan_data, list(laplace_mode = mode)), init = init_mode,
    seed = 123, chains = 4, parallel_chains = 4,
    iter_warmup = 1000, iter_sampling = 1000, refresh = 500,
    adapt_delta = if (mode == 1) 0.95 else 0.9)
}

ref <- fits$reference$summary(variables = pop)
jnt <- fits$joint$summary(variables = pop)
comparison <- tibble(
  parameter = ref$variable, true = unlist(true[c(
    "b1_pop","b2_pop","b1_sd","b2_sd","base_01","cf_lvl_01","cf_vel_01",
    "base_03","cf_lvl_03","cf_vel_03","s01","s03","rho")]),
  ref_mean = ref$mean, joint_mean = jnt$mean,
  diff = abs(ref$mean - jnt$mean), ref_sd = ref$sd, joint_sd = jnt$sd,
  diff_over_mcse = abs(ref$mean - jnt$mean) / (ref$sd / sqrt(pmin(ref$ess_bulk, jnt$ess_bulk))),
  sd_ratio = jnt$sd / ref$sd)
cat("\n=== population comparison (reference HMC vs joint-Laplace) ===\n")
print(comparison, n = nrow(comparison))

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
key <- comparison |> filter(parameter %in% c("cf_lvl_01","cf_vel_01","cf_lvl_03","cf_vel_03","s01","s03","rho_frailty"))
cat(sprintf("\nMax diff/MCSE: %.2f | min SD ratio: %.2f\n", max_ratio, min_sd_ratio))
cat("Key coupling/frailty params (HMC vs joint):\n"); print(key |> select(parameter, true, ref_mean, joint_mean, diff_over_mcse))
if (!converged) {
  cat("\nNO-GO: joint-Laplace did not converge clean (solver 1 should hold if log-concave).\n")
} else if (max_ratio < 5 && min_sd_ratio > 0.8) {
  cat("\nPASS: joint SLD+MS+frailty surrogate recovers ALL params (coefs, sigmas, rho) vs full HMC.\n")
} else {
  cat("\nFAIL: joint surrogate converged but a param differs (check the diff_over_mcse column).\n")
}
cat("\nPhase 1d joint-frailty surrogate gate complete.\n")
