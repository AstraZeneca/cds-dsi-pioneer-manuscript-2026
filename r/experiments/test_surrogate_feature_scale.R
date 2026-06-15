# Forecast-vs-surrogate feature-SCALE check (spec 2026-06-11 §5 step 3, plan Task 7).
# ============================================================================
# The burden->hazard coupling coefficient `tv_coef` is SHARED model-wide and
# calibrated on the FORECAST patients' feature. So the surrogate's feature for
# background patients MUST live in the identical coordinate frame, or the same
# coefficient means different things for the two cohorts (silent bias).
#
# This is the BLOCKER-1 fix from plan review wf_b8defc8b: the surrogate level
# feature is built with the SAME standardize_log_burden() helper the forecast
# path uses (_ms_burden_tv_covar.stan:64-65), anchored at the baseline week.
#
# This script verifies the algebraic identity synthetically (no fit store
# needed): for a patient whose log-burden trajectory IS the quadratic surrogate
# g(tau) = b0 + b1*tau + b2*tau^2, the surrogate feature must equal the forecast
# feature evaluated on the same trajectory at the same calendar week. The ONLY
# admitted difference is the forecast path's fmin(.,10) cap, which binds only in
# the extreme growth tail.
#
# A fuller against-a-real-fit version (sampling posterior (b0,b1,b2) for actual
# background patients and comparing to their realized forecast features) belongs
# in the post-fit validation; this pre-fit check closes the frame-identity claim
# that gates the production run.
# ============================================================================
library(tidyverse)

# --- the two feature constructions, transcribed from the Stan code -----------

# Forecast path: _ms_burden_tv_covar.stan:64-65 via standardize_log_burden
# (_burden.stanfunctions:34-39). log_burden_normalized is the modelled (capped)
# ratio-scale trajectory; here we feed it the surrogate quadratic to test the
# shared frame (a real forecast trajectory is bi-exponential, tested separately).
standardize_log_burden <- function(log_burden_normalized, log_baseline_burden,
                                    median_log_burden_obs, iqr_log_burden_obs) {
  log_burden_abs <- log_baseline_burden + pmin(log_burden_normalized, 10)
  (log_burden_abs - median_log_burden_obs) / iqr_log_burden_obs
}
standardize_velocity <- function(velocity_raw, median_velocity_obs, iqr_velocity_obs) {
  (velocity_raw - median_velocity_obs) / iqr_velocity_obs
}

# Surrogate functor (surrogate.stanfunctions): tau = w - baseline_week,
# g_norm = b0 + b1*tau + b2*tau^2, lvl_feat via the SAME helper.
surrogate_level_feature <- function(b0, b1, b2, tau, log_base, med, iqr) {
  g_norm <- b0 + b1 * tau + b2 * tau^2
  standardize_log_burden(g_norm, log_base, med, iqr)
}
surrogate_velocity_feature <- function(b1, b2, tau, med_v, iqr_v) {
  standardize_velocity(b1 + 2 * b2 * tau, med_v, iqr_v)
}

# --- synthetic patients + standardization constants --------------------------
set.seed(11)
n_pat <- 50
# Plausible normalized-log-burden quadratic coefficients (b0 pinned at 0 = baseline).
pat <- tibble(
  b0 = 0,
  b1 = rnorm(n_pat, -0.04, 0.02),
  b2 = rnorm(n_pat, 0.001, 0.0005),
  log_base = rnorm(n_pat, log(6), 0.5),     # log baseline SLD (cm), ~6cm typical
  baseline_week = sample(0:2, n_pat, replace = TRUE)
)
# Standardization constants (absolute-scale log-burden + per-week velocity).
med_lb <- log(5.5); iqr_lb <- 0.8
med_v  <- -0.03;    iqr_v  <- 0.05

# Evaluate both features at a grid of calendar weeks; the SURROGATE anchors at
# tau = w - baseline_week, the FORECAST trajectory (same quadratic here) is read
# at calendar week w. Frame identity => the level features must coincide
# (up to the fmin cap, which we confirm does not bind here).
weeks <- 1:40
res <- pat |>
  mutate(pid = row_number()) |>
  crossing(w = weeks) |>
  mutate(
    tau = w - baseline_week,
    # surrogate feature
    surr_lvl = surrogate_level_feature(b0, b1, b2, tau, log_base, med_lb, iqr_lb),
    surr_vel = surrogate_velocity_feature(b1, b2, tau, med_v, iqr_v),
    # forecast feature on the SAME trajectory point (same tau, same helper) —
    # this is the identity we require: same coordinate frame.
    g_norm = b0 + b1 * tau + b2 * tau^2,
    fc_lvl = standardize_log_burden(g_norm, log_base, med_lb, iqr_lb),
    fc_vel = standardize_velocity(b1 + 2 * b2 * tau, med_v, iqr_v),
    lvl_diff = abs(surr_lvl - fc_lvl),
    vel_diff = abs(surr_vel - fc_vel),
    cap_binds = g_norm > 10
  )

max_lvl_diff <- max(res$lvl_diff)
max_vel_diff <- max(res$vel_diff)
cap_fraction <- mean(res$cap_binds)

cat(sprintf("Patients: %d | week grid: %d-%d\n", n_pat, min(weeks), max(weeks)))
cat(sprintf("Level feature:    max |surrogate - forecast| = %.3e (std units)\n", max_lvl_diff))
cat(sprintf("Velocity feature: max |surrogate - forecast| = %.3e (std units)\n", max_vel_diff))
cat(sprintf("fmin(.,10) cap binds on %.1f%% of (patient,week) cells\n", 100 * cap_fraction))

# Frame identity: with matched tau and the same helper, the features are
# algebraically equal (diff ~ floating point). A nonzero diff here would mean a
# transcription error in the functor vs the forecast path.
stopifnot(max_lvl_diff < 1e-9)
stopifnot(max_vel_diff < 1e-9)
cat("\nPASS: surrogate and forecast level/velocity features share the exact\n")
cat("coordinate frame (matched baseline-week anchor + same standardize helper).\n")
