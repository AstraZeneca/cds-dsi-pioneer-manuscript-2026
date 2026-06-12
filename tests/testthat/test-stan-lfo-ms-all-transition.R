library(testthat)
library(here)
library(stringr)
library(posterior)
source(here("tests/testthat/helper-stan.R"))

# Build a cohort with all transition types, plus a cutoff set beyond every
# event so re-censoring is a no-op and LFO log-lik must equal full log-lik.
#
# Patients 1-4 exercise the 0->1 / 0->2 / 1->2 transitions; patients 5-6 add
# the 0->3 dropout and 3->2 off-trial-death transitions so that turning on
# enable_03/enable_32 is non-vacuous. Assessment visits are at weeks 6, 12, 18
# for every patient; all visit-observed events (progression, dropout) fall on a
# visit week so the visit-gated 0->1 path lines up with the hazard accounting.
#
# `gated` toggles enable_ms_visit_gated_01 (+ enable_03/enable_32), letting the
# SAME cohort be run on either the continuous path or the production config.
make_cohort <- function(cutoff_week, gated = FALSE) {
  n <- 6L; max_t <- 30L
  # Distinct, mildly informative hazards per transition.
  lcs <- function(v) matrix(-v, nrow = n, ncol = max_t)
  list(
    N = n, MAX_T = max_t,
    # Patient 1: 0->1 progression at wk 12 (visit week), alive (state 1)
    # Patient 2: 0->2 direct death at wk 8 (state 2, censored_01 = 1)
    # Patient 3: 0->1 (wk 6, visit) -> 1->2 death, sojourn 5 (os event wk 11) (state 2)
    # Patient 4: admin-censored in state 0 at wk 20
    # Patient 5: 0->3 dropout at wk 12 (visit week), alive off-trial (state 3)
    # Patient 6: 0->3 dropout at wk 6 (visit week) -> 3->2 off-trial death, sojourn 4 (wk 10) (state 3)
    ms_final_state = c(1L, 2L, 2L, 0L, 3L, 3L),
    ms_time_01     = c(12L, 0L, 6L, 20L, 0L, 0L),
    ms_censored_01 = c(0L, 1L, 0L, 1L, 1L, 1L),
    ms_time_02     = c(12L, 8L, 6L, 20L, 0L, 0L),
    ms_time_12     = c(0L, 0L, 5L, 0L, 0L, 0L),
    ms_time_03     = c(0L, 0L, 0L, 0L, 12L, 6L),
    ms_time_32     = c(0L, 0L, 0L, 0L, 0L, 4L),
    ms_os_event_12 = c(0L, 0L, 11L, 0L, 0L, 0L),
    interval_censored     = c(0L, 0L, 0L, 0L, 0L, 0L),
    ms_prog_deterministic = c(1L, 0L, 1L, 0L, 0L, 0L),
    cutoff_visit_week = rep(cutoff_week, n),
    cutoff_cal_week   = rep(cutoff_week, n),
    patient_visit_pos = c(1L, 4L, 7L, 10L, 13L, 16L, 19L),  # 3 visits each
    t_patient_visits  = rep(c(6L, 12L, 18L), n),
    log_cond_surv_01 = lcs(0.10), log_cond_surv_02 = lcs(0.05),
    log_cond_surv_12_s = lcs(0.07), log_cond_surv_12_t = lcs(0.07),
    log_cond_surv_03 = lcs(0.02), log_cond_surv_32 = lcs(0.03),
    enable_01 = 1L, enable_02 = 1L, enable_12 = 1L, ms_time_scale_12 = 1L,
    enable_03 = if (gated) 1L else 0L,
    enable_32 = if (gated) 1L else 0L,
    enable_visit_gated_01 = if (gated) 1L else 0L
  )
}

run <- function(data) {
  fit <- test_stan_function(
    here("tests/testthat/stan/test_lfo_ms_all_transition.stan"), data)
  posterior::as_draws_df(fit$draws())
}

test_that("late cutoff: LFO all-transition log-lik == full multistate_lpmf", {
  d <- run(make_cohort(cutoff_week = 30L))  # beyond every event
  expect_equal(get_stan_val(d, "lfo_ll"), get_stan_val(d, "full_ll"),
               tolerance = 1e-6)
})

test_that("0->2 death is actually fit (early cutoff changes log-lik)", {
  late  <- run(make_cohort(cutoff_week = 30L))
  early <- run(make_cohort(cutoff_week = 7L))   # censors P1 progression, P3 death
  # If only 0->1 were fit, 0->2 hazard terms would be absent; with all
  # transitions, the re-censored cohort yields a different (valid) log-lik.
  expect_false(isTRUE(all.equal(get_stan_val(late, "lfo_ll"),
                                get_stan_val(early, "lfo_ll"),
                                tolerance = 1e-6)))
  expect_true(is.finite(get_stan_val(early, "lfo_ll")))
})

test_that("production config (visit-gated 0->1 + 0->3 + 3->2): LFO == full at late cutoff", {
  # Mirrors the publication OOS-rerun flags: enable_ms_02/12/03/32 = ON,
  # ms_time_scale_12 = 1, enable_ms_visit_gated_01 = 1. This exercises the
  # sum_at_visits_below visit-week branch of multistate_lpmf (a DIFFERENT code
  # path than the continuous tests above) and the 0->3 / 3->2 transitions.
  d <- run(make_cohort(cutoff_week = 30L, gated = TRUE))  # beyond every event
  lfo  <- get_stan_val(d, "lfo_ll")
  full <- get_stan_val(d, "full_ll")
  expect_true(is.finite(lfo))
  expect_true(is.finite(full))
  expect_equal(lfo, full, tolerance = 1e-6)
})

test_that("progression-only PFS censors death but operational PFS counts it", {
  # Guards the LFO OOS RECIST stamp against regressing to operational PFS.
  # Scenario: 2-patient RNG-free call to calculate_all_patients_endpoints_rng.
  #   - Patient 1: Direct death (0->2) at week 10 without progression.
  #     * Operational PFS: death IS a PFS event -> sample_right_censored[1] == 0
  #     * Progression-only PFS: no 0->1 event -> sample_prog_right_censored[1] == 1
  #       (the OOS RECIST stamp cannot write PD because progression is censored)
  #   - Patient 2: Observed 0->1 progression at week 12.
  #     * Both operational and progression-only PFS report an event (both censored==0).
  #
  # This asymmetry is the whole point of the fix: operational PFS includes death,
  # but the OOS stamp must NOT mark PD for a patient who died without progression.

  fit <- test_stan_function(
    here("tests/testthat/stan/test_oos_prog_only_all.stan"),
    list(dummy = 0L)
  )
  draws <- posterior::as_draws_df(fit$draws())

  # All assertions are inside the Stan harness; n_failures==0 means all passed
  expect_equal(get_stan_val(draws, "n_failures"), 0)
})
