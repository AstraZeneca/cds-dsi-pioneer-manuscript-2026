# tests/testthat/test-stan-multistate-ic.R
library(testthat)
library(cmdstanr)
library(here)
source(here("tests/testthat/helper-stan.R"))

test_that("gap=0 IC likelihood matches no-IC likelihood", {
  # Single state-2 patient progressed then died, no gap
  # T_d=12, time_12=8, T_death=20
  log_surv_val <- -0.1

  data <- list(
    n_patients    = 1,
    max_t         = 30,
    max_sojourn_t = 20,
    final_state   = 2L,
    time_01       = 12L,
    time_02       = 0L,
    time_12       = 8L,
    time_03       = 0L,
    time_32       = 0L,
    censored_01   = 0L,
    prog_deterministic = 0L,
    ms_ic_gap_01  = 0L,     # no IC
    enable_03     = 0L,
    n_total_visits  = 1L,
    t_patient_visits = 1L,
    patient_visit_pos = c(1L, 2L),
    log_surv_val  = log_surv_val
  )

  fit <- test_stan_function(
    "tests/testthat/stan/test_multistate_ic.stan", data
  )
  draws <- fit$draws(format = "df")
  expect_equal(draws$ll_ic[1], draws$ll_no_ic[1], tolerance = 1e-8)
})

test_that("IC marginalization matches hand-computed expected value", {
  # State 2 patient: T_c=4, T_d=6 (gap=2), T_death=10
  # time_12 = T_death - T_d = 4
  # Constant log-conditional-survival = -0.1 for all transitions
  c_val <- -0.1
  c_h   <- log1p(-exp(c_val))  # log(1 - exp(-0.1)) = log(hazard prob)

  # Stan call uses enable_01=1, enable_02=1, enable_03=0.
  # base_surv: survived 01 and 02 risks weeks 1..T_c=4 (2 transitions × 4 weeks)
  base_surv <- 2 * 4 * c_val  # = -0.8

  # For each candidate s in {5, 6} (k=1,2 with T_c=4, T_d=6):
  # alpha(s): survived 0->1 from T_c+1 to s-1, event at s
  # beta(s):  survived 02 from T_c+1 to s-1 (enable_03=0, no 03 contribution)
  # gamma(s): sojourn event at tau_s = time_12 + (T_d - s)

  # k=1, s=5: alpha = c_h, beta = 0, tau_s = 4+1 = 5
  # sojourn event at 5: sum(lcs[1:4]) + log1m_exp(lcs[5]) = 4*c_val + c_h
  term1 <- c_h + 0 + (4 * c_val + c_h)

  # k=2, s=6: alpha = c_val + c_h, beta = c_val (only 02), tau_s = 4+0 = 4
  # sojourn event at 4: sum(lcs[1:3]) + log1m_exp(lcs[4]) = 3*c_val + c_h
  term2 <- (c_val + c_h) + (1 * c_val) + (3 * c_val + c_h)

  expected_ll <- base_surv + matrixStats::logSumExp(c(term1, term2))

  data <- list(
    n_patients    = 1L,
    max_t         = 30L,
    max_sojourn_t = 15L,
    final_state   = 2L,
    time_01       = 6L,
    time_02       = 0L,
    time_12       = 4L,
    time_03       = 0L,
    time_32       = 0L,
    censored_01   = 0L,
    prog_deterministic = 0L,
    ms_ic_gap_01  = 2L,
    enable_03     = 0L,
    n_total_visits  = 1L,
    t_patient_visits = 1L,
    patient_visit_pos = c(1L, 2L),
    log_surv_val  = c_val
  )

  fit <- test_stan_function(
    "tests/testthat/stan/test_multistate_ic.stan", data
  )
  draws <- fit$draws(format = "df")
  expect_equal(draws$ll_ic[1], expected_ll, tolerance = 1e-6)
})

test_that("IC disabled for deterministic progression (prog_deterministic=1)", {
  # When prog_deterministic=1, the IC gap must be ignored regardless of ms_ic_gap_01.
  # ll_ic (gap=2, prog_det=1) should equal ll_no_ic from the same fit
  # (ll_no_ic always zeroes gaps, so it gives the gap=0 path).
  # With the fix, both paths use T_c = T_d = 6 and produce identical likelihoods.
  c_val <- -0.1
  data_ic <- list(
    n_patients = 1L, max_t = 30L, max_sojourn_t = 15L,
    final_state = 2L, time_01 = 6L, time_02 = 0L, time_12 = 4L,
    time_03 = 0L, time_32 = 0L,
    censored_01 = 0L,
    prog_deterministic = 1L, ms_ic_gap_01 = 2L,
    enable_03 = 0L,
    n_total_visits  = 1L,
    t_patient_visits = 1L,
    patient_visit_pos = c(1L, 2L),
    log_surv_val = c_val
  )

  fit_ic <- test_stan_function("tests/testthat/stan/test_multistate_ic.stan", data_ic)
  draws  <- fit_ic$draws(format = "df")

  ll_ic   <- draws$ll_ic[1]
  ll_noop <- draws$ll_no_ic[1]
  expect_equal(ll_ic, ll_noop, tolerance = 1e-8)
})

test_that("State-0 (admin-censored) contributes 0->3 survival at visit weeks", {
  # A patient right-censored in state 0 at week 10, weekly visits
  # With enable_03, state-0 patients contribute 0→3 survival at all visit weeks
  # below time_01 + 1 = 11 (i.e., visits 1–10 are all included)
  c_val <- -0.1
  visits <- 1:10

  data <- list(
    n_patients = 1L, max_t = 15L, max_sojourn_t = 15L,
    final_state = 0L,
    time_01 = 10L, time_02 = 10L, time_12 = 0L,
    time_03 = 0L, time_32 = 0L,
    censored_01 = 1L,
    prog_deterministic = 0L,
    ms_ic_gap_01 = 0L,
    enable_03 = 1L,
    n_total_visits = length(visits),
    t_patient_visits = visits,
    patient_visit_pos = c(1L, length(visits) + 1L),
    log_surv_val = c_val
  )

  fit <- test_stan_function("tests/testthat/stan/test_multistate_ic.stan", data)
  draws <- fit$draws(format = "df")

  # Expected: 0→1 survival (10 weeks) + 0→2 survival (10 weeks)
  #         + 0→3 survival at 10 visit weeks (all below time_01 + 1 = 11)
  expected_ll <- 3 * 10 * c_val  # 0->1, 0->2, and 0->3
  expect_equal(draws$ll_ic[1], expected_ll, tolerance = 1e-6)
})

test_that("State-3 0->3 survival sums only at visit weeks", {
  # Patient drops out at week 12 (visit), visits at weeks 6 and 12
  # 0->3 hazard should accumulate only at week 6 (survival) and week 12 (event)
  c_val <- -0.1
  c_h   <- log1p(-exp(c_val))  # log(1 - exp(-0.1))
  visits <- c(6L, 12L)

  data <- list(
    n_patients = 1L, max_t = 20L, max_sojourn_t = 5L,
    final_state = 3L,
    time_01 = 0L, time_02 = 0L, time_12 = 0L,
    time_03 = 12L, time_32 = 0L,
    censored_01 = 1L,
    prog_deterministic = 0L,
    ms_ic_gap_01 = 0L,
    enable_03 = 1L,
    n_total_visits = length(visits),
    t_patient_visits = visits,
    patient_visit_pos = c(1L, length(visits) + 1L),
    log_surv_val = c_val
  )

  fit <- test_stan_function("tests/testthat/stan/test_multistate_ic.stan", data)
  draws <- fit$draws(format = "df")

  # 0->1 survival: weeks 1..11   = 11 * c_val
  # 0->2 survival: weeks 1..11   = 11 * c_val
  # 0->3 survival: only week 6   = 1  * c_val  (NOT all of weeks 1..11)
  # 0->3 event at week 12        = c_h
  expected_ll <- 11 * c_val + 11 * c_val + c_val + c_h
  expect_equal(draws$ll_ic[1], expected_ll, tolerance = 1e-6)
})

test_that("State-1 competing-risk 0->3 sums only at visit weeks before progression", {
  # Patient progresses at week 12, visits at weeks 6 and 12
  # 0->3 competing risk: only week 6 contributes (< 12)
  c_val <- -0.1
  c_h   <- log1p(-exp(c_val))
  visits <- c(6L, 12L)

  data <- list(
    n_patients = 1L, max_t = 20L, max_sojourn_t = 10L,
    final_state = 1L,
    time_01 = 12L, time_02 = 0L, time_12 = 5L,
    time_03 = 0L, time_32 = 0L,
    censored_01 = 0L,
    prog_deterministic = 0L,
    ms_ic_gap_01 = 0L,
    enable_03 = 1L,
    n_total_visits = length(visits),
    t_patient_visits = visits,
    patient_visit_pos = c(1L, length(visits) + 1L),
    log_surv_val = c_val
  )

  fit <- test_stan_function("tests/testthat/stan/test_multistate_ic.stan", data)
  draws <- fit$draws(format = "df")

  # gap=0 no-IC path, T_c = T_d = 12
  # 0->1: sum(lcs_01[1:T_c]) + event at T_d = 12*c_val + c_h
  # 0->2: sum(lcs_02[1:T_c])                = 12*c_val
  # 0->3: visits below T_d=12: only week 6  = 1*c_val
  # sojourn censored 1..5                   = 5*c_val
  expected_ll <- (12 * c_val + c_h) + 12 * c_val + c_val + 5 * c_val
  expect_equal(draws$ll_ic[1], expected_ll, tolerance = 1e-6)
})
