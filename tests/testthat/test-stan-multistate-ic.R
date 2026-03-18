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
    censored_01   = 0L,
    censored_02   = 1L,
    censored_12   = 0L,
    prog_deterministic = 0L,
    ms_ic_gap_01  = 0L,     # no IC
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
    censored_01   = 0L,
    censored_02   = 1L,
    censored_12   = 0L,
    prog_deterministic = 0L,
    ms_ic_gap_01  = 2L,
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
  data_ic <- list(n_patients=1L, max_t=30L, max_sojourn_t=15L,
    final_state=2L, time_01=6L, time_02=0L, time_12=4L,
    censored_01=0L, censored_02=1L, censored_12=0L,
    prog_deterministic=1L, ms_ic_gap_01=2L, log_surv_val=c_val)

  fit_ic <- test_stan_function("tests/testthat/stan/test_multistate_ic.stan", data_ic)
  draws  <- fit_ic$draws(format="df")

  ll_ic   <- draws$ll_ic[1]
  ll_noop <- draws$ll_no_ic[1]
  expect_equal(ll_ic, ll_noop, tolerance = 1e-8)
})
