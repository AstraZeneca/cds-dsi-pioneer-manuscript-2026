library(testthat)
library(here)

# Helper: sentinel large time for disabled transitions
BIG <- 999L

make_ms_data <- function(n, max_t, event_time, censored, final_state,
                         lcs_01, lcs_02 = NULL, time_02 = NULL, censored_02 = NULL) {
  if (is.null(lcs_02))      lcs_02      <- matrix(-0.05, nrow = n, ncol = max_t)
  if (is.null(time_02))     time_02     <- rep(max_t, n)
  if (is.null(censored_02)) censored_02 <- rep(1L, n)

  list(
    N = n, MAX_T = max_t,
    event_time_01 = event_time,
    censored_01   = censored,
    log_cond_surv_01 = lcs_01,
    final_state = final_state,
    time_02     = time_02,
    censored_02 = censored_02,
    time_12     = rep(0L, n),
    censored_12 = rep(1L, n),
    log_cond_surv_02   = lcs_02,
    log_cond_surv_12_s = matrix(-0.05, n, max_t),
    log_cond_surv_12_t = matrix(-0.05, n, max_t),
    prog_deterministic = rep(0L, n),
    time_03     = rep(BIG, n),
    censored_32 = rep(1L, n),
    log_cond_surv_03  = matrix(-0.05, n, max_t),
    log_cond_surv_32  = matrix(-0.05, n, max_t)
  )
}

test_that("calc_ms_single_transition_loglik: detection-week convention", {
  n <- 3L; max_t <- 8L
  lcs <- matrix(-0.08, nrow = n, ncol = max_t)
  # Patient 1: event at week 4 (last_surv_week = 3)
  # Patient 2: censored at week 6 (last_surv_week = 6)
  # Patient 3: event at week 1 (last_surv_week = 0 -> only hazard)
  event_time <- c(4L, 6L, 1L)
  censored   <- c(0L, 1L, 0L)

  data <- make_ms_data(n, max_t, event_time, censored,
                       final_state = c(1L, 0L, 1L), lcs_01 = lcs)
  fit <- test_stan_function(
    here("tests", "testthat", "stan", "test_multistate_loglik_all.stan"),
    data
  )
  d <- posterior::as_draws_df(fit$draws())

  expected <- r_calc_ms_stl(event_time, censored, lcs)
  for (i in seq_len(n)) {
    expect_equal(get_stan_val(d, "st_llik", i), expected[i], tolerance = 1e-6,
      label = paste0("st_llik[", i, "]"))
  }
})

test_that("multistate_lpmf 01-only fast path == sum(calc_ms_single_transition_loglik)", {
  n <- 4L; max_t <- 6L
  lcs <- matrix(-0.1, nrow = n, ncol = max_t)
  event_time <- c(2L, 4L, 6L, 3L)
  censored   <- c(0L, 1L, 1L, 0L)

  data <- make_ms_data(n, max_t, event_time, censored,
                       final_state = c(1L, 0L, 0L, 1L), lcs_01 = lcs)
  fit <- test_stan_function(
    here("tests", "testthat", "stan", "test_multistate_loglik_all.stan"),
    data
  )
  d <- posterior::as_draws_df(fit$draws())

  expected_sum <- r_multistate_lpmf_01only(event_time, censored, lcs)
  expect_equal(get_stan_val(d, "ms_lpmf_01only"), expected_sum, tolerance = 1e-6)
  # Also verify fast path == sum(st_llik)
  stan_st_sum <- sum(sapply(seq_len(n), function(i) get_stan_val(d, "st_llik", i)))
  expect_equal(get_stan_val(d, "ms_lpmf_01only"), stan_st_sum, tolerance = 1e-6)
})

test_that("multistate_lpmf 01+02: state-0 censored patient accumulates both survival terms", {
  n <- 1L; max_t <- 5L
  lcs_01 <- matrix(-0.1, nrow = n, ncol = max_t)
  lcs_02 <- matrix(-0.05, nrow = n, ncol = max_t)
  cens_time <- 4L

  data <- make_ms_data(
    n, max_t,
    event_time = cens_time, censored = 1L,
    final_state = 0L,
    lcs_01 = lcs_01, lcs_02 = lcs_02,
    time_02 = cens_time, censored_02 = 1L
  )
  fit <- test_stan_function(
    here("tests", "testthat", "stan", "test_multistate_loglik_all.stan"),
    data
  )
  d <- posterior::as_draws_df(fit$draws())

  expected <- r_multistate_state0_cens(cens_time, lcs_01[1, ], lcs_02[1, ])
  expect_equal(get_stan_val(d, "ms_lpmf_01_02"), expected, tolerance = 1e-6)
})

test_that("multistate_lpmf 01+02: progressed patient (state 1) accumulates 01+02 survival + 01 hazard", {
  n <- 1L; max_t <- 8L
  lcs_01 <- matrix(-0.1, nrow = n, ncol = max_t)
  lcs_02 <- matrix(-0.06, nrow = n, ncol = max_t)
  t01 <- 5L  # progressed at week 5

  data <- make_ms_data(
    n, max_t,
    event_time = t01, censored = 0L,
    final_state = 1L,
    lcs_01 = lcs_01, lcs_02 = lcs_02,
    time_02 = BIG, censored_02 = 1L
  )
  # Set censored_12 = 1 (no post-progression death observed)
  data$time_12     <- 0L
  data$censored_12 <- 1L

  fit <- test_stan_function(
    here("tests", "testthat", "stan", "test_multistate_loglik_all.stan"),
    data
  )
  d <- posterior::as_draws_df(fit$draws())

  # State 1 contribution: 01 survival[1..t01-1] + 02 survival[1..t01-1] + 01 hazard at t01
  expected <- sum(lcs_01[1, 1:(t01 - 1)]) +
              sum(lcs_02[1, 1:(t01 - 1)]) +
              log1p(-exp(lcs_01[1, t01]))
  expect_equal(get_stan_val(d, "ms_lpmf_01_02"), expected, tolerance = 1e-6)
})
