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
    ms_ic_gap_01       = rep(1L, n),  # 1-week gap: T_c = T_d - 1, matches detection-week convention
    # Dummy visit arrays — enable_03=0 in all tests, so these are never accessed
    N_visits           = n,
    t_patient_visits   = rep(1L, n),
    patient_visit_pos  = seq_len(n + 1L),
    time_03     = rep(BIG, n),
    censored_32 = rep(1L, n),
    log_cond_surv_03  = matrix(-0.05, n, max_t),
    log_cond_surv_32  = matrix(-0.05, n, max_t),
    weight      = rep(1.0, n)
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

test_that("multistate_lpmf: weight=0.5 halves patient contribution", {
  n <- 2L; max_t <- 8L
  lcs <- matrix(-0.08, nrow = n, ncol = max_t)
  data_full <- make_ms_data(n, max_t, c(4L, 6L), c(0L, 1L),
                            final_state = c(1L, 0L), lcs_01 = lcs)
  data_half <- modifyList(data_full, list(weight = c(0.5, 1.0)))

  fit_full <- test_stan_function(
    here("tests", "testthat", "stan", "test_multistate_loglik_all.stan"),
    data_full
  )
  fit_half <- test_stan_function(
    here("tests", "testthat", "stan", "test_multistate_loglik_all.stan"),
    data_half
  )

  d_full <- posterior::as_draws_df(fit_full$draws())
  d_half <- posterior::as_draws_df(fit_half$draws())

  full_ll  <- get_stan_val(d_full, "ms_lpmf_01only")
  half_ll  <- get_stan_val(d_half, "ms_lpmf_01only")
  p1_ll    <- get_stan_val(d_full, "st_llik", 1)
  expect_equal(full_ll - half_ll, 0.5 * p1_ll, tolerance = 1e-6)
})

test_that("multistate_lpmf: weight=0.5 halves patient contribution in full illness-death loop", {
  n <- 2L; max_t <- 8L
  lcs <- matrix(-0.08, nrow = n, ncol = max_t)
  # Both patients progress (state 1): exercises the full multi-transition loop
  data_full <- make_ms_data(n, max_t,
    event_time = c(4L, 6L), censored = c(0L, 0L),
    final_state = c(1L, 1L), lcs_01 = lcs,
    lcs_02 = matrix(-0.02, nrow = n, ncol = max_t),
    time_02 = c(8L, 8L), censored_02 = c(1L, 1L)
  )
  data_half <- modifyList(data_full, list(weight = c(0.5, 1.0)))

  fit_full <- test_stan_function(
    here("tests", "testthat", "stan", "test_multistate_loglik_all.stan"),
    data_full
  )
  fit_half <- test_stan_function(
    here("tests", "testthat", "stan", "test_multistate_loglik_all.stan"),
    data_half
  )

  d_full <- posterior::as_draws_df(fit_full$draws())
  d_half <- posterior::as_draws_df(fit_half$draws())

  # ms_lpmf_01_02 uses the general loop (both 0->1 and 0->2 enabled)
  full_ll <- get_stan_val(d_full, "ms_lpmf_01_02")
  half_ll <- get_stan_val(d_half, "ms_lpmf_01_02")

  # The difference should equal 0.5 * patient_1's contribution
  # Compute patient 1's contribution with weight=1 by running with weight=(1,0)
  data_p1_only <- modifyList(data_full, list(weight = c(1.0, 0.0)))
  fit_p1 <- test_stan_function(
    here("tests", "testthat", "stan", "test_multistate_loglik_all.stan"),
    data_p1_only
  )
  d_p1 <- posterior::as_draws_df(fit_p1$draws())
  p1_only_ll <- get_stan_val(d_p1, "ms_lpmf_01_02")

  expect_equal(full_ll - half_ll, 0.5 * p1_only_ll, tolerance = 1e-6)
})
