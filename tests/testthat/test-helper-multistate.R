library(testthat)

test_that("r_calc_ms_stl: censored patient accumulates survival terms only", {
  lcs    <- matrix(-0.1, nrow = 1, ncol = 8)
  result <- r_calc_ms_stl(6L, 1L, lcs)
  # censored at week 6: last_surv = 6, survival = 6 * (-0.1)
  expect_equal(result[1], 6 * (-0.1), tolerance = 1e-9)
})

test_that("r_calc_ms_stl: event patient accumulates survival[1..t-1] + hazard[t]", {
  lcs    <- matrix(-0.1, nrow = 1, ncol = 8)
  result <- r_calc_ms_stl(4L, 0L, lcs)
  # event at week 4: last_surv = 3, survival = 3 * (-0.1), hazard at col 4
  expected <- 3 * (-0.1) + log1p(-exp(-0.1))
  expect_equal(result[1], expected, tolerance = 1e-9)
})

test_that("r_calc_ms_stl: event at week 1 — no survival accumulation, only hazard", {
  lcs    <- matrix(-0.2, nrow = 1, ncol = 5)
  result <- r_calc_ms_stl(1L, 0L, lcs)
  # last_surv = 0, no survival term; hazard at col 1
  expected <- log1p(-exp(-0.2))
  expect_equal(result[1], expected, tolerance = 1e-9)
})

test_that("r_calc_ms_stl: multi-patient vectorized output", {
  n     <- 3L; max_t <- 8L
  lcs   <- matrix(-0.08, nrow = n, ncol = max_t)
  event <- c(4L, 6L, 1L)
  cens  <- c(0L, 1L, 0L)

  result <- r_calc_ms_stl(event, cens, lcs)

  expect_length(result, 3L)
  # Patient 1: event at 4, last_surv=3
  expect_equal(result[1], 3 * (-0.08) + log1p(-exp(-0.08)), tolerance = 1e-9)
  # Patient 2: censored at 6
  expect_equal(result[2], 6 * (-0.08), tolerance = 1e-9)
  # Patient 3: event at 1, last_surv=0
  expect_equal(result[3], log1p(-exp(-0.08)), tolerance = 1e-9)
})

test_that("r_multistate_lpmf_01only equals sum(r_calc_ms_stl)", {
  n     <- 4L; max_t <- 6L
  lcs   <- matrix(-0.1, nrow = n, ncol = max_t)
  event <- c(2L, 4L, 6L, 3L)
  cens  <- c(0L, 1L, 1L, 0L)

  expected_sum <- sum(r_calc_ms_stl(event, cens, lcs))
  result       <- r_multistate_lpmf_01only(event, cens, lcs)
  expect_equal(result, expected_sum, tolerance = 1e-9)
})

test_that("r_multistate_state0_cens accumulates both transition survival terms", {
  t_cens <- 4L
  lcs_01 <- rep(-0.1, 8)
  lcs_02 <- rep(-0.05, 8)
  result <- r_multistate_state0_cens(t_cens, lcs_01, lcs_02)
  # sum(-0.1, ...) + sum(-0.05, ...) over 4 weeks
  expected <- 4 * (-0.1) + 4 * (-0.05)
  expect_equal(result, expected, tolerance = 1e-9)
})

test_that("r_multistate_state0_cens: t_cens=0 gives lp=0", {
  result <- r_multistate_state0_cens(0L, rep(-0.1, 5), rep(-0.05, 5))
  expect_equal(result, 0.0, tolerance = 1e-12)
})
