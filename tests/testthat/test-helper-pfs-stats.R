library(testthat)

# ---------------------------------------------------------------------------
# r_calculate_log_marginal_exit_prob
# ---------------------------------------------------------------------------

test_that("r_calculate_log_marginal_exit_prob: scalar and multi-step values", {
  lcs <- c(-0.1, -0.2, -0.3)

  result <- r_calculate_log_marginal_exit_prob(lcs)

  # t = 1: log1p(-exp(-0.1))
  expect_equal(result[1], log1p(-exp(-0.1)), tolerance = 1e-12,
    label = "lmep t=1")

  # t = 2: log1p(-exp(-0.2)) + (-0.1)
  expect_equal(result[2], log1p(-exp(-0.2)) + (-0.1), tolerance = 1e-12,
    label = "lmep t=2")

  # t = 3: log1p(-exp(-0.3)) + (-0.1) + (-0.2)
  expect_equal(result[3], log1p(-exp(-0.3)) + (-0.1) + (-0.2), tolerance = 1e-12,
    label = "lmep t=3")

  # output length matches input length
  expect_length(result, 3L)
})

test_that("r_calculate_log_marginal_exit_prob: single-element input", {
  result <- r_calculate_log_marginal_exit_prob(c(-0.5))
  expect_equal(result[1], log1p(-exp(-0.5)), tolerance = 1e-12)
  expect_length(result, 1L)
})

# ---------------------------------------------------------------------------
# r_survival_quantiles
# ---------------------------------------------------------------------------

test_that("r_survival_quantiles: median of 5 evenly-spaced patients is calculable", {
  surv  <- c(2L, 4L, 6L, 8L, 10L)
  last  <- 12L
  res   <- r_survival_quantiles(surv, last, c(0.5))

  expect_equal(res$cannot_calculate[1], 0L,
    label = "median calculable when all times <= last_surv_time")
  # N=5, p=0.5: k loops while 0.5 >= k/5 -> k ends at 3; pos = 0.5*4+1 = 3; d = 3-2 = 1
  # quantile = sorted[2] + 1*(sorted[3]-sorted[2]) = 4 + 1*(6-4) = 6
  expect_equal(res$quantiles[1], 6.0, tolerance = 1e-10,
    label = "median value for c(2,4,6,8,10)")
})

test_that("r_survival_quantiles: lower quartile of 5 patients", {
  surv <- c(2L, 4L, 6L, 8L, 10L)
  last <- 12L
  res  <- r_survival_quantiles(surv, last, c(0.25))
  expect_equal(res$cannot_calculate[1], 0L)
  # p=0.25, N=5: k loops while 0.25 >= k/5 -> k ends at 2
  # pos = 0.25*4+1 = 2; d = 2-1 = 1; quantile = sorted[1]+1*(sorted[2]-sorted[1]) = 2+1*(4-2) = 4
  expect_equal(res$quantiles[1], 4.0, tolerance = 1e-10,
    label = "lower quartile value")
})

test_that("r_survival_quantiles: cannot_calculate=1 when all times beyond last_surv_time", {
  surv <- c(5L, 6L, 7L, 8L, 9L)
  last <- 3L
  res  <- r_survival_quantiles(surv, last, c(0.5))
  expect_equal(res$cannot_calculate[1], 1L,
    label = "cannot_calculate when sorted[k-1] > last_surv_time")
  expect_equal(res$quantiles[1], 0.0,
    label = "quantile set to 0.0 when not calculable")
})

test_that("r_survival_quantiles: multiple quantiles returned correctly", {
  surv <- c(2L, 4L, 6L, 8L, 10L)
  last <- 12L
  res  <- r_survival_quantiles(surv, last, c(0.25, 0.5))
  expect_length(res$quantiles, 2L)
  expect_length(res$cannot_calculate, 2L)
  expect_equal(res$cannot_calculate, c(0L, 0L))
})

# ---------------------------------------------------------------------------
# r_km_quantiles
# ---------------------------------------------------------------------------

test_that("r_km_quantiles: median from decreasing KM curve", {
  # km_survival: S(0..5), length 6, decreasing
  km  <- c(1.0, 0.8, 0.6, 0.4, 0.2, 0.1)
  res <- r_km_quantiles(km, c(0.5))

  expect_equal(res$cannot_calculate[1], 0L,
    label = "median calculable from km curve")
  # S drops from 0.6 to 0.4 between t=2 (index 3) and t=3 (index 4)
  # weight = (0.6 - 0.5) / (0.6 - 0.4) = 0.5; quantile = (4-2) + 0.5 = 2.5
  expect_equal(res$quantiles[1], 2.5, tolerance = 1e-10,
    label = "km median interpolated between t=2 and t=3")
})

test_that("r_km_quantiles: p > km_survival[1] → quantile=0, cannot_calc=0", {
  km  <- c(0.8, 0.6, 0.4, 0.2, 0.1)
  res <- r_km_quantiles(km, c(0.9))
  expect_equal(res$cannot_calculate[1], 0L,
    label = "p > S(0): cannot_calc = 0 (time = 0)")
  expect_equal(res$quantiles[1], 0.0,
    label = "p > S(0): quantile = 0")
})

test_that("r_km_quantiles: p never crossed → cannot_calc=1", {
  # KM stays above 0.1 throughout
  km  <- c(1.0, 0.8, 0.6, 0.4, 0.2)
  res <- r_km_quantiles(km, c(0.05))
  expect_equal(res$cannot_calculate[1], 1L,
    label = "p=0.05 never reached: cannot_calc stays 1")
})

test_that("r_km_quantiles: both quartiles returned and calculable", {
  km  <- c(1.0, 0.8, 0.6, 0.4, 0.2, 0.1)
  res <- r_km_quantiles(km, c(0.25, 0.5))
  expect_equal(res$cannot_calculate, c(0L, 0L))
  expect_length(res$quantiles, 2L)
  # upper quartile (p=0.25): S drops from 0.4 to 0.2 between t=3 (idx 4) and t=4 (idx 5)
  # weight = (0.4 - 0.25)/(0.4 - 0.2) = 0.75; quantile = (5-2) + 0.75 = 3.75
  expect_equal(res$quantiles[1], 3.75, tolerance = 1e-10,
    label = "upper quartile (p=0.25)")
})

# ---------------------------------------------------------------------------
# r_calc_km_pfs_n
# ---------------------------------------------------------------------------

test_that("r_calc_km_pfs_n: boundary and interpolation cases", {
  km <- c(1.0, 0.8, 0.6, 0.4, 0.2, 0.1)   # indices 1..6 = t=0..5

  # n = 0 → km[1]
  expect_equal(r_calc_km_pfs_n(km, 0), 1.0,
    label = "n=0 returns km[1]")

  # n >= max_t (5) → km[6]
  expect_equal(r_calc_km_pfs_n(km, 5), 0.1,
    label = "n=max_t returns km[end]")
  expect_equal(r_calc_km_pfs_n(km, 10), 0.1,
    label = "n > max_t returns km[end]")

  # fractional: n=2.5 → (1-0.5)*km[3] + 0.5*km[4] = 0.5*0.6 + 0.5*0.4 = 0.5
  expect_equal(r_calc_km_pfs_n(km, 2.5), 0.5, tolerance = 1e-10,
    label = "n=2.5 interpolates between km[3] and km[4]")

  # integer: n=2 → km[3] = 0.6 (no interpolation)
  expect_equal(r_calc_km_pfs_n(km, 2), 0.6, tolerance = 1e-10,
    label = "n=2 returns km[3]")
})
