library(testthat)

test_that("r_find_first: basic find_first semantics", {
  # found at position 3
  expect_equal(r_find_first(c(3L, 3L, 4L), c(4L), 1L), 3L)
  # found at position 1
  expect_equal(r_find_first(c(4L, 3L, 3L), c(4L), 1L), 1L)
  # not found
  expect_equal(r_find_first(c(3L, 3L, 3L), c(4L), 1L), 0L)
  # run of 2 required
  expect_equal(r_find_first(c(3L, 4L, 4L, 3L), c(4L), 2L), 2L)
  # run of 2, not found (only one 4)
  expect_equal(r_find_first(c(3L, 4L, 3L, 4L), c(4L), 2L), 0L)
  # multiple values in 'what'
  expect_equal(r_find_first(c(3L, 1L, 3L), c(1L, 2L), 1L), 2L)
})

test_that("r_map_idx_to_week: maps obs vs forecast correctly", {
  curr   <- c(4L, 8L, 12L)
  fore   <- c(16L, 20L)
  max_t  <- 52L
  # idx in obs range
  expect_equal(r_map_idx_to_week(1L, curr, fore, max_t), 4L)
  expect_equal(r_map_idx_to_week(3L, curr, fore, max_t), 12L)
  # idx in forecast range
  expect_equal(r_map_idx_to_week(4L, curr, fore, max_t), 16L)
  expect_equal(r_map_idx_to_week(5L, curr, fore, max_t), 20L)
  # idx == 0 → censoring sentinel
  expect_equal(r_map_idx_to_week(0L, curr, fore, max_t), 52L)
  # out of forecast range → censoring sentinel
  expect_equal(r_map_idx_to_week(6L, curr, fore, max_t), 52L)
})

test_that("r_find_first_week: week and right_censored", {
  arr  <- c(3L, 3L, 4L)
  curr <- c(4L, 8L, 12L)
  fore <- c(16L, 20L)
  mt   <- 52L

  # found in obs → week = curr_visits[3] = 12, rc = 0
  res <- r_find_first_week(arr, c(4L), 1L, curr, fore, mt)
  expect_equal(res$week, 12L)
  expect_equal(res$right_censored, 0L)

  # not found → week = max_all_t, rc = 1
  res <- r_find_first_week(c(3L, 3L, 3L), c(4L), 1L, curr, fore, mt)
  expect_equal(res$week, 52L)
  expect_equal(res$right_censored, 1L)

  # found in forecast (idx=4 > n_obs=3) → week = fore[1] = 16
  arr4 <- c(3L, 3L, 3L, 4L)
  res  <- r_find_first_week(arr4, c(4L), 1L, curr, fore, mt)
  expect_equal(res$week, 16L)
  expect_equal(res$right_censored, 0L)

  # empty curr_visits → idx maps directly to forecast
  res <- r_find_first_week(c(3L, 4L), c(4L), 1L, integer(0), c(8L, 16L), mt)
  expect_equal(res$week, 16L)   # idx=2, forecast_idx=2
  expect_equal(res$right_censored, 0L)
})

test_that("r_find_first_forecast_week: delegates with empty curr_visits", {
  # found at forecast position 1
  res <- r_find_first_forecast_week(c(4L), c(4L), 1L, c(8L), 52L)
  expect_equal(res$week, 8L)
  expect_equal(res$right_censored, 0L)

  # not found
  res <- r_find_first_forecast_week(c(3L, 3L), c(4L), 1L, c(8L, 16L), 52L)
  expect_equal(res$week, 52L)
  expect_equal(res$right_censored, 1L)
})
