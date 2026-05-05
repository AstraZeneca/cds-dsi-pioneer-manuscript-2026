library(testthat)

# helper-util.R is auto-loaded by testthat

test_that("r_count_positive counts elements > 0 correctly", {
  expect_equal(r_count_positive(c(3, -1, 0, 5, 2)), 3L)
  expect_equal(r_count_positive(c(-1, -2, -3)), 0L)
  expect_equal(r_count_positive(c(0, 0, 0)), 0L)
  expect_equal(r_count_positive(c(1, 2, 3)), 3L)
  expect_equal(r_count_positive(c(1L)), 1L)
})

test_that("r_which returns correct indices for mask > 0", {
  expect_equal(r_which(c(1, 0, 1, 0, 1)), c(1L, 3L, 5L))
  expect_equal(r_which(c(0, 0, 0)), integer(0))
  expect_equal(r_which(c(1, 1, 1)), c(1L, 2L, 3L))
})

test_that("r_which inverse returns indices where mask <= 0", {
  expect_equal(r_which(c(1, 0, 1, 0, 1), inverse = TRUE), c(2L, 4L))
  expect_equal(r_which(c(0, 0, 0), inverse = TRUE), c(1L, 2L, 3L))
})

test_that("r_rep_each repeats each element correctly", {
  expect_equal(r_rep_each(c(1L, 2L), 3L), c(1L, 1L, 1L, 2L, 2L, 2L))
  expect_equal(r_rep_each(c(10L, 20L), 2L), c(10L, 10L, 20L, 20L))
  expect_equal(r_rep_each(c(5L), 4L), c(5L, 5L, 5L, 5L))
})

test_that("r_months_to_weeks converts 12 months to ~52.18 weeks", {
  expect_equal(r_months_to_weeks(12L), 365.25 / 7, tolerance = 1e-6)
  expect_equal(r_months_to_weeks(0L), 0)
})

test_that("r_calendar_date_to_study_date_scalar computes study date correctly", {
  expect_equal(r_calendar_date_to_study_date_scalar(100L, 110L), 11L)
  expect_equal(r_calendar_date_to_study_date_scalar(1L, 1L), 1L)
  expect_equal(r_calendar_date_to_study_date_scalar(50L, 50L), 1L)
})

test_that("r_calendar_date_to_study_date_vec vectorises over first_vec", {
  first_vec <- c(100L, 105L, 108L)
  expect_equal(
    r_calendar_date_to_study_date_vec(first_vec, 110L),
    c(11L, 6L, 3L)
  )
})

test_that("r_id2idx maps ids to 1-based contiguous indices", {
  expect_equal(r_id2idx(c(3L, 1L, 2L)), c(3L, 1L, 2L))
  expect_equal(r_id2idx(c(5L, 3L, 4L)), c(3L, 1L, 2L))
  expect_equal(r_id2idx(c(1L, 1L, 1L)), c(1L, 1L, 1L))
})
