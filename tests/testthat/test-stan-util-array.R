library(testthat)
library(cmdstanr)
library(posterior)

# helper-stan.R and helper-util.R are auto-loaded by testthat

test_that("Stan util array functions: count_positive, which, rep_each, months_to_weeks, calendar_date_to_study_date, id2idx", {
  arr_int       <- c(3L, -1L, 0L, 5L, 2L)
  mask          <- c(1L, 0L, 1L, 0L, 1L)
  to_repeat     <- c(10L, 20L)
  id_arr        <- c(3L, 1L, 4L, 1L, 5L)
  search_arr    <- c(1L, 3L, 2L, 3L, 3L, 1L)

  stan_data <- list(
    N              = 5L,
    N_REPEAT       = 2L,
    arr_int        = arr_int,
    mask           = mask,
    to_repeat      = to_repeat,
    repeats        = 2L,
    mon            = 12L,
    first_calendar = 100L,
    calendar_date  = 110L,
    id_arr         = id_arr,
    N_ALL          = 6L,
    search_arr     = search_arr,
    search_what    = 3L,
    n_succ         = 1L
  )

  fit      <- test_stan_function(
    stan_file = here::here("tests/testthat/stan/test_util_array_all.stan"),
    data      = stan_data
  )
  draws_df <- posterior::as_draws_df(fit$draws())
  get_val  <- function(var, ...) get_stan_val(draws_df, var, ...)

  # --- count_positive ---
  expect_equal(
    get_val("cnt_pos"),
    r_count_positive(arr_int),
    label = "cnt_pos"
  )
  # r_count_positive(c(3,-1,0,5,2)) == 3
  expect_equal(get_val("cnt_pos"), 3, label = "cnt_pos == 3")

  # --- which (normal) ---
  expect_equal(get_val("which_idx", 1), 1, label = "which_idx[1]")
  expect_equal(get_val("which_idx", 2), 3, label = "which_idx[2]")
  expect_equal(get_val("which_idx", 3), 5, label = "which_idx[3]")

  # --- which (inverse) ---
  expect_equal(get_val("which_inv_idx", 1), 2, label = "which_inv_idx[1]")
  expect_equal(get_val("which_inv_idx", 2), 4, label = "which_inv_idx[2]")

  # --- rep_each ---
  # to_repeat = c(10, 20), repeats = 2 → c(10, 10, 20, 20)
  expect_equal(get_val("rep_each_out", 1), 10, label = "rep_each_out[1]")
  expect_equal(get_val("rep_each_out", 2), 10, label = "rep_each_out[2]")
  expect_equal(get_val("rep_each_out", 3), 20, label = "rep_each_out[3]")
  expect_equal(get_val("rep_each_out", 4), 20, label = "rep_each_out[4]")

  # --- months_to_weeks ---
  expect_equal(
    get_val("weeks"),
    r_months_to_weeks(12L),
    tolerance = 1e-6,
    label = "weeks ≈ r_months_to_weeks(12)"
  )

  # --- calendar_date_to_study_date (scalar) ---
  # first_calendar=100, calendar_date=110 → 110 - 100 + 1 = 11
  expect_equal(get_val("study_date"), 11, label = "study_date == 11")
  expect_equal(
    get_val("study_date"),
    r_calendar_date_to_study_date_scalar(100L, 110L),
    label = "study_date matches R oracle"
  )

  # --- calendar_date_to_study_date (vectorised over arr_int as first_dates) ---
  # arr_int = c(3,-1,0,5,2), calendar_date = 110
  # study_dates_vec[i] = 110 - arr_int[i] + 1
  expected_vec <- r_calendar_date_to_study_date_vec(arr_int, 110L)
  for (i in seq_along(arr_int)) {
    expect_equal(
      get_val("study_dates_vec", i),
      expected_vec[i],
      label = paste0("study_dates_vec[", i, "]")
    )
  }

  # --- id2idx ---
  # id_arr = c(3,1,4,1,5), min=1 → subtract 0 → c(3,1,4,1,5) - 1 + 1 = c(3,1,4,1,5)
  expected_idx <- r_id2idx(id_arr)
  expect_equal(get_val("idx_arr", 1), expected_idx[1], label = "idx_arr[1]")
  expect_equal(get_val("idx_arr", 2), expected_idx[2], label = "idx_arr[2]")
  expect_equal(get_val("idx_arr", 3), expected_idx[3], label = "idx_arr[3]")
  expect_equal(get_val("idx_arr", 4), expected_idx[4], label = "idx_arr[4]")
  expect_equal(get_val("idx_arr", 5), expected_idx[5], label = "idx_arr[5]")
})
