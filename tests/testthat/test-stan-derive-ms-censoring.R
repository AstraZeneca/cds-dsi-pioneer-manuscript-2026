library(testthat)
library(here)
library(stringr)

# ============================================================================
# Tests for derive_ms_censoring_indicators()
#
# Covers all six patient patterns:
#   admin_censored   state=0 time_01=0 time_32=0  → c02=1 c12=1 c32=1
#   true_dropout     state=0 time_01=0 time_32=0  → c02=1 c12=1 c32=1
#   progressed_alive state=1 time_01>0 time_32=0  → c02=1 c12=1 c32=1
#   progressed_died  state=2 time_01>0 time_32=0  → c02=1 c12=0 c32=1
#   died_on_trial    state=2 time_01=0 time_32=0  → c02=0 c12=1 c32=1
#   died_off_trial   state=3 time_01=0 time_32>0  → c02=1 c12=1 c32=0
# ============================================================================

run_censoring_indicators <- function(final_state, time_01, time_32) {
  test_stan_function(
    here("tests", "testthat", "stan", "test_derive_ms_censoring_all.stan"),
    list(
      N              = length(final_state),
      ms_final_state = final_state,
      ms_time_01     = time_01,
      ms_time_32     = time_32
    )
  )
}

# Helper: extract censored_XY[i] from draws
get_cens <- function(d, xy, i) get_stan_val(d, paste0("censored_", xy), i)

test_that("derive_ms_censoring_indicators: all six patient patterns correct", {
  # One patient per pattern (order matches table in function doc)
  final_state <- c(
    0L,  # admin_censored
    0L,  # true_dropout      (same as admin_censored at data level)
    1L,  # progressed_alive
    2L,  # progressed_died   (0→1→2: time_01 > 0)
    2L,  # died_on_trial     (0→2:   time_01 == 0)
    3L   # died_off_trial
  )
  time_01 <- c(0L, 0L, 5L, 5L, 0L, 0L)
  time_32 <- c(0L, 0L, 0L, 0L, 0L, 8L)

  fit <- run_censoring_indicators(final_state, time_01, time_32)
  d   <- posterior::as_draws_df(fit$draws())

  # Expected: (censored_02, censored_12, censored_32) per patient
  expected <- list(
    admin_censored   = c(1L, 1L, 1L),
    true_dropout     = c(1L, 1L, 1L),
    progressed_alive = c(1L, 1L, 1L),
    progressed_died  = c(1L, 0L, 1L),
    died_on_trial    = c(0L, 1L, 1L),
    died_off_trial   = c(1L, 1L, 0L)
  )

  patterns <- names(expected)
  for (i in seq_along(patterns)) {
    pat  <- patterns[i]
    exp  <- expected[[pat]]
    for (j in seq_along(c("02", "12", "32"))) {
      xy <- c("02", "12", "32")[j]
      expect_equal(
        get_cens(d, xy, i), exp[j],
        label = sprintf("%s: censored_%s", pat, xy)
      )
    }
  }
})

test_that("derive_ms_censoring_indicators: censored_02 is 0 only for direct-death patients", {
  # Only final_state==2 with time_01==0 should have censored_02==0
  final_state <- c(2L, 2L, 1L, 0L)
  time_01     <- c(0L, 5L, 5L, 0L)  # first is direct death, rest are not
  time_32     <- c(0L, 0L, 0L, 0L)

  fit <- run_censoring_indicators(final_state, time_01, time_32)
  d   <- posterior::as_draws_df(fit$draws())

  expect_equal(get_cens(d, "02", 1L), 0L, label = "state2+t01=0 → event")
  expect_equal(get_cens(d, "02", 2L), 1L, label = "state2+t01>0 → censored")
  expect_equal(get_cens(d, "02", 3L), 1L, label = "state1       → censored")
  expect_equal(get_cens(d, "02", 4L), 1L, label = "state0       → censored")
})

test_that("derive_ms_censoring_indicators: censored_12 is 0 only for post-progression death", {
  final_state <- c(2L, 2L, 1L, 0L)
  time_01     <- c(5L, 0L, 5L, 0L)  # first is post-progression death
  time_32     <- c(0L, 0L, 0L, 0L)

  fit <- run_censoring_indicators(final_state, time_01, time_32)
  d   <- posterior::as_draws_df(fit$draws())

  expect_equal(get_cens(d, "12", 1L), 0L, label = "state2+t01>0 → event")
  expect_equal(get_cens(d, "12", 2L), 1L, label = "state2+t01=0 → censored")
  expect_equal(get_cens(d, "12", 3L), 1L, label = "state1       → censored")
  expect_equal(get_cens(d, "12", 4L), 1L, label = "state0       → censored")
})

test_that("derive_ms_censoring_indicators: censored_32 is 0 only when off-trial death observed", {
  final_state <- c(3L, 3L, 0L, 2L)
  time_01     <- c(0L, 0L, 0L, 0L)
  time_32     <- c(8L, 0L, 0L, 0L)  # only first has observed off-trial death

  fit <- run_censoring_indicators(final_state, time_01, time_32)
  d   <- posterior::as_draws_df(fit$draws())

  expect_equal(get_cens(d, "32", 1L), 0L, label = "time_32>0 → event")
  expect_equal(get_cens(d, "32", 2L), 1L, label = "state3+time_32=0 → censored")
  expect_equal(get_cens(d, "32", 3L), 1L, label = "state0+time_32=0 → censored")
  expect_equal(get_cens(d, "32", 4L), 1L, label = "state2+time_32=0 → censored")
})

test_that("derive_ms_censoring_indicators: mutually exclusive — no patient fires both 02 and 12", {
  # A patient can never be both a direct death (02 event) and a post-progression death (12 event)
  final_state <- c(2L, 2L, 0L, 1L, 3L)
  time_01     <- c(0L, 7L, 0L, 4L, 0L)
  time_32     <- c(0L, 0L, 0L, 0L, 3L)

  fit <- run_censoring_indicators(final_state, time_01, time_32)
  d   <- posterior::as_draws_df(fit$draws())

  for (i in seq_along(final_state)) {
    c02 <- get_cens(d, "02", i)
    c12 <- get_cens(d, "12", i)
    expect_true(
      c02 == 1L || c12 == 1L,
      label = sprintf("patient %d: at most one of censored_02/12 can be 0", i)
    )
  }
})
