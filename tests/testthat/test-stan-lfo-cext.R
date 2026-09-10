library(testthat)
library(cmdstanr)
source(here::here("tests/testthat/helper-stan.R"))

# Regression tests for LFO correctness bugs fixed in b0a30e99.
# All tests use a minimal 2-patient fixture (1 eval-trial + 1 historical).

test_that("C-EXT: historical patient's cutoff_last_visit_idx is not truncated at cutoff", {
  # Fixture: eval_trial=1, patient 1 is trial 1, patient 2 is historical (trial 2).
  # Cutoff calendar_day=107. Without the C-EXT override, patient 2 would be
  # truncated at study-day 107 (which exceeds all their visits anyway for this
  # setup, but the important property is that the override sets the idx to the
  # *last* visit in the patient block unconditionally).
  fit <- test_stan_function(
    "tests/testthat/stan/test_lfo_cext_all.stan",
    data = list(lfo_eval_trial = 1L)
  )
  draws_df <- posterior::as_draws_df(fit$draws())

  # Patient 2 (historical): visits are at global indices 6-10. The C-EXT override
  # must set cutoff_last_visit_idx[2] = 10 (full visit end), not a truncated index.
  expect_equal(
    get_stan_val(draws_df, "cext_historical_cutoff_idx"),
    10,
    label = "historical patient cutoff_last_visit_idx should be its last visit index"
  )

  # Patient 1 (eval-trial): visits at global indices 1-5. Cutoff calendar_day=107
  # → study-day cutoff = 107-100+1 = 8 → last visit on or before week 8 is visit
  # at study-day 8 (idx=4 within patient block, global idx=4).
  expect_equal(
    get_stan_val(draws_df, "cext_eval_trial_cutoff_idx"),
    4,
    label = "eval-trial patient cutoff_last_visit_idx should be truncated at cutoff"
  )
})

test_that("B: lfo_testing_patient_idx is 0 for historical patients, non-zero for eval-trial", {
  fit <- test_stan_function(
    "tests/testthat/stan/test_lfo_cext_all.stan",
    data = list(lfo_eval_trial = 1L)
  )
  draws_df <- posterior::as_draws_df(fit$draws())

  # Eval-trial patient 1 should be at slot 1 in the OOS output vectors
  expect_equal(
    get_stan_val(draws_df, "b_eval_trial_testing_idx"),
    1,
    label = "eval-trial patient should have testing index 1"
  )

  # Historical patient 2 is not in the eval trial; reverse-lookup must return 0
  expect_equal(
    get_stan_val(draws_df, "b_historical_testing_idx"),
    0,
    label = "historical patient should have testing index 0 (excluded from OOS scoring)"
  )
})

test_that("H3: GQ array cell [1,2] is zero (not NaN) at last cutoff after zero-init", {
  # At the global last cutoff, only [1,1] is written; [1,2] must be 0 from pre-init,
  # not NaN from uninitialized Stan memory.
  fit <- test_stan_function(
    "tests/testthat/stan/test_lfo_cext_all.stan",
    data = list(lfo_eval_trial = 1L)
  )
  draws_df <- posterior::as_draws_df(fit$draws())

  expect_equal(
    get_stan_val(draws_df, "h3_cell_written"),
    -1.5,
    label = "written cell [1,1] should hold the assigned value"
  )
  expect_equal(
    get_stan_val(draws_df, "h3_cell_unwritten"),
    0.0,
    label = "unwritten cell [1,2] must be 0.0 after zero-init, not NaN"
  )
})
