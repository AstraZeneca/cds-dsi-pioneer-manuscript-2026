library(testthat)
library(here)

# Regression guard for BUG 2 (forecast warp re-anchored at the cutoff instead of the
# baseline week). The forecast must be the exact analytic continuation of the in-sample
# arm on ONE baseline-anchored clock. With baseline < cutoff < horizon (all distinct),
# this test can tell a baseline-anchored forecast from a cutoff-anchored one — the
# discrimination test_gr_decay_plateau.stan cannot make (it pins times[1] = baseline = 0).

stan_file <- here::here("tests/testthat/stan/test_gr_decay_continuity.stan")

base_data <- list(
  baseline_week = 0,
  cutoff_week = 24,        # last observed visit, well after baseline
  n_forecast = 50L,
  horizon_week = 520,      # ~10 yr: long lever so re-acceleration shows clearly
  init_decrease = log(5),
  init_growth = log(2),
  decrease_rate = 0.03,
  growth_rate = 0.08,
  kappa = 0.05
)

test_that("baseline-anchored forecast is the exact analytic continuation at the seam and beyond", {
  fit <- test_stan_function(stan_file, base_data)
  draws <- posterior::as_draws_df(fit$draws())

  seam_correct  <- get_stan_val(draws, "seam_growth_correct")
  cutoff_state  <- get_stan_val(draws, "cutoff_state_growth")
  max_abs_err   <- get_stan_val(draws, "max_abs_err_correct")

  # (a) Seam continuity: first forecast point (anchor) == in-sample cutoff state.
  expect_equal(seam_correct, cutoff_state, tolerance = 1e-9)

  # (b) The whole correct forecast growth arm matches init + growth_rate*phi(week - baseline).
  expect_lt(max_abs_err, 1e-6)
})

test_that("the cutoff-anchored (BUG 2) forecast is strictly higher — the test discriminates the bug", {
  fit <- test_stan_function(stan_file, base_data)
  draws <- posterior::as_draws_df(fit$draws())

  final_correct  <- get_stan_val(draws, "final_growth_correct")
  final_buggy    <- get_stan_val(draws, "final_growth_buggy")
  final_analytic <- get_stan_val(draws, "final_growth_analytic")
  buggy_gap      <- get_stan_val(draws, "buggy_minus_correct_final")

  # Correct forecast tracks the analytic continuation; buggy one does not.
  expect_equal(final_correct, final_analytic, tolerance = 1e-6)

  # Re-anchoring at the cutoff resets the warp clock (phi'(0)=1), so the growth rate
  # re-accelerates => buggy final growth arm is STRICTLY HIGHER than correct. This is the
  # assertion that would FAIL if BUG 2 were reintroduced in the production tv_factor.
  expect_gt(buggy_gap, 1e-3)
  expect_gt(final_buggy, final_correct)
})
