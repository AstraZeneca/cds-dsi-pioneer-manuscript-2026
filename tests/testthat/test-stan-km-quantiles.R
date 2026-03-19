# Tests for estimate_kaplan_meier, km_quantiles, and calc_km_pfs_n (Issue 79)
#
# Exercises all three KM utility functions in pfs.stanfunctions:
#   - estimate_kaplan_meier()  (4 scenarios: mixed, all-events, all-censored, tied)
#   - km_quantiles()           (median, Q1/Q3, cannot-calculate, multi-quantile)
#   - calc_km_pfs_n()          (exact integer times, fractional interpolation)
#
# All expected values are derived by hand from the KM step-function formula.
# The Stan model is compiled with fixed_param = TRUE so no sampling occurs.
# The single output variable n_failures counts scenario mismatches.

library(testthat)
library(cmdstanr)
library(here)

test_that("km_and_quantiles: all scenarios produce expected outputs (n_failures == 0)", {
  stan_file <- here(
    "tests", "testthat", "stan",
    "test_km_and_quantiles_all.stan"
  )

  fit <- test_stan_function(
    stan_file = stan_file,
    data      = list(dummy = 0L),
    seed      = 42L
  )

  draws_df   <- posterior::as_draws_df(fit$draws())
  n_failures <- as.integer(draws_df$n_failures[[1]])

  expect_equal(
    n_failures,
    0L,
    label = "n_failures",
    info  = paste0(
      "One or more KM/quantile scenarios produced unexpected outputs. ",
      "See Stan print() messages above for which scenarios failed."
    )
  )
})
