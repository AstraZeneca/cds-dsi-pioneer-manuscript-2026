# Tests for compute_trial_cif() in pfs.stanfunctions (Issue 79)
#
# Exercises all four patient classification branches:
#   0→2 direct death  : rc=0 AND os_censored=0 AND pfs==os
#   0→1 progression   : rc=0 AND NOT 0→2
#   0→3 dropout       : rc=1 AND pfs <= max_t
#   fully censored    : rc=1 AND pfs > max_t  (no CIF contribution)
#
# Three 4-patient scenarios:
#   CIF-A: one of each cause + fully censored → total final = 3/4
#   CIF-B: two 0→1 at different times → cumulative 0→1 reaches 2/4 → total = 1
#   CIF-C: all fully censored → all CIFs stay zero throughout
#
# The Stan model is compiled with fixed_param = TRUE so no sampling occurs.

library(testthat)
library(cmdstanr)
library(here)

test_that("compute_trial_cif: all scenarios produce expected CIF values (n_failures == 0)", {
  stan_file <- here(
    "tests", "testthat", "stan",
    "test_cif_computation_all.stan"
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
    info  = str_c(
      "One or more CIF scenarios produced unexpected outputs. ",
      "See Stan print() messages above for which scenarios failed."
    )
  )
})
