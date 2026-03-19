# Tests for classify-first endpoint functions (Issue 79)
#
# Exercises all four pure functions in pfs.stanfunctions:
#   - classify_spop_exit()    (9 scenarios)
#   - derive_spop_pfs()       (7 scenarios)
#   - classify_sample_exit()  (7 scenarios)
#   - derive_sample_pfs()     (7 scenarios)
#
# And the deterministic (non-RNG) branches of the two OS derivation functions
# in multistate.stanfunctions:
#   - derive_spop_os_rng()    (4 deterministic scenarios)
#   - derive_sample_os_rng()  (6 deterministic scenarios)
#
# The Stan model is compiled with fixed_param = TRUE so no sampling occurs.
# All expected outputs are computed by hand from the function logic.
# The single output variable n_failures counts scenario mismatches.

library(testthat)
library(cmdstanr)
library(here)

test_that("classify_first_exit: all 40 scenarios produce expected outputs (n_failures == 0)", {
  stan_file <- here(
    "tests", "testthat", "stan",
    "test_classify_first_exit_all.stan"
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
      "One or more classify-first scenarios produced unexpected outputs. ",
      "See Stan print() messages above for which scenarios failed."
    )
  )
})
