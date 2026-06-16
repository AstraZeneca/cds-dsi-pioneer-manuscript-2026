library(testthat)
library(cmdstanr)
library(here)

test_that("standardize_velocity and central_difference_row produce expected values (n_failures == 0)", {
  stan_file <- here("tests", "testthat", "stan", "test_velocity_covar_all.stan")

  fit <- test_stan_function(
    stan_file = stan_file,
    data      = list(dummy = 0L),
    seed      = 42L
  )

  draws_df   <- posterior::as_draws_df(fit$draws())
  n_failures <- as.integer(draws_df$n_failures[[1]])

  expect_equal(n_failures, 0L, label = "n_failures")
})
