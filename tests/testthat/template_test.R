# Template for a new Stan test harness in R
# Copy, rename, and edit as needed

library(cmdstanr)
library(testthat)
library(here)

# Set up Stan model and data
stan_file <- here("tests/testthat/stan/template_test.stan")

data <- list(
  # Fill in with your test data
  N = 3,
  x = c(1.0, 2.0, 3.0)
)

source(here("tests/testthat/helper-stan.R"))

# Run Stan model using the helper
fit <- test_stan_function(
  stan_file = stan_file,
  data = data
)

gq <- fit$draws(format = "list")

# Example test: replace with your own checks
test_that("template test runs and outputs expected values", {
  expect_equal(length(gq$y), data$N)
  # Add more checks here
})
