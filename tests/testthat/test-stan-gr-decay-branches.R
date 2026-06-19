library(testthat)
library(here)
library(stringr)

test_that("cumsum (phi-difference) form equals direct form for constant rate", {
  for (k in c(1e-10, 0.01, 0.05, 0.2)) {
    data <- list(n_t = 40L, growth_rate = 0.08, kappa = k, init_g = log(3))
    fit <- test_stan_function(
      here::here("tests/testthat/stan/test_gr_decay_warp_branches.stan"), data
    )
    draws <- posterior::as_draws_df(fit$draws())
    direct <- vapply(1:40, function(i) get_stan_val(draws, "direct", i), numeric(1))
    cform  <- vapply(1:40, function(i) get_stan_val(draws, "cumsum_form", i), numeric(1))
    expect_equal(direct, cform, tolerance = 1e-9,
                 info = paste("kappa =", k))
  }
})

test_that("kappa -> 0 recovers the linear growth arm exactly", {
  data <- list(n_t = 20L, growth_rate = 0.08, kappa = 1e-12, init_g = log(3))
  fit <- test_stan_function(
    here::here("tests/testthat/stan/test_gr_decay_warp_branches.stan"), data
  )
  draws <- posterior::as_draws_df(fit$draws())
  direct <- vapply(1:20, function(i) get_stan_val(draws, "direct", i), numeric(1))
  expected <- data$init_g + data$growth_rate * (0:19)
  expect_equal(direct, expected, tolerance = 1e-6)
})
