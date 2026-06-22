library(testthat)
library(here)
library(stringr)  # str_c() used inside get_stan_val() helper

phi_ref <- function(t, k) if (abs(k * t) < 1e-12) t else (1 - exp(-k * t)) / k

test_that("growth_warp matches phi and row_vector overload agrees with scalar", {
  data <- list(n_t = 5L, t = c(1, 4, 12, 52, 200), kappa = 0.05)
  fit <- test_stan_function(
    here::here("tests/testthat/stan/test_growth_warp.stan"), data
  )
  draws <- posterior::as_draws_df(fit$draws())
  phi <- vapply(1:5, function(k) get_stan_val(draws, "phi", k), numeric(1))
  phi_rv <- vapply(1:5, function(k) get_stan_val(draws, "phi_rv", k), numeric(1))
  expect_equal(phi, vapply(data$t, phi_ref, numeric(1), k = data$kappa), tolerance = 1e-6)
  expect_equal(phi, phi_rv, tolerance = 1e-10)
})

test_that("growth_warp recovers linear growth as kappa -> 0", {
  data <- list(n_t = 3L, t = c(1, 12, 52), kappa = 1e-10)
  fit <- test_stan_function(
    here::here("tests/testthat/stan/test_growth_warp.stan"), data
  )
  draws <- posterior::as_draws_df(fit$draws())
  phi <- vapply(1:3, function(k) get_stan_val(draws, "phi", k), numeric(1))
  expect_equal(phi, data$t, tolerance = 1e-6)
})

test_that("growth_warp plateaus at 1/kappa as t grows", {
  k <- 0.1
  data <- list(n_t = 1L, t = array(5000), kappa = k)
  fit <- test_stan_function(
    here::here("tests/testthat/stan/test_growth_warp.stan"), data
  )
  draws <- posterior::as_draws_df(fit$draws())
  expect_equal(get_stan_val(draws, "phi", 1), 1 / k, tolerance = 1e-4)
})

test_that("phi-differences telescope to phi(t_n)", {
  data <- list(n_t = 6L, t = c(2, 5, 9, 20, 60, 150), kappa = 0.03)
  fit <- test_stan_function(
    here::here("tests/testthat/stan/test_growth_warp.stan"), data
  )
  draws <- posterior::as_draws_df(fit$draws())
  expect_equal(
    get_stan_val(draws, "telescoped"),
    get_stan_val(draws, "phi_last"),
    tolerance = 1e-10
  )
})
