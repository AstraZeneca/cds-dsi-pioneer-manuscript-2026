library(testthat)
library(here)
library(stringr)

test_that("pop-only config: kappa constant across patients", {
  np <- 6L
  data <- list(
    n_forecast_patients = np, n_covar = 2L,
    enable_gr_decay = 1L, enable_pop_cov_gr_decay = 0L,
    forecast_patient_idx = 1:np,
    Q_covar_design_matrix = matrix(rnorm(np * 2), np, 2),
    gr_decay_log_loc_pop = array(log(0.02)),
    gr_decay_coef_qr_pop = numeric(0)
  )
  fit <- test_stan_function(
    here::here("tests/testthat/stan/test_gr_decay_reduction.stan"), data
  )
  draws <- posterior::as_draws_df(fit$draws())
  vals <- vapply(1:np, function(i) get_stan_val(draws, "gr_decay_log_loc_patient", i), numeric(1))
  # Constant across patients (pop-only) and equal to the intercept. 1e-6 absorbs the
  # CSV round-trip precision loss on the passed-in real (sibling tests use 1e-6 too).
  expect_equal(vals, rep(vals[1], np), tolerance = 1e-12)       # strictly constant across patients
  expect_equal(vals, rep(log(0.02), np), tolerance = 1e-6)      # equals the intercept
})

test_that("flag off: kappa log-loc is zero (value unused downstream)", {
  np <- 4L
  data <- list(
    n_forecast_patients = np, n_covar = 0L,
    enable_gr_decay = 0L, enable_pop_cov_gr_decay = 0L,
    forecast_patient_idx = 1:np,
    Q_covar_design_matrix = matrix(0, np, 0),
    gr_decay_log_loc_pop = numeric(0),
    gr_decay_coef_qr_pop = numeric(0)
  )
  fit <- test_stan_function(
    here::here("tests/testthat/stan/test_gr_decay_reduction.stan"), data
  )
  draws <- posterior::as_draws_df(fit$draws())
  vals <- vapply(1:np, function(i) get_stan_val(draws, "gr_decay_log_loc_patient", i), numeric(1))
  expect_equal(vals, rep(0, np), tolerance = 1e-12)
})
