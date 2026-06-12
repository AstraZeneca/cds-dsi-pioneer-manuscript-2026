library(testthat)
library(here)

# Hand-computed expectations:
# state row = (log_decrease, log_growth); baseline multiplies inside log.
# 2-arg:  log(baseline) + log(exp(ld) + exp(lg))
# 3-arg:  log(baseline) + log(exp(ld) + exp(lg) + exp(static_log_level))
test_that("calc_log_burden_mean: 3-arg with -inf static equals 2-arg", {
  ld <- c(log(40), log(20))
  lg <- c(log(10), log(50))
  baseline <- 100
  data <- list(
    n_visits = 2L,
    patient_states = cbind(ld, lg),
    baseline = baseline,
    static_log_level = -Inf
  )
  fit <- test_stan_function(
    here::here("tests/testthat/stan/test_calc_log_burden_mean_all.stan"),
    data
  )
  out_2 <- as.numeric(fit$draws("out_2arg", format = "draws_matrix"))
  out_3 <- as.numeric(fit$draws("out_3arg", format = "draws_matrix"))
  expect_equal(out_3, out_2, tolerance = 1e-10)
})

test_that("calc_log_burden_mean: finite static adds a constant compartment", {
  ld <- c(log(40), log(20))
  lg <- c(log(10), log(50))
  baseline <- 100
  static_log_level <- log(30)  # pi_static-scaled constant compartment
  data <- list(
    n_visits = 2L,
    patient_states = cbind(ld, lg),
    baseline = baseline,
    static_log_level = static_log_level
  )
  fit <- test_stan_function(
    here::here("tests/testthat/stan/test_calc_log_burden_mean_all.stan"),
    data
  )
  out_3 <- as.numeric(fit$draws("out_3arg", format = "draws_matrix"))
  expected <- log(baseline) + log(exp(ld) + exp(lg) + exp(static_log_level))
  # Tolerance is set by Stan's CSV output precision (~6-7 sig figs), not the
  # math: out_3 comes from the fit's text-serialized draws, expected from R.
  expect_equal(out_3, expected, tolerance = 1e-6)
})
