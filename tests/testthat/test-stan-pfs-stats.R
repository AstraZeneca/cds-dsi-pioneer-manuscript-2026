library(testthat)
library(cmdstanr)
library(posterior)

stan_file <- here::here("tests", "testthat", "stan", "test_pfs_stats_all.stan")

test_that("pfs stats Stan functions match R oracles: main cases", {
  log_cond_surv  <- c(-0.1, -0.2, -0.3)
  surv_time      <- c(2L, 4L, 6L, 8L, 10L)
  last_surv_time <- 12L
  sq_p           <- c(0.25, 0.5)
  km_survival    <- c(1.0, 0.8, 0.6, 0.4, 0.2, 0.1)
  km_p           <- c(0.25, 0.5)
  pfs_n_query    <- 2.5

  stan_data <- list(
    T_lmep         = length(log_cond_surv),
    log_cond_surv  = log_cond_surv,
    N_sq           = length(surv_time),
    surv_time      = surv_time,
    last_surv_time = last_surv_time,
    P_sq           = length(sq_p),
    sq_p           = sq_p,
    T_km           = length(km_survival),
    km_survival    = km_survival,
    P_km           = length(km_p),
    km_p           = km_p,
    pfs_n_query    = pfs_n_query
  )

  fit      <- test_stan_function(stan_file, data = stan_data)
  draws_df <- posterior::as_draws_df(fit$draws())
  gv       <- function(var, ...) get_stan_val(draws_df, var, ...)

  # ---- calculate_log_marginal_exit_prob ----
  r_lmep <- r_calculate_log_marginal_exit_prob(log_cond_surv)

  expect_equal(gv("lmep", 1), r_lmep[1], tolerance = 1e-6,
    label = "lmep[1]")
  expect_equal(gv("lmep", 2), r_lmep[2], tolerance = 1e-6,
    label = "lmep[2]")
  expect_equal(gv("lmep", 3), r_lmep[3], tolerance = 1e-6,
    label = "lmep[3]")

  # ---- survival_quantiles ----
  r_sq <- r_survival_quantiles(surv_time, last_surv_time, sq_p)

  expect_equal(gv("sq_quantiles", 1), r_sq$quantiles[1], tolerance = 1e-6,
    label = "sq_quantiles[1] (p=0.25)")
  expect_equal(gv("sq_quantiles", 2), r_sq$quantiles[2], tolerance = 1e-6,
    label = "sq_quantiles[2] (p=0.5)")
  expect_equal(gv("sq_cannot_calc", 1), r_sq$cannot_calculate[1],
    label = "sq_cannot_calc[1]")
  expect_equal(gv("sq_cannot_calc", 2), r_sq$cannot_calculate[2],
    label = "sq_cannot_calc[2]")

  # ---- survival_median == survival_quantiles(p=0.5) ----
  r_sq_med <- r_survival_quantiles(surv_time, last_surv_time, c(0.5))
  expect_equal(gv("surv_med"), r_sq_med$quantiles[1], tolerance = 1e-6,
    label = "surv_med matches survival_quantiles(p=0.5)")
  expect_equal(gv("surv_med_flag"), r_sq_med$cannot_calculate[1],
    label = "surv_med_flag matches survival_quantiles(p=0.5) flag")

  # ---- km_quantiles ----
  r_km <- r_km_quantiles(km_survival, km_p)

  expect_equal(gv("km_q", 1), r_km$quantiles[1], tolerance = 1e-6,
    label = "km_q[1] (p=0.25)")
  expect_equal(gv("km_q", 2), r_km$quantiles[2], tolerance = 1e-6,
    label = "km_q[2] (p=0.5)")
  expect_equal(gv("km_cannot_calc", 1), r_km$cannot_calculate[1],
    label = "km_cannot_calc[1]")
  expect_equal(gv("km_cannot_calc", 2), r_km$cannot_calculate[2],
    label = "km_cannot_calc[2]")

  # ---- km_median == km_quantiles(p=0.5) ----
  r_km_med <- r_km_quantiles(km_survival, c(0.5))
  expect_equal(gv("km_med"), r_km_med$quantiles[1], tolerance = 1e-6,
    label = "km_med matches km_quantiles(p=0.5)")
  expect_equal(gv("km_med_flag"), r_km_med$cannot_calculate[1],
    label = "km_med_flag matches km_quantiles(p=0.5) flag")

  # ---- calc_km_pfs_n ----
  r_pfs_n <- r_calc_km_pfs_n(km_survival, pfs_n_query)
  expect_equal(gv("pfs_n"), r_pfs_n, tolerance = 1e-6,
    label = "calc_km_pfs_n(km_survival, 2.5)")
})

test_that("pfs stats Stan functions: edge case surv_time all beyond last_surv_time", {
  log_cond_surv  <- c(-0.1)
  surv_time_edge <- c(5L, 6L, 7L, 8L, 9L)
  last_surv_edge <- 3L
  sq_p_edge      <- c(0.5)
  km_survival    <- c(1.0, 0.8, 0.5, 0.2)
  km_p           <- c(0.5)
  pfs_n_query    <- 1.0

  stan_data <- list(
    T_lmep         = length(log_cond_surv),
    log_cond_surv  = log_cond_surv,
    N_sq           = length(surv_time_edge),
    surv_time      = surv_time_edge,
    last_surv_time = last_surv_edge,
    P_sq           = length(sq_p_edge),
    sq_p           = sq_p_edge,
    T_km           = length(km_survival),
    km_survival    = km_survival,
    P_km           = length(km_p),
    km_p           = km_p,
    pfs_n_query    = pfs_n_query
  )

  fit      <- test_stan_function(stan_file, data = stan_data)
  draws_df <- posterior::as_draws_df(fit$draws())
  gv       <- function(var, ...) get_stan_val(draws_df, var, ...)

  # survival_quantiles: all times > last_surv_time → cannot_calculate = 1
  expect_equal(gv("sq_cannot_calc", 1), 1L,
    label = "sq_cannot_calc[1] = 1 when all surv_time > last_surv_time")
  expect_equal(gv("sq_quantiles", 1), 0.0,
    label = "sq_quantiles[1] = 0 when not calculable")
})
