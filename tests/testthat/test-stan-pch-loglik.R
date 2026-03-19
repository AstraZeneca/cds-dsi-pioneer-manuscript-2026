library(testthat)
library(here)

test_that("calc_pch_loglik: right-censored patients", {
  n <- 3L; max_t <- 10L
  lcs       <- matrix(-0.05, nrow = n, ncol = max_t)
  last_surv <- c(3L, 7L, 10L)
  rc        <- c(1L, 1L, 1L)
  ic        <- c(0L, 0L, 0L)

  stan_data <- list(
    N                        = n,
    MAX_T                    = max_t,
    N_EXIT                   = 1L,
    last_surv_week           = last_surv,
    exit_event               = rep(1L, n),
    right_censored           = rc,
    interval_censored        = ic,
    ignore_interval_censoring = 0L,
    log_cond_prob_surv       = array(lcs, dim = c(1L, n, max_t)),
    start_from               = rep(1L, n),
    end_at                   = rep(max_t, n)
  )
  fit <- test_stan_function(
    here("tests", "testthat", "stan", "test_pch_loglik_all.stan"),
    stan_data
  )
  d <- posterior::as_draws_df(fit$draws())

  expected <- r_calc_pch_loglik_simple(last_surv, rc, lcs)
  for (i in seq_len(n)) {
    expect_equal(
      get_stan_val(d, "lp_full", i), expected[i], tolerance = 1e-6,
      label = paste0("lp_full[", i, "] right-censored")
    )
    expect_equal(
      get_stan_val(d, "lp_single", i), expected[i], tolerance = 1e-6,
      label = paste0("lp_single[", i, "] right-censored")
    )
  }
})

test_that("calc_pch_loglik: event patients (no censoring)", {
  n <- 2L; max_t <- 8L
  lcs       <- matrix(-0.1, nrow = n, ncol = max_t)
  last_surv <- c(3L, 6L)
  rc        <- c(0L, 0L)
  ic        <- c(0L, 0L)

  stan_data <- list(
    N                        = n,
    MAX_T                    = max_t,
    N_EXIT                   = 1L,
    last_surv_week           = last_surv,
    exit_event               = rep(1L, n),
    right_censored           = rc,
    interval_censored        = ic,
    ignore_interval_censoring = 0L,
    log_cond_prob_surv       = array(lcs, dim = c(1L, n, max_t)),
    start_from               = rep(1L, n),
    end_at                   = rep(max_t, n)
  )
  fit <- test_stan_function(
    here("tests", "testthat", "stan", "test_pch_loglik_all.stan"),
    stan_data
  )
  d <- posterior::as_draws_df(fit$draws())

  expected <- r_calc_pch_loglik_simple(last_surv, rc, lcs)
  for (i in seq_len(n)) {
    expect_equal(
      get_stan_val(d, "lp_full", i), expected[i], tolerance = 1e-5,
      label = paste0("lp_full[", i, "] event")
    )
  }
})

test_that("calc_pch_loglik: mixed right-censored and event patients", {
  n <- 3L; max_t <- 6L
  lcs       <- matrix(-0.08, nrow = n, ncol = max_t)
  last_surv <- c(2L, 4L, 6L)
  rc        <- c(1L, 0L, 1L)
  ic        <- c(0L, 0L, 0L)

  stan_data <- list(
    N                        = n,
    MAX_T                    = max_t,
    N_EXIT                   = 1L,
    last_surv_week           = last_surv,
    exit_event               = rep(1L, n),
    right_censored           = rc,
    interval_censored        = ic,
    ignore_interval_censoring = 0L,
    log_cond_prob_surv       = array(lcs, dim = c(1L, n, max_t)),
    start_from               = rep(1L, n),
    end_at                   = rep(max_t, n)
  )
  fit <- test_stan_function(
    here("tests", "testthat", "stan", "test_pch_loglik_all.stan"),
    stan_data
  )
  d <- posterior::as_draws_df(fit$draws())

  expected <- r_calc_pch_loglik_simple(last_surv, rc, lcs)
  for (i in seq_len(n)) {
    expect_equal(
      get_stan_val(d, "lp_full", i), expected[i], tolerance = 1e-5,
      label = paste0("lp_full[", i, "] mixed")
    )
  }
})

test_that("pch_lpmf equals sum(calc_pch_loglik)", {
  n <- 3L; max_t <- 6L
  lcs       <- matrix(-0.08, nrow = n, ncol = max_t)
  last_surv <- c(2L, 4L, 6L)
  rc        <- c(1L, 0L, 1L)
  ic        <- c(0L, 0L, 0L)

  stan_data <- list(
    N                        = n,
    MAX_T                    = max_t,
    N_EXIT                   = 1L,
    last_surv_week           = last_surv,
    exit_event               = rep(1L, n),
    right_censored           = rc,
    interval_censored        = ic,
    ignore_interval_censoring = 0L,
    log_cond_prob_surv       = array(lcs, dim = c(1L, n, max_t)),
    start_from               = rep(1L, n),
    end_at                   = rep(max_t, n)
  )
  fit <- test_stan_function(
    here("tests", "testthat", "stan", "test_pch_loglik_all.stan"),
    stan_data
  )
  d <- posterior::as_draws_df(fit$draws())

  expected_sum <- sum(r_calc_pch_loglik_simple(last_surv, rc, lcs))
  expect_equal(get_stan_val(d, "lpmf_val"), expected_sum, tolerance = 1e-5)

  # Also verify lp_full components sum to lpmf_val
  stan_sum <- sum(sapply(seq_len(n), function(i) get_stan_val(d, "lp_full", i)))
  expect_equal(get_stan_val(d, "lpmf_val"), stan_sum, tolerance = 1e-7)
})

test_that("lp_no_end matches lp_full when end_at = max_t for all patients", {
  n <- 2L; max_t <- 5L
  lcs       <- matrix(-0.12, nrow = n, ncol = max_t)
  last_surv <- c(3L, 5L)
  rc        <- c(0L, 1L)
  ic        <- c(0L, 0L)

  stan_data <- list(
    N                        = n,
    MAX_T                    = max_t,
    N_EXIT                   = 1L,
    last_surv_week           = last_surv,
    exit_event               = rep(1L, n),
    right_censored           = rc,
    interval_censored        = ic,
    ignore_interval_censoring = 0L,
    log_cond_prob_surv       = array(lcs, dim = c(1L, n, max_t)),
    start_from               = rep(1L, n),
    end_at                   = rep(max_t, n)
  )
  fit <- test_stan_function(
    here("tests", "testthat", "stan", "test_pch_loglik_all.stan"),
    stan_data
  )
  d <- posterior::as_draws_df(fit$draws())

  for (i in seq_len(n)) {
    expect_equal(
      get_stan_val(d, "lp_no_end", i),
      get_stan_val(d, "lp_full", i),
      tolerance = 1e-9,
      label = paste0("lp_no_end[", i, "] == lp_full[", i, "]")
    )
  }
})
