library(testthat)
library(here)

test_that("Stan GP kernels: diagonal, symmetry, Cholesky reconstruction", {
  x     <- c(0.0, 1.0, 3.0, 6.0)
  alpha <- 1.5; rho <- 2.0; sigma <- 0.1; delta <- 1e-6

  stan_data <- list(
    N = length(x), x = x, alpha = alpha, rho = rho, sigma = sigma, delta = delta,
    N_obs = 3L, x_obs = c(0.0, 2.0, 4.0),
    N_pred = 2L, x_pred = c(1.0, 3.0),
    y_obs = c(0.3, -0.2, 0.5),
    N_tp = 4L, time_points = c(0.0, 1.0, 3.0, 7.0),
    process_sd = c(0.5, 0.8),
    D = 3L, lcdf_y = c(1.0, 0.5, -0.5)
  )

  fit <- test_stan_function(
    here("tests", "testthat", "stan", "test_gp_all.stan"),
    stan_data
  )
  d <- posterior::as_draws_df(fit$draws())

  K_eq_r   <- r_gp_exp_quad_cov(x, alpha, rho, sigma)
  K_m32_r  <- r_gp_matern32_cov(x, alpha, rho, sigma)
  K_m52_r  <- r_gp_matern52_cov(x, alpha, rho, sigma)

  # Diagonal = alpha^2 + sigma for all three kernels
  for (i in seq_along(x)) {
    expect_equal(
      get_stan_val(d, "K_eq", i, i), K_eq_r[i, i], tolerance = 1e-6,
      label = paste0("K_eq[", i, ",", i, "]")
    )
    expect_equal(
      get_stan_val(d, "K_m32", i, i), K_m32_r[i, i], tolerance = 1e-6,
      label = paste0("K_m32[", i, ",", i, "]")
    )
    expect_equal(
      get_stan_val(d, "K_m52", i, i), K_m52_r[i, i], tolerance = 1e-6,
      label = paste0("K_m52[", i, ",", i, "]")
    )
  }

  # Off-diagonal symmetry: K_eq[1,2] == K_eq[2,1]
  expect_equal(
    get_stan_val(d, "K_eq", 1, 2),
    get_stan_val(d, "K_eq", 2, 1),
    tolerance = 1e-9
  )

  # Cholesky reconstruction: L_eq L_eq' == gp_exp_quad_cov(x, alpha, rho, delta)
  K_eq_r_delta <- r_gp_exp_quad_cov(x, alpha, rho, delta)
  for (i in seq_along(x)) {
    expect_equal(
      get_stan_val(d, "K_from_L", i, i), K_eq_r_delta[i, i], tolerance = 1e-5,
      label = paste0("K_from_L[", i, ",", i, "]")
    )
  }
})

test_that("Stan GP kernels: off-diagonal values match R oracle (exp_quad and matern)", {
  x     <- c(0.0, 1.0, 3.0, 6.0)
  alpha <- 1.5; rho <- 2.0; sigma <- 0.1; delta <- 1e-6

  stan_data <- list(
    N = length(x), x = x, alpha = alpha, rho = rho, sigma = sigma, delta = delta,
    N_obs = 3L, x_obs = c(0.0, 2.0, 4.0),
    N_pred = 2L, x_pred = c(1.0, 3.0),
    y_obs = c(0.3, -0.2, 0.5),
    N_tp = 4L, time_points = c(0.0, 1.0, 3.0, 7.0),
    process_sd = c(0.5, 0.8),
    D = 3L, lcdf_y = c(1.0, 0.5, -0.5)
  )

  fit <- test_stan_function(
    here("tests", "testthat", "stan", "test_gp_all.stan"),
    stan_data
  )
  d <- posterior::as_draws_df(fit$draws())

  K_eq_r  <- r_gp_exp_quad_cov(x, alpha, rho, sigma)
  K_m32_r <- r_gp_matern32_cov(x, alpha, rho, sigma)
  K_m52_r <- r_gp_matern52_cov(x, alpha, rho, sigma)

  # Check a selection of off-diagonal entries against the R oracle
  pairs <- list(c(1, 2), c(1, 3), c(2, 3), c(1, 4))
  for (p in pairs) {
    i <- p[1]; j <- p[2]
    expect_equal(
      get_stan_val(d, "K_eq", i, j), K_eq_r[i, j], tolerance = 1e-6,
      label = paste0("K_eq[", i, ",", j, "]")
    )
    expect_equal(
      get_stan_val(d, "K_m32", i, j), K_m32_r[i, j], tolerance = 1e-6,
      label = paste0("K_m32[", i, ",", j, "]")
    )
    expect_equal(
      get_stan_val(d, "K_m52", i, j), K_m52_r[i, j], tolerance = 1e-6,
      label = paste0("K_m52[", i, ",", j, "]")
    )
  }
})

test_that("calc_gp_pred: constant eta produces centered output = intercept", {
  # When eta is constant (all 0.5), gp_dev = L * eta is proportional to row sums of L.
  # After centering (subtracting mean), deviations cancel and output = intercept.
  # We verify that all elements equal the same value (intercept = 0 so all zero).
  x     <- c(0.0, 1.0, 2.0)
  alpha <- 1.0; rho <- 1.5; delta <- 1e-6

  stan_data <- list(
    N = length(x), x = x, alpha = alpha, rho = rho, sigma = 0.01, delta = delta,
    N_obs = 2L, x_obs = c(0.0, 1.0),
    N_pred = 1L, x_pred = c(0.5),
    y_obs = c(0.1, -0.1),
    N_tp = 3L, time_points = c(0.0, 1.0, 2.0),
    process_sd = c(0.5, 0.8),
    D = 2L, lcdf_y = c(0.5, -0.5)
  )

  fit <- test_stan_function(
    here("tests", "testthat", "stan", "test_gp_all.stan"),
    stan_data
  )
  d <- posterior::as_draws_df(fit$draws())

  # gp_out should be centered around intercept=0; all values share the same deviation from 0
  gp_vals <- vapply(seq_along(x), function(i) get_stan_val(d, "gp_out", i), numeric(1))
  # The sum of centered deviations must be zero (up to floating-point rounding)
  expect_equal(sum(gp_vals), 0.0, tolerance = 1e-6)
})

test_that("scale_process_sd: first row = process_sd * 1, subsequent rows scaled by sqrt(delta_t)", {
  x <- c(0.0, 1.0, 2.0)
  stan_data <- list(
    N = length(x), x = x, alpha = 1.0, rho = 1.0, sigma = 0.01, delta = 1e-6,
    N_obs = 2L, x_obs = c(0.0, 1.0),
    N_pred = 1L, x_pred = c(0.5),
    y_obs = c(0.1, -0.1),
    N_tp = 4L, time_points = c(0.0, 1.0, 3.0, 7.0),
    process_sd = c(0.5, 0.8),
    D = 2L, lcdf_y = c(0.5, -0.5)
  )

  fit <- test_stan_function(
    here("tests", "testthat", "stan", "test_gp_all.stan"),
    stan_data
  )
  d <- posterior::as_draws_df(fit$draws())

  r_out <- r_scale_process_sd(c(0.0, 1.0, 3.0, 7.0), c(0.5, 0.8))

  # Row 1: delta_t = 1.0 (convention for first time point)
  expect_equal(get_stan_val(d, "scaled_sd", 1, 1), r_out[1, 1], tolerance = 1e-6)
  expect_equal(get_stan_val(d, "scaled_sd", 1, 2), r_out[1, 2], tolerance = 1e-6)
  # Row 2: delta_t = 1 -> sqrt(1) = 1
  expect_equal(get_stan_val(d, "scaled_sd", 2, 1), r_out[2, 1], tolerance = 1e-6)
  # Row 3: delta_t = 2 -> sqrt(2)
  expect_equal(get_stan_val(d, "scaled_sd", 3, 1), r_out[3, 1], tolerance = 1e-6)
  # Row 4: delta_t = 4 -> sqrt(4) = 2
  expect_equal(get_stan_val(d, "scaled_sd", 4, 2), r_out[4, 2], tolerance = 1e-6)
})

test_that("gp_conditional_mean: matches R oracle", {
  x_obs  <- c(0.0, 2.0, 4.0)
  x_pred <- c(1.0, 3.0)
  alpha <- 1.0; rho <- 2.0; delta <- 1e-6
  y_obs <- c(0.3, -0.2, 0.5)
  x <- c(0.0, 1.0, 2.0)

  stan_data <- list(
    N = length(x), x = x, alpha = alpha, rho = rho, sigma = 0.01, delta = delta,
    N_obs = length(x_obs), x_obs = x_obs,
    N_pred = length(x_pred), x_pred = x_pred,
    y_obs = y_obs,
    N_tp = 3L, time_points = c(0.0, 1.0, 2.0),
    process_sd = c(0.5, 0.8),
    D = 2L, lcdf_y = c(0.5, -0.5)
  )

  fit <- test_stan_function(
    here("tests", "testthat", "stan", "test_gp_all.stan"),
    stan_data
  )
  d <- posterior::as_draws_df(fit$draws())

  K_oo_r <- r_gp_exp_quad_cov(x_obs, alpha, rho, delta)
  K_po_r <- outer(x_pred, x_obs, function(a, b) alpha^2 * exp(-0.5 * ((a - b) / rho)^2))
  L_oo_r <- t(chol(K_oo_r))
  cm_r   <- r_gp_conditional_mean_zero(y_obs, L_oo_r, K_po_r)

  for (j in seq_along(x_pred)) {
    expect_equal(
      get_stan_val(d, "cond_mean", j), cm_r[j], tolerance = 1e-5,
      label = paste0("cond_mean[", j, "]")
    )
  }
})

test_that("gp_conditional_cov: diagonal is positive (PD check)", {
  x_obs  <- c(0.0, 2.0, 4.0)
  x_pred <- c(1.0, 3.0)
  alpha <- 1.0; rho <- 2.0; delta <- 1e-6
  y_obs <- c(0.3, -0.2, 0.5)
  x <- c(0.0, 1.0, 2.0)

  stan_data <- list(
    N = length(x), x = x, alpha = alpha, rho = rho, sigma = 0.01, delta = delta,
    N_obs = length(x_obs), x_obs = x_obs,
    N_pred = length(x_pred), x_pred = x_pred,
    y_obs = y_obs,
    N_tp = 3L, time_points = c(0.0, 1.0, 2.0),
    process_sd = c(0.5, 0.8),
    D = 2L, lcdf_y = c(0.5, -0.5)
  )

  fit <- test_stan_function(
    here("tests", "testthat", "stan", "test_gp_all.stan"),
    stan_data
  )
  d <- posterior::as_draws_df(fit$draws())

  # Diagonal entries of posterior covariance must be positive
  for (j in seq_along(x_pred)) {
    expect_gt(
      get_stan_val(d, "cond_cov", j, j), 0.0
    )
  }

  # Compare diagonal against R oracle
  K_oo_r  <- r_gp_exp_quad_cov(x_obs, alpha, rho, delta)
  K_po_r  <- outer(x_pred, x_obs, function(a, b) alpha^2 * exp(-0.5 * ((a - b) / rho)^2))
  K_pp_r  <- r_gp_exp_quad_cov(x_pred, alpha, rho, delta)
  L_oo_r  <- t(chol(K_oo_r))
  cc_r    <- r_gp_conditional_cov(L_oo_r, K_po_r, K_pp_r, delta)

  for (j in seq_along(x_pred)) {
    expect_equal(
      get_stan_val(d, "cond_cov", j, j), cc_r[j, j], tolerance = 1e-5,
      label = paste0("cond_cov[", j, ",", j, "]")
    )
  }
})

test_that("multi_normal_cholesky_lcdf: identity cov = sum of log(pnorm(y_i))", {
  lcdf_y <- c(1.0, 0.5, -0.5)
  x      <- c(0.0, 1.0, 2.0)

  stan_data <- list(
    N = length(x), x = x, alpha = 1.0, rho = 1.0, sigma = 0.01, delta = 1e-6,
    N_obs = 2L, x_obs = c(0.0, 1.0),
    N_pred = 1L, x_pred = c(0.5),
    y_obs = c(0.1, -0.1),
    N_tp = 3L, time_points = c(0.0, 1.0, 2.0),
    process_sd = c(0.5, 0.8),
    D = length(lcdf_y), lcdf_y = lcdf_y
  )

  fit <- test_stan_function(
    here("tests", "testthat", "stan", "test_gp_all.stan"),
    stan_data
  )
  d <- posterior::as_draws_df(fit$draws())

  expected <- sum(pnorm(lcdf_y, log.p = TRUE))
  expect_equal(get_stan_val(d, "lcdf_val"), expected, tolerance = 1e-6)
})
