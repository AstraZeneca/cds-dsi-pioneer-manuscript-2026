library(testthat)

test_that("r_gp_exp_quad_cov: diagonal = alpha^2 + sigma, symmetric, PD", {
  x <- c(0.0, 1.0, 3.0); alpha <- 1.5; rho <- 2.0; sigma <- 0.1
  K <- r_gp_exp_quad_cov(x, alpha, rho, sigma)
  expect_equal(K[1, 1], alpha^2 + sigma, tolerance = 1e-9)
  expect_equal(K, t(K))
  expect_true(all(eigen(K, only.values = TRUE)$values > 0))
})

test_that("r_gp_matern32_cov: diagonal = alpha^2 + sigma, off-diagonal < diagonal", {
  x <- c(0.0, 1.0); alpha <- 1.0; rho <- 1.0; sigma <- 0.01
  K <- r_gp_matern32_cov(x, alpha, rho, sigma)
  expect_equal(K[1, 1], alpha^2 + sigma, tolerance = 1e-9)
  expect_true(K[1, 2] < K[1, 1])
})

test_that("r_gp_matern52_cov: symmetric and PD", {
  x <- c(0.0, 0.5, 2.0); alpha <- 0.8; rho <- 1.5; sigma <- 0.05
  K <- r_gp_matern52_cov(x, alpha, rho, sigma)
  expect_equal(K, t(K))
  expect_true(all(eigen(K, only.values = TRUE)$values > 0))
})

test_that("r_scale_process_sd: first row = process_sd * 1, delta_t applied", {
  tp <- c(0.0, 1.0, 3.0, 7.0)
  psd <- c(0.5, 0.8)
  out <- r_scale_process_sd(tp, psd)
  expect_equal(out[1, ], psd * 1.0, tolerance = 1e-9)
  expect_equal(out[2, ], psd * sqrt(1.0), tolerance = 1e-9)
  expect_equal(out[3, ], psd * sqrt(2.0), tolerance = 1e-9)
  expect_equal(out[4, ], psd * sqrt(4.0), tolerance = 1e-9)
})

test_that("r_gp_conditional_mean_zero: matches direct formula", {
  x_obs  <- c(0.0, 2.0, 4.0)
  x_pred <- c(1.0, 3.0)
  alpha <- 1.0; rho <- 2.0; delta <- 1e-6
  y_obs <- c(0.3, -0.2, 0.5)
  K_oo <- r_gp_exp_quad_cov(x_obs, alpha, rho, delta)
  K_po <- outer(x_pred, x_obs, function(a, b) alpha^2 * exp(-0.5 * ((a - b) / rho)^2))
  L_oo <- t(chol(K_oo))
  result <- r_gp_conditional_mean_zero(y_obs, L_oo, K_po)
  direct <- as.vector(K_po %*% solve(K_oo, y_obs))
  expect_equal(result, direct, tolerance = 1e-8)
})

test_that("r_multi_normal_cholesky_lcdf: identity cov = product of std normals", {
  y <- c(1.0, 0.5, -0.5)
  mu <- rep(0.0, 3L)
  L  <- diag(3L)
  result   <- r_multi_normal_cholesky_lcdf(y, mu, L)
  expected <- sum(pnorm(y, log.p = TRUE))
  expect_equal(result, expected, tolerance = 1e-9)
})
