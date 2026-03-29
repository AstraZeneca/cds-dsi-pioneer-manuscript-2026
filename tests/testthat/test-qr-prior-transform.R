library(purrr)
library(tibble)

test_that("transform_priors_to_qr_space round-trips correctly", {
  source(here::here("r/priors.R"))

  set.seed(42)

  # Create a design matrix with correlated covariates
  n <- 100
  p <- 4
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  # Introduce correlation: col 2 is correlated with col 1
  X[, 2] <- X[, 1] * 0.7 + X[, 2] * 0.3
  colnames(X) <- str_c("x", seq_len(p))

  # Elicited priors in original space
  coef_mean_orig <- c(0.5, -0.3, 0.0, 0.8)
  coef_sd_orig <- c(0.2, 0.15, 0.3, 0.25)

  # Transform to QR space
  qr_result <- transform_priors_to_qr_space(coef_mean_orig, coef_sd_orig, X)

  # Compute R_stan the same way Stan does
  R_stan <- qr.R(qr(X)) / sqrt(n - 1)

  # Verify R_stan matches

  expect_equal(qr_result$R_stan, R_stan)

  # Verify mean transformation: mu_qr = R * mu_orig
  expect_equal(qr_result$coef_mean_qr, as.vector(R_stan %*% coef_mean_orig))

  # Verify SD transformation: diag of R * diag(sd^2) * R'
  expected_var_qr <- diag(R_stan %*% diag(coef_sd_orig^2) %*% t(R_stan))
  expect_equal(qr_result$coef_sd_qr, sqrt(expected_var_qr))

  # Back-transform from QR space: beta = R^{-1} * beta_qr
  beta_qr <- qr_result$coef_mean_qr
  beta_recovered <- as.vector(solve(R_stan, beta_qr))

  expect_equal(beta_recovered, coef_mean_orig, tolerance = 1e-10)
})

test_that("transform_priors_to_qr_space handles isotropic prior (identity case)", {
  source(here::here("r/priors.R"))

  set.seed(123)
  n <- 50
  p <- 3
  # Orthogonal design matrix (no rotation effect on isotropic prior)
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(X) <- str_c("x", seq_len(p))

  # Isotropic prior: N(0, I) — should change under rotation
  coef_mean <- rep(0, p)
  coef_sd <- rep(1, p)

  qr_result <- transform_priors_to_qr_space(coef_mean, coef_sd, X)

  # Mean of 0 stays 0 under any linear transform

  expect_equal(qr_result$coef_mean_qr, rep(0, p), tolerance = 1e-10)

  # SDs may change (R is not orthogonal in general), but that's expected
  expect_true(all(qr_result$coef_sd_qr > 0))
})

test_that("get_tumor_priors uses QR-transformed priors when design matrix provided", {
  source(here::here("r/priors.R"))

  set.seed(42)
  n <- 50
  p <- 3
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(X) <- str_c("x", seq_len(p))

  coef_elicited <- list(
    coef_mean = c(0.5, -0.3, 0.0),
    coef_sd = c(0.2, 0.15, 0.3)
  )

  stan_data <- list(
    n_levels = 2L,
    n_covar = p,
    n_time_varying_covar = 0L,
    n_time_invariant_covar = p
  )

  # With design matrix: priors should be transformed
  priors_with <- get_tumor_priors(stan_data, coef_elicited, X)
  # Without design matrix: priors should be original
  priors_without <- get_tumor_priors(stan_data, coef_elicited, NULL)

  # Frac and init priors should differ when design matrix is provided
  expect_false(
    isTRUE(all.equal(priors_with$frac_coef_qr_pop_mean,
                     priors_without$frac_coef_qr_pop_mean))
  )
  expect_false(
    isTRUE(all.equal(priors_with$init_coef_qr_pop_mean,
                     priors_without$init_coef_qr_pop_mean))
  )

  # TR priors should remain unchanged (isotropic N(0,1))
  expect_equal(
    as.vector(priors_with$tr_coef_qr_pop_mean),
    as.vector(priors_without$tr_coef_qr_pop_mean)
  )
})
