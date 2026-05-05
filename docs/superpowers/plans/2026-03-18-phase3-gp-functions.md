# Plan: Phase 3 — GP Functions Coverage

**Date:** 2026-03-18
**Branch:** karim/tests
**Source file:** `stan/gp.stanfunctions`
**Target assertions:** ~70 new assertions

## Background

`gp.stanfunctions` provides the Gaussian process infrastructure used throughout the state-space
model. All 18 exported functions (across 3 kernel families, 6 `calc_gp_pred` overloads, the
conditional posterior pair, and the multivariate normal utilities) are currently untested.

## Functions to cover

| Group | Functions |
|-------|-----------|
| Kernel wrappers (with nugget) | `gp_exp_quad_cov`, `gp_matern32_cov`, `gp_matern52_cov` |
| Cholesky wrappers | `gp_exp_quad_cholesky_cov`, `gp_matern32_cholesky_cov`, `gp_matern52_cholesky_cov` |
| GP predictor | `calc_gp_pred` (6 overloads: `vector`/`row_vector` × intercept/scalar/none, plus Kronecker) |
| Process SD scaling | `scale_process_sd` |
| Conditional posterior | `gp_conditional_mean` (2 overloads), `gp_conditional_cov`, `gp_conditional` (2 overloads) |
| Multivariate normal | `multi_normal_cholesky_lcdf` (2 overloads) |

`multi_normal_rng` is excluded (RNG; verifiable only probabilistically — covered separately if needed).

---

## Task 1 — R helper `helper-gp.R`

**File:** `tests/testthat/helper-gp.R`
**Purpose:** R oracle functions that mirror the Stan implementations exactly.

```r
# Kernel: base Stan built-in + nugget sigma on diagonal
r_gp_exp_quad_cov <- function(x, alpha, rho, sigma) {
  n <- length(x)
  K <- matrix(0, n, n)
  for (i in seq_len(n)) for (j in seq_len(n)) {
    K[i, j] <- alpha^2 * exp(-0.5 * ((x[i] - x[j]) / rho)^2)
  }
  diag(K) <- diag(K) + sigma
  K
}

r_gp_matern32_cov <- function(x, alpha, rho, sigma) {
  n <- length(x)
  K <- matrix(0, n, n)
  for (i in seq_len(n)) for (j in seq_len(n)) {
    r <- sqrt(3) * abs(x[i] - x[j]) / rho
    K[i, j] <- alpha^2 * (1 + r) * exp(-r)
  }
  diag(K) <- diag(K) + sigma
  K
}

r_gp_matern52_cov <- function(x, alpha, rho, sigma) {
  n <- length(x)
  K <- matrix(0, n, n)
  for (i in seq_len(n)) for (j in seq_len(n)) {
    r <- sqrt(5) * abs(x[i] - x[j]) / rho
    K[i, j] <- alpha^2 * (1 + r + r^2/3) * exp(-r)
  }
  diag(K) <- diag(K) + sigma
  K
}

r_scale_process_sd <- function(time_points, process_sd) {
  n <- length(time_points)
  out <- matrix(0, n, 2)
  for (t in seq_len(n)) {
    delta_t <- if (t > 1) time_points[t] - time_points[t - 1] else 1.0
    out[t, ] <- process_sd * sqrt(delta_t)
  }
  out
}

r_gp_conditional_mean <- function(mu_obs, mu_pred, y_obs, L_K_obs, K_pred_obs) {
  # Solves L_K \ (y - mu_obs) via forward/backward substitution, then L_K' \
  alpha <- backsolve(t(L_K_obs), forwardsolve(L_K_obs, y_obs - mu_obs))
  mu_pred + K_pred_obs %*% alpha
}

r_gp_conditional_cov <- function(L_K_obs, K_pred_obs, K_pred_pred, delta) {
  V <- forwardsolve(L_K_obs, t(K_pred_obs))   # L_K \ K(X,X*)
  K_pred_pred - t(V) %*% V + diag(delta, nrow(K_pred_pred))
}

r_multi_normal_cholesky_lcdf <- function(y, mu, L_Sigma) {
  z <- forwardsolve(L_Sigma, y - mu)
  sum(pnorm(z, log.p = TRUE))
}
```

**Test file:** `tests/testthat/test-helper-gp.R`
**Assertions (~20):**
- `r_gp_exp_quad_cov`: diagonal = alpha^2 + sigma; off-diag decays with distance
- `r_gp_matern32_cov`: diagonal = alpha^2 + sigma; symmetry check
- `r_gp_matern52_cov`: symmetry; positive-definite (all eigenvalues > 0)
- `r_scale_process_sd`: first row = process_sd × 1.0; second row = process_sd × sqrt(delta_t)
- `r_gp_conditional_mean`: matches direct formula mu_pred + K_xstar_x inv(K_xx) (y-mu)
- `r_gp_conditional_cov`: PD check; diagonal equals K_pred_pred diag - V^T V + delta
- `r_multi_normal_cholesky_lcdf`: independent case: sum of log(pnorm(z_i)); correlated case vs
  `mvtnorm::pmvnorm`

---

## Task 2 — Stan harness `test_gp_kernels_all.stan`

**File:** `tests/testthat/stan/test_gp_kernels_all.stan`

```stan
functions {
  #include gp.stanfunctions
}
data {
  int<lower=1> N;
  array[N] real x;
  real alpha;
  real rho;
  real sigma;
  real delta;
  // For conditional GP test
  int<lower=1> N_obs;
  int<lower=1> N_pred;
  array[N_obs] real x_obs;
  array[N_pred] real x_pred;
  vector[N_obs] y_obs;
  // For scale_process_sd
  int<lower=1> N_tp;
  array[N_tp] real time_points;
  vector[2] process_sd;
  // For multi_normal_cholesky_lcdf
  int<lower=1> D;
  vector[D] lcdf_y;
  real lcdf_mu;
}
generated quantities {
  // Kernel tests
  matrix[N, N] K_eq  = gp_exp_quad_cov(x, alpha, rho, sigma);
  matrix[N, N] K_m32 = gp_matern32_cov(x, alpha, rho, sigma);
  matrix[N, N] K_m52 = gp_matern52_cov(x, alpha, rho, sigma);
  // Cholesky tests
  matrix[N, N] L_eq  = gp_exp_quad_cholesky_cov(x, alpha, rho, delta);
  matrix[N, N] L_m32 = gp_matern32_cholesky_cov(x, alpha, rho, delta);
  matrix[N, N] L_m52 = gp_matern52_cholesky_cov(x, alpha, rho, delta);
  // Verify L L' reconstructs K
  matrix[N, N] K_from_L_eq = L_eq * L_eq';
  // calc_gp_pred (vector overload with explicit intercept)
  vector[N] eta_vec = rep_vector(0.5, N);
  vector[N] gp_pred_vec = calc_gp_pred(x, 1.0, alpha, rho, delta, eta_vec);
  // scale_process_sd
  matrix[N_tp, 2] scaled_sd = scale_process_sd(time_points, process_sd);
  // gp_conditional_mean + cov
  matrix[N_obs, N_obs] K_obs_obs = gp_exp_quad_cov(x_obs, alpha, rho, delta);
  matrix[N_pred, N_obs] K_pred_obs = gp_exp_quad_cov(x_pred, x_obs, alpha, rho);
  matrix[N_pred, N_pred] K_pred_pred = gp_exp_quad_cov(x_pred, alpha, rho, delta);
  matrix[N_obs, N_obs] L_obs = cholesky_decompose(K_obs_obs);
  vector[N_pred] cond_mean = gp_conditional_mean(y_obs, L_obs, K_pred_obs);
  matrix[N_pred, N_pred] cond_cov = gp_conditional_cov(L_obs, K_pred_obs, K_pred_pred, delta);
  // multi_normal_cholesky_lcdf
  matrix[D, D] L_id = diag_matrix(rep_vector(1.0, D));
  real lcdf_val = multi_normal_cholesky_lcdf(lcdf_y, lcdf_mu, L_id);
}
```

**Test file:** `tests/testthat/test-stan-gp-kernels.R`

```r
library(testthat)

test_that("Stan GP kernel functions match R oracles", {
  x     <- c(0.0, 1.0, 3.0, 6.0)
  alpha <- 1.5; rho <- 2.0; sigma <- 0.1; delta <- 1e-6
  x_obs  <- c(0.0, 2.0, 4.0)
  x_pred <- c(1.0, 3.0)
  y_obs  <- c(0.3, -0.2, 0.5)

  stan_data <- list(
    N = length(x), x = x, alpha = alpha, rho = rho, sigma = sigma, delta = delta,
    N_obs = length(x_obs), N_pred = length(x_pred),
    x_obs = x_obs, x_pred = x_pred, y_obs = y_obs,
    N_tp = 4L, time_points = c(0.0, 1.0, 3.0, 7.0),
    process_sd = c(0.5, 0.8),
    D = 3L, lcdf_y = c(1.0, 0.5, -0.5), lcdf_mu = 0.0
  )

  fit <- test_stan_function("test_gp_kernels_all", stan_data)
  d   <- posterior::as_draws_df(fit$draws())

  # --- Kernel: diagonal = alpha^2 + sigma ---
  K_eq_r <- r_gp_exp_quad_cov(x, alpha, rho, sigma)
  for (i in seq_along(x)) {
    expect_equal(get_stan_val(d, "K_eq", i, i), K_eq_r[i, i], tolerance = 1e-6,
      label = paste0("K_eq[", i, ",", i, "]"))
  }

  # --- Cholesky: L L' ≈ K ---
  L_eq_r <- chol(r_gp_exp_quad_cov(x, alpha, rho, delta))
  K_from_L_r <- t(L_eq_r) %*% L_eq_r   # R chol is upper-tri; Stan is lower-tri
  for (i in 1:length(x)) {
    expect_equal(get_stan_val(d, "K_from_L_eq", i, i),
                 r_gp_exp_quad_cov(x, alpha, rho, delta)[i, i],
                 tolerance = 1e-5, label = paste0("K_from_L_eq[", i, ",", i, "]"))
  }

  # --- scale_process_sd: first row = process_sd × 1 ---
  expect_equal(get_stan_val(d, "scaled_sd", 1, 1), 0.5 * 1.0, tolerance = 1e-6)
  expect_equal(get_stan_val(d, "scaled_sd", 1, 2), 0.8 * 1.0, tolerance = 1e-6)
  expect_equal(get_stan_val(d, "scaled_sd", 2, 1), 0.5 * sqrt(1.0), tolerance = 1e-6)
  expect_equal(get_stan_val(d, "scaled_sd", 3, 1), 0.5 * sqrt(2.0), tolerance = 1e-6)

  # --- gp_conditional_mean ---
  K_oo <- r_gp_exp_quad_cov(x_obs, alpha, rho, delta)
  K_po <- outer(x_pred, x_obs, function(a, b) alpha^2 * exp(-0.5 * ((a - b)/rho)^2))
  L_oo <- t(chol(K_oo))
  cond_mean_r <- r_gp_conditional_mean(rep(0, 3), rep(0, 2), y_obs, L_oo, K_po)
  expect_equal(get_stan_val(d, "cond_mean", 1), cond_mean_r[1], tolerance = 1e-5)
  expect_equal(get_stan_val(d, "cond_mean", 2), cond_mean_r[2], tolerance = 1e-5)

  # --- multi_normal_cholesky_lcdf: identity covariance = product of std normals ---
  expected_lcdf <- sum(pnorm(c(1.0, 0.5, -0.5), log.p = TRUE))
  expect_equal(get_stan_val(d, "lcdf_val"), expected_lcdf, tolerance = 1e-6)
})
```

**Assertion count:** ~25 Stan assertions

---

## Summary

| Task | File | Assertions |
|------|------|-----------|
| 1 | `helper-gp.R` + `test-helper-gp.R` | ~20 R |
| 2 | `test_gp_kernels_all.stan` + `test-stan-gp-kernels.R` | ~25 Stan |
| **Total** | | **~45 new assertions** |
