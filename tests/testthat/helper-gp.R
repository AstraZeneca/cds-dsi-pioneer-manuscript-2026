# tests/testthat/helper-gp.R
#
# Oracle functions for gp.stanfunctions tests

r_gp_exp_quad_cov <- function(x, alpha, rho, sigma) {
  n <- length(x)
  K <- matrix(0.0, n, n)
  for (i in seq_len(n)) for (j in seq_len(n))
    K[i, j] <- alpha^2 * exp(-0.5 * ((x[i] - x[j]) / rho)^2)
  diag(K) <- diag(K) + sigma
  K
}

r_gp_matern32_cov <- function(x, alpha, rho, sigma) {
  n <- length(x)
  K <- matrix(0.0, n, n)
  for (i in seq_len(n)) for (j in seq_len(n)) {
    r <- sqrt(3) * abs(x[i] - x[j]) / rho
    K[i, j] <- alpha^2 * (1 + r) * exp(-r)
  }
  diag(K) <- diag(K) + sigma
  K
}

r_gp_matern52_cov <- function(x, alpha, rho, sigma) {
  n <- length(x)
  K <- matrix(0.0, n, n)
  for (i in seq_len(n)) for (j in seq_len(n)) {
    r <- sqrt(5) * abs(x[i] - x[j]) / rho
    K[i, j] <- alpha^2 * (1 + r + r^2 / 3) * exp(-r)
  }
  diag(K) <- diag(K) + sigma
  K
}

r_scale_process_sd <- function(time_points, process_sd) {
  # Returns [n_time x 2] matrix: row t = process_sd * sqrt(delta_t[t])
  # delta_t[1] = 1.0 (convention matching Stan implementation)
  n <- length(time_points)
  out <- matrix(0.0, n, 2L)
  for (t in seq_len(n)) {
    delta_t <- if (t > 1L) time_points[t] - time_points[t - 1L] else 1.0
    out[t, ] <- process_sd * sqrt(delta_t)
  }
  out
}

r_gp_conditional_mean_zero <- function(y_obs, L_K_obs, K_pred_obs) {
  # Zero-mean version: K_pred_obs * inv(K_obs) * y_obs
  # = K_pred_obs * L_K_obs^{-T} * L_K_obs^{-1} * y_obs
  alpha <- backsolve(t(L_K_obs), forwardsolve(L_K_obs, y_obs))
  as.vector(K_pred_obs %*% alpha)
}

r_gp_conditional_cov <- function(L_K_obs, K_pred_obs, K_pred_pred, delta) {
  V <- forwardsolve(L_K_obs, t(K_pred_obs))
  K_pred_pred - t(V) %*% V + diag(delta, nrow(K_pred_pred))
}

r_multi_normal_cholesky_lcdf <- function(y, mu, L_Sigma) {
  z <- forwardsolve(L_Sigma, y - mu)
  sum(pnorm(z, log.p = TRUE))
}
