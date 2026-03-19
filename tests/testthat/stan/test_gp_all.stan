functions {
  #include "gp.stanfunctions"
}
data {
  int<lower=1> N;
  array[N] real x;
  real alpha;
  real rho;
  real sigma;
  real delta;
  // conditional GP
  int<lower=1> N_obs;
  int<lower=1> N_pred;
  array[N_obs] real x_obs;
  array[N_pred] real x_pred;
  vector[N_obs] y_obs;
  // scale_process_sd
  int<lower=1> N_tp;
  array[N_tp] real time_points;
  vector[2] process_sd;
  // lcdf
  int<lower=1> D;
  vector[D] lcdf_y;
}
generated quantities {
  // Kernels (with nugget sigma)
  matrix[N, N] K_eq   = gp_exp_quad_cov(x, alpha, rho, sigma);
  matrix[N, N] K_m32  = gp_matern32_cov(x, alpha, rho, sigma);
  matrix[N, N] K_m52  = gp_matern52_cov(x, alpha, rho, sigma);
  // Cholesky factors (with nugget delta)
  matrix[N, N] L_eq   = gp_exp_quad_cholesky_cov(x, alpha, rho, delta);
  matrix[N, N] L_m32  = gp_matern32_cholesky_cov(x, alpha, rho, delta);
  matrix[N, N] L_m52  = gp_matern52_cholesky_cov(x, alpha, rho, delta);
  // Verify L L' == K (for exp_quad with nugget delta)
  matrix[N, N] K_from_L = L_eq * L_eq';
  // calc_gp_pred (vector eta, scalar intercept = 0.0)
  vector[N] eta = rep_vector(0.5, N);
  vector[N] gp_out = calc_gp_pred(x, 0.0, alpha, rho, delta, eta);
  // scale_process_sd
  matrix[N_tp, 2] scaled_sd = scale_process_sd(time_points, process_sd);
  // GP conditional mean + cov (zero-mean overloads)
  matrix[N_obs, N_obs] K_oo   = gp_exp_quad_cov(x_obs, alpha, rho, delta);
  // Stan built-in 4-arg cross-covariance: gp_exp_quad_cov(x1, x2, alpha, rho) -> matrix[N_pred, N_obs]
  matrix[N_pred, N_obs] K_po  = gp_exp_quad_cov(x_pred, x_obs, alpha, rho);
  matrix[N_pred, N_pred] K_pp = gp_exp_quad_cov(x_pred, alpha, rho, delta);
  matrix[N_obs, N_obs] L_oo   = cholesky_decompose(K_oo);
  vector[N_pred] cond_mean    = gp_conditional_mean(y_obs, L_oo, K_po);
  matrix[N_pred, N_pred] cond_cov = gp_conditional_cov(L_oo, K_po, K_pp, delta);
  // multi_normal_cholesky_lcdf with identity covariance and scalar mu = 0.0
  // Stan requires pipe syntax for user-defined functions ending in _lcdf
  matrix[D, D] L_id  = diag_matrix(rep_vector(1.0, D));
  real lcdf_val = multi_normal_cholesky_lcdf(lcdf_y | 0.0, L_id);
}
