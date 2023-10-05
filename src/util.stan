vector calc_gp_pred(array[] real x, real intercept, real alpha, real rho, real delta, vector eta) {
  int n_x = size(x);
  matrix[n_x, n_x] K = gp_exp_quad_cov(x, alpha, rho) + diag_matrix(rep_vector(delta, n_x));
  matrix[n_x, n_x] L_K = cholesky_decompose(K);
  
  return intercept + L_K * eta;
}  