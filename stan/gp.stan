// /** Calculate Gaussian process variance-covariance matrix. 
// *
//   * @param x Proximity measures
// * @param alpha GP variance parameter
// * @param rho GP Smoothness/scale parameter
// * @return Variance-covariance matrix
// */
//   matrix calc_gp_vcov(array[] real x, real alpha, real rho) {
//     return gp_exp_quad_cov(x, alpha, rho);
//   }

matrix gp_exp_quad_cov(array[] real x, real alpha, real rho, real sigma) {
  return gp_exp_quad_cov(x, alpha, rho) + diag_matrix(rep_vector(sigma, size(x)));
}

matrix gp_matern32_cov(array[] real x, real alpha, real rho, real sigma) {
  return gp_matern32_cov(x, alpha, rho) + diag_matrix(rep_vector(sigma, size(x)));
}

matrix gp_matern52_cov(array[] real x, real alpha, real rho, real sigma) {
  return gp_matern52_cov(x, alpha, rho) + diag_matrix(rep_vector(sigma, size(x)));
}

/** Calculate Gaussian process Cholesky variance-covariance matrix. 
*
  * @param x Proximity measures
* @param alpha GP variance parameter
* @param rho GP Smoothness/scale parameter
* @param delta Variance or small epsilon to add to ensure proper matrix
* @return Variance-covariance matrix
*/
matrix gp_exp_quad_cholesky_cov(array[] real x, real alpha, real rho, real delta) {
  return cholesky_decompose(gp_exp_quad_cov(x, alpha, rho, delta));
}

matrix gp_matern32_cholesky_cov(array[] real x, real alpha, real rho, real delta) {
  return cholesky_decompose(gp_matern32_cov(x, alpha, rho, delta));
}

matrix gp_matern52_cholesky_cov(array[] real x, real alpha, real rho, real delta) {
  return cholesky_decompose(gp_matern52_cov(x, alpha, rho, delta));
}

/** Calculate one dimensional GP predictor.
*
  * @param x Proximity measures
* @param intercept GP mean
* @param alpha GP variance parameter
* @param rho GP Smoothness/scale parameter
* @param delta Variance or small epsilon to add to ensure proper matrix
* @param eta Standard normal (raw) parameters
* @return GP values for the given `x` 
*/
  vector calc_gp_pred(array[] real x, real intercept, real alpha, real rho, real delta, vector eta) {
    int n_x = size(x);
    matrix[n_x, n_x] L_K = gp_exp_quad_cholesky_cov(x, alpha, rho, delta); 
    
    return intercept + L_K * eta;
  }  

row_vector calc_gp_pred(array[] real x, real intercept, real alpha, real rho, real delta, row_vector eta) {
  int n_x = size(x);
  matrix[n_x, n_x] L_K = gp_exp_quad_cholesky_cov(x, alpha, rho, delta); 
  
  return intercept + eta * L_K';
}

tuple(vector, matrix) gp_conditional(
  vector mu_obs,             // Mean for the observed 
  vector mu_pred,            // Mean for the predicted 
  vector y_obs,              // Observations
  matrix K_obs_obs,          // Cov between observed points
  matrix K_pred_obs,         // Cross cov
  matrix K_pred_pred,        // Cov between prediction points
  real delta
) {
  int n_obs = size(mu_obs);
  int n_pred = size(mu_pred);
  
  matrix[n_obs, n_obs] L_K = cholesky_decompose(K_obs_obs);
  vector[n_obs] K_div_y_obs = mdivide_left_tri_low(L_K, y_obs - mu_obs); // inverse(tri(L_K)) * y

  K_div_y_obs = mdivide_right_tri_low(K_div_y_obs', L_K)'; // (inverse(tri(L_K)) * y)' * inverse(L_K))'

  matrix[n_obs, n_pred] v_pred = mdivide_left_tri_low(L_K, K_pred_obs'); // inverse(L_K) * K(X,X*)

// Just walking through these calculations to ensure it's doing the right thing.
  // N(mu_pred + K(X,X*)' * (inverse(tri(L_K)) * (y - mu_obs))' * inverse(L_K))', K(X*,X*) - (inverse(L_K) * K(X,X*))' * inverse(L_K) * K(X,X*))
  // N(mu_pred + K(X*,X) * inverse(L_K)' * inverse(L_K) * y, K(X*,X*) - K(X,X*)' * inverse(L_K)' * inverse(L_K) * K(X,X*))
// N(mu_pred + K(X*,X) * inverse(L_K'L_K) * y, K(X*,X*) - K(X*,X) * inverse(L_K'L_K) * K(X,X*))
// N(mu_pred + K(X*,X) * inverse(K(X,X)) * y, K(X*,X*) - K(X*,X) * inverse(K(X,X)) * K(X,X*)) <-- Correct!
  vector[n_pred] mu_cond = mu_pred + K_pred_obs * K_div_y_obs;
  matrix[n_pred, n_pred] K_pred_missing = K_pred_pred - v_pred' * v_pred + diag_matrix(rep_vector(delta, n_pred));
 
  return(mu_cond, K_pred_missing); 
}

/** This is the calculation needed to extrapolate a GP that is fit using observed y and x. We are predicting for x*.
 * For details, see Rasmussen' and Williams' "Gaussian Processes for Machine Learning".
 *
 * @param x_pred Proxmity measures to predict for
 * @param y Observed outcomes
 * @param x Observed proxmity measures
 * @param K_obs GP variance-covariance matrix for observed `(x, y)`
 * @param alpha GP variance parameter
 * @param rho GP Smoothness/scale parameter
 * @param delta Variance or small epsilon to add to ensure proper matrix
 * @return Predicted GP values conditional on observed data (interpolated from) 
 */
// vector multi_normal_rng(vector y_obs, array[] real x_obs, array[] real x_pred, matrix K_obs, matrix K_pred, real alpha, real rho, real delta) {
//   return gp_pred_rng(zeros_vector(size(x_obs)), zeros_vector(size(x_pred)), y_obs, x_obs, x_pred, K_obs, K_pred, alpha, rho, delta);
// }

// vector multi_normal_rng(real mu_obs, real mu_pred, vector y_obs, array[] real x_obs, array[] real x_pred, matrix K_obs, matrix K_pred, real alpha, real rho, real delta) {
vector multi_normal_rng(real mu_obs, real mu_pred, vector y_obs, matrix K_obs, matrix K_pred_obs, matrix K_pred) { //, real alpha, real rho, real delta) {
  return multi_normal_rng(rep_vector(mu_obs, rows(y_obs)), rep_vector(mu_pred, rows(K_pred)), y_obs, K_obs, K_pred_obs, K_pred);
}

// vector multi_normal_rng(real mu, vector y_obs, array[] real x_obs, array[] real x_pred, matrix K_obs, matrix K_pred, real alpha, real rho, real delta) {
vector multi_normal_rng(real mu, vector y_obs, matrix K_obs, matrix K_pred_obs, matrix K_pred) { //, real alpha, real rho, real delta) {
  return multi_normal_rng(mu, mu, y_obs, K_obs, K_pred_obs, K_pred);
}
 
vector multi_normal_rng(vector mu_obs, vector mu_pred, vector y_obs, matrix K_obs, matrix K_pred_obs, matrix K_pred) { // , real alpha, real rho, real delta) {
  int n_pred = rows(K_pred);
  
  vector[n_pred] mu_cond;
  matrix[n_pred, n_pred] Sigma_cond;
  
  (mu_cond, Sigma_cond) = gp_conditional(mu_obs, mu_pred, y_obs, K_obs, K_pred_obs, K_pred, 1e-5);
 
  // return multi_normal_cholesky_rng(mu_cond, cholesky_decompose(Sigma_cond));
  return multi_normal_rng(mu_cond, Sigma_cond);
}

// vector gp_pred_rng(array[] real x_pred, vector y, array[] real x, matrix K_obs, matrix K_missing, real alpha, real rho) {
//   return gp_pred_rng(x_pred, y, x, K_obs, K_missing, alpha, rho, 0);
// }
// 
// vector gp_pred_rng(array[] real x_pred, vector y, array[] real x, matrix K_obs, real alpha, real rho, real delta) {
//   return gp_pred_rng(x_pred, y, x, gp_exp_quad_cov(x_pred, alpha, rho), alpha, rho, delta);
// }
// 
// vector gp_pred_rng(array[] real x_pred, vector y, array[] real x, matrix K_obs, real alpha, real rho) {
//   return gp_pred_rng(x_pred, y, x, K_obs, alpha, rho, 0);
// }

// vector calc_gp_pred(array[] real x, real intercept, real alpha, real rho, vector eta) {
//   return calc_gp_pred(x, intercept, alpha, 1e-9, eta);
// }  

/**
 * Calculate log of multivariate cholesky normal CDF 
 * 
 * @param y Vector at which to evaluate the CDF
 * @param mu Mean vector
 * @param Sigma Covariance matrix
 * @return Log of multivariate normal CDF evaluated at y
 */
real multi_normal_cholesky_lcdf(vector y, vector mu, matrix L_Sigma) {
  int K = rows(y);
  // These are now independent standard normal random variables
  vector[K] z = mdivide_left_tri_low(L_Sigma, y - mu);
  
  return std_normal_lcdf(z);
}

real multi_normal_cholesky_lcdf(vector y, real mu, matrix L_Sigma) {
  return multi_normal_cholesky_lcdf(y | rep_vector(mu, size(y)), L_Sigma);
}

real multi_normal_lcdf(vector y, vector mu_obs, vector mu_pred, vector y_cond, matrix K_obs, matrix K_pred_obs, matrix K_pred, real delta) {
  int n_obs = rows(y);
  int n_pred = rows(mu_pred);
  
  vector[n_pred] mu_cond;
  matrix[n_pred, n_pred] Sigma_cond;
  (mu_cond, Sigma_cond) = gp_conditional(mu_obs, mu_pred, y_cond, K_obs, K_pred_obs, K_pred, delta);
  
  matrix[n_pred, n_pred] L_Sigma_cond = cholesky_decompose(Sigma_cond);
  
  return multi_normal_cholesky_lcdf(y | mu_cond, L_Sigma_cond);
}