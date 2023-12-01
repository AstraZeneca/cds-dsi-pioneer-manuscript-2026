matrix calc_gp_cholesky_vcov(array[] real x, real alpha, real rho, real delta) {
  int n_x = size(x);
  matrix[n_x, n_x] K = gp_exp_quad_cov(x, alpha, rho) + diag_matrix(rep_vector(delta, n_x));
  return cholesky_decompose(K);
}

// matrix calc_gp_cholesky_vcov(array[] real x, real alpha, real rho) {
//   return calc_gp_cholesky_vcov(x, alpha, rho, 1e-9); 
// }

// Calculate one dimensional GP predictor
vector calc_gp_pred(array[] real x, real intercept, real alpha, real rho, real delta, vector eta) {
  int n_x = size(x);
  matrix[n_x, n_x] L_K = calc_gp_cholesky_vcov(x, alpha, rho, delta); 
  
  return intercept + L_K * eta;
}  

// vector calc_gp_pred(array[] real x, real intercept, real alpha, real rho, vector eta) {
//   return calc_gp_pred(x, intercept, alpha, 1e-9, eta);
// }  

// Scale tumor sizes by the standard deviation of all non-zero tumors (a size of zero means the tumor doesn't exist yet/anymore).
tuple(real, real, vector) standardize_nonzero_tumor_sizes(vector tumor_size) {
  int n_tumor_measures = rows(tumor_size);
  array[n_tumor_measures] int nonzero_tumor_idx;
  int measured_pos = 1;
  real tumor_mean;
  real tumor_sd;
  
  for (t in 1:n_tumor_measures) {
    if (tumor_size[t] > 0) {
      nonzero_tumor_idx[measured_pos] = t;
      measured_pos += 1;
    }
  }
  
  tumor_mean = mean(tumor_size[nonzero_tumor_idx[:(measured_pos - 1)]]);
  tumor_sd = sd(tumor_size[nonzero_tumor_idx[:(measured_pos - 1)]]);
  
  return (tumor_mean, tumor_sd, (tumor_size - tumor_mean) / tumor_sd); 
}

array[] real pfs_quantiles_from_prob(vector exit_prob, array[] real p) {
  int max_t = rows(exit_prob);
  int n_p = size(p);
  array[n_p] int sorted_p_idx = sort_indices_asc(p);
  vector[max_t + 1] cumul_prob = append_row(0.0, cumulative_sum(exit_prob)); 
  array[n_p] real q;
  int pfs = 1;
  
  for (p_index in 1:n_p) {
    real curr_p = p[sorted_p_idx][p_index];
    real q_part;
    
    while (cumul_prob[pfs + 1] < curr_p) {
      pfs += 1;
    }
    
    q_part = (curr_p - cumul_prob[pfs]) / (cumul_prob[pfs + 1] - cumul_prob[pfs]);
    
    q[sorted_p_idx[p_index]] = (pfs - 1) * (1 - q_part) + pfs * q_part; 
  }
  
  return q;
}