// Calculate the last observed measure for each patient. 
array[] int get_max_t(array[] int t_measure, array[] int n_measures, array[] int n_patient_tumors) {
  int n_patients = size(n_patient_tumors);
  int t_pos = 1;
  int tumor_pos = 1;
  array[n_patients] int max_t;
    
  for (i in 1:n_patients) {
    int tumor_end = tumor_pos + n_patient_tumors[i] - 1;
    int t_end = t_pos + sum(n_measures[tumor_pos:tumor_end]) - 1; 
    
    max_t[i] = max(t_measure[t_pos:t_end]);
    
    t_pos = t_end + 1;
    tumor_pos = tumor_end + 1;
  }
  
  return max_t;
}

matrix calc_gp_vcov(array[] real x, real alpha, real rho, real delta) {
  int n_x = size(x);
  return gp_exp_quad_cov(x, alpha, rho) + diag_matrix(rep_vector(delta, n_x));
}

matrix calc_gp_cholesky_vcov(array[] real x, real alpha, real rho, real delta) {
  return cholesky_decompose(calc_gp_vcov(x, alpha, rho, delta));
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


vector gp_pred_rng(array[] real x_pred, vector y, array[] real x, matrix K_obs, real alpha, real rho, real delta) {
  int n_obs = rows(y);
  int n_pred = size(x_pred);
  
  matrix[n_obs, n_obs] L_K = cholesky_decompose(K_obs);
  vector[n_obs] K_div_y_obs = mdivide_left_tri_low(L_K, y); // inverse(tri(L_K)) * y
  
  K_div_y_obs = mdivide_right_tri_low(K_div_y_obs', L_K)'; // (inverse(tri(L_K)) * y)' * inverse(L_K))'
  
  matrix[n_obs, n_pred] K_x_obs_x_pred = gp_exp_quad_cov(x, x_pred, alpha, rho); // K(X,X*)
  matrix[n_obs, n_pred] v_pred = mdivide_left_tri_low(L_K, K_x_obs_x_pred); // inverse(L_K) * K(X,X*)
  
  // Just walking through these calculations to ensure it's doing the right thing. 
  // N(K(X,X*)' * (inverse(tri(L_K)) * y)' * inverse(L_K))', K(X*,X*) - (inverse(L_K) * K(X,X*))' * inverse(L_K) * K(X,X*))
  // N(K(X*,X) * inverse(L_K)' * inverse(L_K) * y, K(X*,X*) - K(X,X*)' * inverse(L_K)' * inverse(L_K) * K(X,X*))
  // N(K(X*,X) * inverse(L_K'L_K) * y, K(X*,X*) - K(X*,X) * inverse(L_K'L_K) * K(X,X*))
  // N(K(X*,X) * inverse(K(X,X)) * y, K(X*,X*) - K(X*,X) * inverse(K(X,X)) * K(X,X*)) <-- Correct!
  return multi_normal_rng(K_x_obs_x_pred' * K_div_y_obs, gp_exp_quad_cov(x_pred, alpha, rho) - v_pred' * v_pred + diag_matrix(rep_vector(delta, n_pred)));
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

// A missing measure is defined as one that lies between a _tumor's_ first assessment to the _patient's_ last assessment.
array[] int calculate_n_missing_measures(array[] int n_measures, array[] int t_measure, array[] int n_patient_tumors) {
  int tumor_pos = 1;
  int t_measure_pos = 1;
  int n_patients = size(n_patient_tumors);
  int n_tumors = size(n_measures);
  array[n_tumors] int n_missing_measures = rep_array(0, n_tumors);
  
  for (i in 1:n_patients) {
    int tumor_end = tumor_pos + n_patient_tumors[i] - 1;
    int first_tumor_t_measure_pos = t_measure_pos;
    int last_tumor_t_measure_end = t_measure_pos + sum(n_measures[tumor_pos:tumor_end]) - 1;
    
    int max_patient_t = max(t_measure[first_tumor_t_measure_pos:last_tumor_t_measure_end]);
    // int max_patient_t = max(t_measure[t_measure_pos:t_measure_end]);
    // int full_patient_measure_width = max_patient_t - min_patient_t + 1;
    
    for (j in 1:n_patient_tumors[i]) {
      int t_measure_end = t_measure_pos + n_measures[tumor_pos] - 1;
      
      int min_tumor_t = min(t_measure[t_measure_pos:t_measure_end]);
      int full_patient_measure_width = max_patient_t - min_tumor_t + 1;
      
      n_missing_measures[tumor_pos] = full_patient_measure_width - n_measures[tumor_pos];
      
      tumor_pos += 1;  
      t_measure_pos = t_measure_end + 1;
    }
    
  }
  
  return n_missing_measures;
}

array[] int calculate_t_missing_measure(
  array[] int n_measures, array[] int n_missing_measures, array[] int t_measure, array[] int n_patient_tumors 
) { 
  int n_patients = size(n_patient_tumors); 
  int tumor_pos = 1;
  int t_measure_pos = 1;
  int t_missing_measure_pos = 1;
  int n_tumors = size(n_measures);
  array[sum(n_missing_measures)] int t_missing_measure;
  
  for (i in 1:n_patients) {
    int tumor_end = tumor_pos + n_patient_tumors[i] - 1;
    int first_tumor_t_measure_pos = t_measure_pos;
    int last_tumor_t_measure_end = t_measure_pos + sum(n_measures[tumor_pos:tumor_end]) - 1;
    
    // int t_measure_end = t_measure_pos + sum(n_measures[tumor_pos:tumor_end]) - 1;
    // int min_patient_t = min(t_measure[t_measure_pos:t_measure_end]);
    int max_patient_t = max(t_measure[first_tumor_t_measure_pos:last_tumor_t_measure_end]);
    
    for (j in 1:n_patient_tumors[i]) {
      int t_measure_end = t_measure_pos + n_measures[tumor_pos] - 1;
      
      int min_tumor_t = min(t_measure[t_measure_pos:t_measure_end]);
      int measures_checked = 0;
      
      for (k in min_tumor_t:max_patient_t) {
        if (measures_checked >= n_measures[tumor_pos] || t_measure[t_measure_pos] > k) {
          if (t_missing_measure_pos > sum(n_missing_measures)) {
            print("i = ", i, ", j = ", j);
          }
          
          t_missing_measure[t_missing_measure_pos] = k;
          t_missing_measure_pos += 1;
        } else {
          t_measure_pos += 1;
          measures_checked += 1;
        }
      }
      
      tumor_pos += 1;
    }
  }
  
  return t_missing_measure;
}

array[] int calc_n_screening_t(array[] int n_patient_tumors, array[] int n_measures, array[] int t_measure) {
  int n_patients = size(n_patient_tumors);
  array[sum(n_patient_tumors)] int n_screening_t = rep_array(0, sum(n_patient_tumors));
  
  int tumor_pos = 1;
  int t_measure_pos = 1;
  
  for (i in 1:n_patients) {
    int tumor_end = tumor_pos + n_patient_tumors[i] - 1;
    
    for (j in 1:n_patient_tumors[i]) {
      int t_measure_end = t_measure_pos + n_measures[tumor_pos + j - 1] - 1;
      
      for (tp in t_measure_pos:t_measure_end) {
        if (t_measure[tp] <= 0) {
          n_screening_t[tumor_pos + j - 1] += 1;
        }
      }
      
      t_measure_pos = t_measure_end + 1;
    }
    
    tumor_pos = tumor_end + 1;
  }
  
  return n_screening_t;
}