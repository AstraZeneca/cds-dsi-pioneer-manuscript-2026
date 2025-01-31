/** Calculate the last observed measure for each patient. 
 *
 * @param t_measure The week each assessment was done.
 * @param n_measures The number of assessments per tumor.
 * @param n_patient_tumors Array with the number of tumors per patient.
 * @return Get the week of last observation per patient.
 */
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

/** Calculate Gaussian process variance-covariance matrix. 
 *
 * @param x Proximity measures
 * @param alpha GP variance parameter
 * @param rho GP Smoothness/scale parameter
 * @param delta Small epsilon to add to ensure proper matrix
 * @return Variance-covariance matrix
 */
matrix calc_gp_vcov(array[] real x, real alpha, real rho, real delta) {
  int n_x = size(x);
  return gp_exp_quad_cov(x, alpha, rho) + diag_matrix(rep_vector(delta, n_x));
}

/** Calculate Gaussian process Cholesky variance-covariance matrix. 
 *
 * @param x Proximity measures
 * @param alpha GP variance parameter
 * @param rho GP Smoothness/scale parameter
 * @param delta Small epsilon to add to ensure proper matrix
 * @return Variance-covariance matrix
 */
matrix calc_gp_cholesky_vcov(array[] real x, real alpha, real rho, real delta) {
  return cholesky_decompose(calc_gp_vcov(x, alpha, rho, delta));
}

/** Calculate one dimensional GP predictor.
 *
 * @param x Proximity measures
 * @param intercept GP mean
 * @param alpha GP variance parameter
 * @param rho GP Smoothness/scale parameter
 * @param delta Small epsilon to add to ensure proper matrix
 * @param eta Standard normal (raw) parameters
 * @return GP values for the given `x` 
 */
vector calc_gp_pred(array[] real x, real intercept, real alpha, real rho, real delta, vector eta) {
  int n_x = size(x);
  matrix[n_x, n_x] L_K = calc_gp_cholesky_vcov(x, alpha, rho, delta); 
  
  return intercept + L_K * eta;
}  

row_vector calc_gp_pred(array[] real x, real intercept, real alpha, real rho, real delta, row_vector eta) {
  int n_x = size(x);
  matrix[n_x, n_x] L_K = calc_gp_cholesky_vcov(x, alpha, rho, delta); 
  
  return intercept + eta * L_K;
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
 * @param delta Small epsilon to add to ensure proper matrix
 * @return Predicted GP values conditional on observed data (interpolated from) 
 */
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

/** Missing measure is defined as one that lies between a _tumor's_ first assessment to the _patient's_ last assessment. Basically,
 * we're counting how many intervals (weeks) we don't have observed assessments of tumor size, for each tumor.
 *
 * @param n_measures The number of assessments per tumor.
 * @param t_measure The week each assessment was done.
 * @param n_patient_tumors Array with the number of tumors per patient.
 * @return Number of missing assessments per tumor
 */
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

/** Return the actual t for which we don't have observed tumor size assessments.
 *
 * @param n_measures The number of assessments per tumor.
 * @param n_missing_measures Number of missing assessments per tumor
 * @param t_measure The week each assessment was done.
 * @param n_patient_tumors Array with the number of tumors per patient.
 * @return The actual weeks in which assessments are unobserved
 */
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
    
    int max_patient_t = max(t_measure[first_tumor_t_measure_pos:last_tumor_t_measure_end]);
    
    for (j in 1:n_patient_tumors[i]) {
      int t_measure_end = t_measure_pos + n_measures[tumor_pos] - 1;
      
      int min_tumor_t = min(t_measure[t_measure_pos:t_measure_end]);
      int measures_checked = 0;
      
      for (k in min_tumor_t:max_patient_t) {
        if (measures_checked >= n_measures[tumor_pos] || t_measure[t_measure_pos] > k) {
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

/** Get number of values in x that are less than or equal to y
 */
int num_leq(array[] int x, int y) {
  int n = 0;
  array[size(x)] int sorted_x = sort_asc(x);
  
  for (i in 1:size(x)) {
    if (sorted_x[i] <= y) {
      n += 1;
    } else {
      break;
    }
  }
  
  return n;
}

/** Identify which elements in a binary array are 0 and which are 1.
 * @param mask Binary array
 * @return tuple(indices of 0 elements, indices of 1 elements)
 */
tuple(array[] int, array[] int) get_mask_idx(array[] int mask) {
  int n = size(mask);
  int n_0 = n - sum(mask);
  array[n] int sorted_idx = sort_indices_asc(mask);

  return(sorted_idx[:n_0], sorted_idx[(n_0 + 1):]); 
}

/** Repeat each value a specific number of times.
 */
array[] int rep_each(array[] int to_repeat, int repeats) {
  int n = size(to_repeat);
  array[n * repeats] int repeated;
  
  int pos = 1;
  
  for (i in 1:n) {
    int end = pos + repeats - 1;
    repeated[pos:end] = rep_array(to_repeat[i], repeats);
    pos = end + 1;
  }
  
  return(repeated);
}

real months_to_weeks(int mon) {
  return mon * 365.25 / (7 * 12);
}

int calendar_date_to_study_date(int first_calendar_date, int calendar_date) {
  return calendar_date - first_calendar_date + 1;
}

array[] int calendar_date_to_study_date(array[] int first_calendar_date, array[] int calendar_date) {
  int n = size(first_calendar_date);
  array[n] int study_date;
  
  for (i in 1:n) {
    study_date[i] = calendar_date_to_study_date(first_calendar_date[i], calendar_date[i]);
  }
  
  return study_date;
}

array[] int calendar_date_to_study_date(array[] int first_calendar_date, int calendar_date) {
  int n = size(first_calendar_date);
  array[n] int study_date;
  
  for (i in 1:n) {
    study_date[i] = calendar_date_to_study_date(first_calendar_date[i], calendar_date);
  }
  
  return study_date;
}

int study_date_to_calendar_date(int first_calendar_date, int study_date) {
  return first_calendar_date + study_date - 1;
}

array[] int study_date_to_calendar_date(array[] int first_calendar_date, array[] int study_date) {
  int n = size(first_calendar_date);
  array[n] int calendar_date;
  
  for (i in 1:n) {
    calendar_date[i] = study_date_to_calendar_date(first_calendar_date[i], study_date[i]);
  }
  
  return calendar_date;
}