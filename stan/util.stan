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
 * @return Variance-covariance matrix
 */
matrix calc_gp_vcov(array[] real x, real alpha, real rho) {
  return gp_exp_quad_cov(x, alpha, rho);
}

matrix calc_gp_vcov(array[] real x, real alpha, real rho, real sigma) {
  return calc_gp_vcov(x, alpha, rho) + diag_matrix(rep_vector(sigma, size(x)));
}

/** Calculate Gaussian process Cholesky variance-covariance matrix. 
 *
 * @param x Proximity measures
 * @param alpha GP variance parameter
 * @param rho GP Smoothness/scale parameter
 * @param delta Variance or small epsilon to add to ensure proper matrix
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
 * @param delta Variance or small epsilon to add to ensure proper matrix
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
  
  return intercept + eta * L_K';
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
vector gp_pred_rng(array[] real x_pred, vector y, array[] real x, matrix K_obs, matrix K_missing, real alpha, real rho, real delta) {
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
  matrix[n_pred, n_pred] K_pred_missing = K_missing - v_pred' * v_pred + diag_matrix(rep_vector(delta, n_pred));
 
  return multi_normal_cholesky_rng(K_x_obs_x_pred' * K_div_y_obs, cholesky_decompose(K_pred_missing));
}

vector gp_pred_rng(array[] real x_pred, vector y, array[] real x, matrix K_obs, matrix K_missing, real alpha, real rho) {
  return gp_pred_rng(x_pred, y, x, K_obs, K_missing, alpha, rho, 0);
}

vector gp_pred_rng(array[] real x_pred, vector y, array[] real x, matrix K_obs, real alpha, real rho, real delta) {
  return gp_pred_rng(x_pred, y, x, gp_exp_quad_cov(x_pred, alpha, rho), alpha, rho, delta);
}

vector gp_pred_rng(array[] real x_pred, vector y, array[] real x, matrix K_obs, real alpha, real rho) {
  return gp_pred_rng(x_pred, y, x, K_obs, alpha, rho, 0);
}

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
  array[n_tumors] int n_missing_measures = zeros_int_array(n_tumors);
  
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

tuple(array[] int, array[] int) calculate_n_missing_visits(array[] int visit_pos, array[] int visits) {
  return calculate_n_missing_visits(visit_pos, visits, 0, 0);
}

tuple(array[] int, array[] int) calculate_n_missing_visits(array[] int visit_pos, array[] int visits, int make_unique, int post_treatment) {
  int n = size(visit_pos) - 1;
  array[n] int n_missing = zeros_int_array(n);
  
  for (i in 1:n) {
    int curr_visit_pos, curr_visit_end;
    (curr_visit_pos, curr_visit_end) = get_pos(visit_pos, i);
   
    if (post_treatment) { 
      while (curr_visit_pos <= curr_visit_end && visits[curr_visit_pos] <= 0) {
        curr_visit_pos += 1;
      }
    }
    
    int max_t_width = max(visits[curr_visit_pos:curr_visit_end]) - min(visits[curr_visit_pos:curr_visit_end]) + 1;
    int n_visits = make_unique ? num_unique(visits[curr_visit_pos:curr_visit_end]) : (curr_visit_end - curr_visit_pos + 1);
    
    n_missing[i] = max_t_width - n_visits;
  }
  
  return (n_missing, create_pos(n_missing));
}

tuple(array[] int, array[] int) get_missing_visits(array[] int visit_pos, array[] int visits, array[] int missing_visit_pos) {
  return get_missing_visits(visit_pos, visits, missing_visit_pos, 0);
}

tuple(array[] int, array[] int) get_missing_visits(array[] int visit_pos, array[] int visits, array[] int missing_visit_pos, int post_treatment) {
  int n = size(visit_pos) - 1;
  array[n] int missing_size = get_pos_size(missing_visit_pos);
  int n_total_missing = sum(missing_size);
  array[n_total_missing] int missing_visits, missing_visits_idx;
  
  
  for (i in 1:n) {
    int missing_pos, missing_end;
    (missing_pos, missing_end) = get_pos(missing_visit_pos, i);
   
    int n_curr_visits = get_pos_size(visit_pos, i), visit_offset = 0; 
    array[n_curr_visits] int sorted_visits = sort_asc(get_int_sub_array(visits, visit_pos, i)); 
    
    if (post_treatment) {
      while (visit_offset < n_curr_visits && sorted_visits[visit_offset + 1] <= 0) {
        visit_offset += 1;
      }
    }
    
    int min_visit = sorted_visits[1 + visit_offset], max_visit = sorted_visits[n_curr_visits];
    
    int visit_count = 0;
    int next_visit = min_visit;
    
    for (m in missing_pos:missing_end) {
      while (next_visit < max_visit && sorted_visits[visit_count + 1 + visit_offset] == next_visit) {
        visit_count += 1;
        next_visit += 1;
      }
      
      missing_visits[m] = next_visit;
      missing_visits_idx[m] = next_visit - min_visit + 1;
      next_visit += 1;
    }
  }
  
  return (missing_visits, missing_visits_idx);
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

int num_unique(array[] int x) {
  return num_unique(x, 1);
}

int num_unique(array[] int x, int post_treatment) {
  int n = size(x);
  int count = 0, last = min(x) - 1;
  array[n] int sorted_x = sort_asc(x);
  
  for (i in 1:n) {
    if ((!post_treatment || sorted_x[i] > 0) && sorted_x[i] > last) {
      count += 1;
      last = sorted_x[i];
    }
  }
  
  return count;
}


array[] int num_unique(array[] int x, array[] int pos) {
  return num_unique(x, pos, 1);
}

array[] int num_unique(array[] int x, array[] int pos, int post_treatment) {
  int n = size(pos) - 1;
  array[n] int count = zeros_int_array(n);
  
  for (i in 1:n) {
    count[i] = num_unique(get_int_sub_array(x, pos, i), post_treatment);
  }
  
  return count;
}

array[] int num_unique(array[] int x, array[] int pos, array[] int sub_pos) {
  return num_unique(x, pos, sub_pos, 1);
}
  
array[] int num_unique(array[] int x, array[] int pos, array[] int sub_pos, int post_treatment) {
  int n = size(pos) - 1;
  array[n] int count = zeros_int_array(n);
  
  for (i in 1:n) {
    int i_pos, i_end;
    (i_pos, i_end) = get_pos(pos, i);
    
    count[i] = num_unique(get_int_sub_array(x, sub_pos, i_pos, i_end), post_treatment);
  }
  
  return count;
}

array[] int unique(array[] int x) {
  return unique(x, 1);
}
  
array[] int unique(array[] int x, int post_treatment) {
  int n = size(x);
  int count = 0, last = min(x) - 1;
  array[n] int sorted_x = sort_asc(x);
  int n_unique = num_unique(x, post_treatment);
  array[n_unique] int unique_x;
  
  for (i in 1:n) {
    
    if ((!post_treatment || sorted_x[i] > 0) && sorted_x[i] > last) {
      count += 1;
      last = sorted_x[i];
      unique_x[count] = last;
    }
  }
  
  return unique_x;
}

// If I call this unique the compiler complains about ambiguity which doesn't make sense ¯\_(ツ)_/¯
tuple(array[] int, array[] int) unique_by_pos(array[] int x, array[] int pos) {
  return unique_by_pos(x, pos, 1);
}
  
tuple(array[] int, array[] int) unique_by_pos(array[] int x, array[] int pos, int post_treatment) {
  int n = size(pos) - 1;
  array[n] int n_unique_x = num_unique(x, pos, post_treatment);
  array[n + 1] int unique_pos = create_pos(n_unique_x);
  array[sum(n_unique_x)] int unique_x;
  
  for (i in 1:n) {
    int n_i = get_pos_size(pos, i);
    array[n_i] int sorted_x_i = sort_asc(get_int_sub_array(x, pos, i));
    int count = 0, last = sorted_x_i[1] - 1;
    
    for (j in 1:n_i) {
      if ((!post_treatment || sorted_x_i[j] > 0) && sorted_x_i[j] > last) {
        last = sorted_x_i[j];
        unique_x[unique_pos[i] + count] = last; 
        count += 1;
      }
    }
  }
  
  return (unique_x, unique_pos);
}

tuple(array[] int, array[] int) unique_by_pos(array[] int x, array[] int pos, array[] int sub_pos) {
  return unique_by_pos(x, pos, sub_pos, 1);
}

tuple(array[] int, array[] int) unique_by_pos(array[] int x, array[] int pos, array[] int sub_pos, int post_treatment) {
  int n = size(pos) - 1;
  array[n] int n_unique_x = num_unique(x, pos, sub_pos, post_treatment);
  array[n + 1] int unique_pos = create_pos(n_unique_x);
  array[sum(n_unique_x)] int unique_x;
  
  for (i in 1:n) {
    int i_pos, i_end;
    (i_pos, i_end) = get_pos(pos, i);
    
    
    int i_unique_pos, i_unique_end;
    (i_unique_pos, i_unique_end) = get_pos(unique_pos, i);
    
    unique_x[i_unique_pos:i_unique_end] = unique(get_int_sub_array(x, sub_pos, i_pos, i_end), post_treatment);
  }
  
  return (unique_x, unique_pos);
}

array[] int id2idx(array[] int id) {
  return id2idx(id, min(id));
}

array[] int id2idx(array[] int id, array[] int pos) {
  int n = size(pos) - 1;
  array[n] int n_p = get_pos_size(pos);
  array[sum(n_p)] int idx;
  
  for (p in 1:n) {
    if (get_pos_size(pos, p) > 0) {
      int p_pos, p_end;
      (p_pos, p_end) = get_pos(pos, p);
    
      idx[p_pos:p_end] = id2idx(id[p_pos:p_end]);
    } 
  }
 
  return idx;
}

array[] int id2idx(array[] int id, int min_id) {
  int n = size(id);
  array[n] int idx;
  
  for (i in 1:n) {
    idx[i] = id[i] - min_id + 1;
  }
  
  return idx;
}

array[] int get_level2level_idx(array[] int hi_level, array[] int low_level) {
  int size_hi = size(hi_level), size_low = size(low_level);
 
  array[size_low] int idx = zeros_int_array(size_low);
  array[size_low] int sorted_low_level = sort_asc(low_level);
  int curr_low_idx = 1;
  
  // print("hi_level = ", hi_level);
  // print("low_level = ", low_level);
  
  for (h in 1:size_hi) {
    // print("hi_level[h] = ", hi_level[h], ", sorted_low_level[curr_low_idx] = ", sorted_low_level[curr_low_idx]);
    
    if (hi_level[h] == sorted_low_level[curr_low_idx]) {
      idx[curr_low_idx] = h;
      
      if (curr_low_idx == size_low) {
        break;
      }
      
      curr_low_idx += 1;
    }
  }
  
  return idx;
}

array[] int get_level2level_idx(array[] int hi_level, array[] int low_level, array[] int low_pos) {
  int n_low = size(low_pos) - 1;
  int size_low = size(low_level);
  array[size_low] int idx = zeros_int_array(size_low);
  
  // print("hi_level = ", hi_level);
  // print("low_level = ", low_level);
  
  for (l in 1:n_low) {
    if (get_pos_size(low_pos, l) > 0) {
      int pos, end;
      (pos, end) = get_pos(low_pos, l);
      
      idx[pos:end] = get_level2level_idx(hi_level, low_level[pos:end]);
      
      // print(l, ": pos = ", get_pos(low_pos, l), ", idx[pos:end] = ", idx[pos:end]);
      // print(l, ": id[pos:end] = ", low_level[pos:end]);
      // print(l, ": hi_level = ", hi_level[idx[pos:end]]);
    } else {
      // print(l, ": pos = ", get_pos(low_pos, l));
    }
  }
  
  return idx;
}

array[] int get_level2level_idx(array[] int hi_level, array[] int hi_pos, array[] int low_level, array[] int low_pos, array[] int low_hi_pos) {
  int n_low = size(low_pos) - 1, n_hi = size(hi_pos) - 1;
  int size_low = size(low_level);
  array[size_low] int idx = zeros_int_array(size_low);
 
  for (h in 1:n_hi) {
    int low_id_from, low_id_to;
    (low_id_from, low_id_to) = get_pos(low_hi_pos, h);
    // print(h, ": low_id_from = ", low_id_from, ", low_id_to = ", low_id_to);
    int low_idx_start, low_idx_end;
    (low_idx_start, low_idx_end) = get_pos(low_pos, low_id_from, low_id_to);
    // print(h, ": low_idx_start = ", low_idx_start , ", low_idx_end = ", low_idx_end);
    
    // print(h, ": get_int_sub_array(hi_level, hi_pos, h) = ", get_int_sub_array(hi_level, hi_pos, h));
    // print(h, ": get_int_sub_array(low_level, low_pos, low_id_from, low_id_to) = ", get_int_sub_array(low_level, low_pos, low_id_from, low_id_to));
    // print(h, ": create_pos(low_pos, low_id_from, low_id_to) = ", create_pos(low_pos, low_id_from, low_id_to));
    
    idx[low_idx_start:low_idx_end] = get_level2level_idx(
      get_int_sub_array(hi_level, hi_pos, h), get_int_sub_array(low_level, low_pos, low_id_from, low_id_to), create_pos(low_pos, low_id_from, low_id_to)
    );
    
    // for (i in 1:(low_id_to - low_id_from + 1)) {
    //   print(h, ", ", i, ": ", get_int_sub_array(get_int_sub_array(low_level, low_pos, low_id_from, low_id_to), create_pos(low_pos, low_id_from, low_id_to), i)); 
    // }
    
    // print(h, ": idx[low_idx_start:low_idx_end] = ", idx[low_idx_start:low_idx_end]);
  }
  
  return idx;
}

array[] int get_idx_dict(array[] int idx) {
  int max_idx = size(idx) > 0 ? max(idx) : 0; 
  array[max_idx] int idx_dict = zeros_int_array(max_idx);
  
  int idx_pos = 1;
  
  for (i in 1:max_idx) {
    if (idx[idx_pos] == i) {
      idx_dict[i] = idx_pos;
      idx_pos += 1;
    }
  }
  
  return idx_dict;
}

tuple(array[] int, array[] int) get_idx_dict(array[] int idx, array[] int pos) {
  int n = size(pos) - 1;
  array[n] int max_idx = get_max(idx, pos);
  array[n + 1] int dict_pos = create_pos(max_idx);
  array[sum(max_idx)] int idx_dict;

  for (p in 1:n) {
      int p_pos, p_end;
      (p_pos, p_end) = get_pos(dict_pos, p);
    if (max_idx[p] > 0) {
      idx_dict[p_pos:p_end] = get_idx_dict(get_int_sub_array(idx, pos, p));
    } 
  } 
  
  return (idx_dict, dict_pos);
}


array[] int get_max(array[] int id, array[] int pos) {
  return get_max(id, pos, 0);
}

array[] int get_max_idx(array[] int id, array[] int pos) {
  return get_max(id, pos, 1);
}

array[] int get_max(array[] int id, array[] int pos, int of_idx) {
  int n = size(pos) - 1;
  array[n] int p_max = zeros_int_array(n);
  
  for (i in 1:n) {
    int n_i = get_pos_size(pos, i);
    
    if (n_i > 0) {
      array[n_i] int id_i = get_int_sub_array(id, pos, i);
      p_max[i] = max(id_i) + (of_idx ? 1 - min(id_i) : 0);
    }
  }
  
  return p_max;
}

array[] int create_pos(array[] int n_x) {
  int n = size(n_x);
  array[n + 1] int pos;
  pos[1] = 1;
  
  for (i in 1:n) {
    pos[i + 1] = pos[i] + n_x[i];  
  } 
  
  assert_equal(pos[n + 1] - 1, sum(n_x));
  
  return pos;
}

array[] int create_pos(array[] int n_x, array[] int sub_pos) {
  int n = size(sub_pos) - 1;
  array[n + 1] int pos;
  pos[1] = 1;
  
  for (i in 1:n) {
    pos[i + 1] = pos[i] + sum(get_int_sub_array(n_x, sub_pos, i));  
  } 
  
  assert_equal(pos[n + 1] - 1, sum(n_x));
  
  return pos;
}

array[] int create_pos(array[] int pos, int from, int to) {
  return create_pos(get_pos_size(pos)[from:to]);
}

tuple(int, int) get_pos(array[] int pos, int from, int to) {
  return (pos[from], pos[to + 1] - 1);
}

tuple(int, int) get_pos(array[] int pos, int n) {
  return get_pos(pos, n, n);
}

int get_pos_size(array[] int pos, int i) {
  return pos[i + 1] - pos[i];
}

array[] int get_pos_size(array[] int pos) {
  int n = size(pos) - 1;
  array[n] int sizes;
  
  for (i in 1:n) {
    sizes[i] = get_pos_size(pos, i);
  }
  
  return sizes;
} 

array[] int get_int_sub_array(array[] int full, array[] int pos, int n) {
  int start, end;
  (start, end) = get_pos(pos, n);

  return full[start:end];
}

array[] int get_int_sub_array(array[] int full, array[] int pos, int from, int to) {
  int from_start, from_end, to_start, to_end;
  (from_start, from_end) = get_pos(pos, from);
  (to_start, to_end) = get_pos(pos, to);

  return full[from_start:to_end];
}

vector get_sub_vector(vector full, array[] int pos, int n) {
  int start, end;
  (start, end) = get_pos(pos, n);

  return full[start:end];
}

array[] int get_min_pos(array[] int x, array[] int pos) {
  int n = size(pos) - 1;
  array[n] int min_pos;
  
  for (i in 1:n) {
    min_pos[i] = get_min_pos(x, pos, i);
  }
  
  return min_pos;
}

array[] int get_max_pos(array[] int x, array[] int pos) {
  int n = size(pos) - 1;
  array[n] int max_pos;
  
  for (i in 1:n) {
    max_pos[i] = get_max_pos(x, pos, i);
  }
  
  return max_pos;
}

int get_min_pos(array[] int x, array[] int pos, int n) {
  return min(get_int_sub_array(x, pos, n));
}

int get_max_pos(array[] int x, array[] int pos, int n) {
  return max(get_int_sub_array(x, pos, n));
}

void assert_equal(int x, int y) {
  if (x != y) {
    fatal_error("Equality assertion failed.");
  }
}

void assert_greater_than_or_equal(int x, int y) {
  if (x > y) {
    fatal_error("Greater than or equal assertion failed.");
  }
}