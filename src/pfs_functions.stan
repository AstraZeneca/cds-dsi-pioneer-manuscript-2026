// Simple regression model for the influence of tumors on surival. 
vector linear_tumor_stimulus(real intercept, vector coef, matrix covar) {
  return intercept + covar * coef; 
} 

// Combine influence of all tumors on survival and calculate probability of survival using a cloglog link function. 
vector calculate_progress_linear_prob(vector log_lambda, real tumor_intercept, vector tumor_coef, matrix tumor_covar) {
    real total_time_invar_tumor_stim = sum(linear_tumor_stimulus(tumor_intercept, tumor_coef, tumor_covar));
    
    return inv_cloglog(log_lambda + total_time_invar_tumor_stim);

}
// Hazard function given a base hazard and time-invariant covariates.  
vector calculate_linear_hazard(vector log_lambda, real tumor_intercept, vector tumor_coef, matrix tumor_covar) {
    real total_time_invar_tumor_stim = sum(linear_tumor_stimulus(tumor_intercept, tumor_coef, tumor_covar));
    
    return exp(log_lambda + total_time_invar_tumor_stim);
}

// Given PFS and tumor measures data, determine interval and right censoring for each patient. 
tuple(array[] int, array[] int) identify_censoring(array[] int pfs, array[] int death_week, array[] int n_patient_tumors, array[] int n_measures, array[] int t_measure) { 
  int n_patients = size(n_patient_tumors);
  array[n_patients] int interval_censored = rep_array(0, n_patients);
  array[n_patients] int right_censored = rep_array(1, n_patients);
  
  int t_pos = 1;
  int tumor_pos = 1;
    
    for (i in 1:n_patients) {
      if (death_week[i] > 0) {
        interval_censored[i] = death_week[i] - pfs[i] - 1; 
        right_censored[i] = 0;
      }
      
      for (j in 1:n_patient_tumors[i]) {
        int t_end = t_pos + n_measures[tumor_pos] - 1; 
        
        // for (t_index in 1:n_measures[tumor_pos]) {
        for (t_index in t_pos:t_end) {
          // if (t_pos + t_index - 1 > size(t_measure)) {
          //   print("i = ", i, ", j = ", j, ", t_index = ", t_index);
          //   print("t_pos = ", t_pos);
          // }
          
          int curr_t = t_measure[t_index]; 
          
          if (death_week[i] == 0 && curr_t > pfs[i]) { 
            // Progression actually happened between pfs[i] and the next measured interval. I'm taking the min here to use closest following
            // t; some tumors might not be observed for all t, so I don't want to arbitrarily use the last one's next t.
            interval_censored[i] = interval_censored[i] > 0 ? min(curr_t - pfs[i] - 1, interval_censored[i]) : curr_t - pfs[i] - 1; 
            
            right_censored[i] *= 0; // Found an observation after pfs[i] for _any_ of the tumors (hence the multiplication)
            
            break;
          } 
        }
        
        t_pos = t_end + 1;
        tumor_pos += 1;
      }
    }
    
    return (interval_censored, right_censored);
}

// Random PFS generator given conditional progression probability and obseration intervals. 
tuple(int, int, int, int) pfs_rng(vector prob, array[] int t) {
  int n_prob = rows(prob);
  int n_t = size(t);
  int actual_pfs = 0;
  int observed_pfs = 0;
  int right_censored = 0;
  int interval_censored = 0;
  int pfs_measure_index = 0;
  
  if (n_prob < t[n_t]) {
    reject("Insufficient number of probabilities provided.");
  }
 
  // Iterate over intervals until disease progression occurs 
  while (actual_pfs < n_prob && !bernoulli_rng(prob[actual_pfs + 1])) {
    actual_pfs += 1;
  }
 
  // Find the index in the measurement t array that corresponds to the interval observed after last progression free interval.
  while (pfs_measure_index < n_t && t[pfs_measure_index + 1] <= actual_pfs) {
    pfs_measure_index += 1;
  }
  
  right_censored = pfs_measure_index >= n_t; // We passed beyond the measurement t array so we must be right censored.
  observed_pfs = pfs_measure_index > 0 ? t[min(pfs_measure_index, n_t)] : 0; // Figure out which of the observed intervals would be the observed PFS
  interval_censored = !right_censored ? t[pfs_measure_index + 1] - observed_pfs - 1 : 0; // Figure how many intervals forward could be the true PFS 
  
  return (interval_censored, right_censored, observed_pfs, actual_pfs);
}

// Survival aggregated over all patients, S(t) = Pr[T > t], t \in {0,..., N} 
vector estimate_kaplan_meier(array[] int pfs, array[] int right_censored, int max_t) {
  vector[max_t + 1] s = rep_vector(1.0, max_t + 1);
  int n_pfs = size(pfs); // How many patients
  array[n_pfs] int sorted_pfs_idx = sort_indices_asc(pfs);
  int pfs_pos = 1;
  int n = n_pfs; // How many patients still haven't seen disease progression. 
  
  // For each time interval in 0..max_t see how many patiented exited and calculate proportion surviving.
  
  for (t in 0:max_t) {
    int ex = 0; // How many saw disease progression (exited) in current interval.
    real prev_s = t > 0 ? s[t] : 1.0;
    
    while ((n > 0) && (pfs_pos <= n_pfs) && (right_censored[sorted_pfs_idx[pfs_pos]] || (pfs[sorted_pfs_idx[pfs_pos]] <= t))) {
      ex += !right_censored[sorted_pfs_idx[pfs_pos]];
      pfs_pos += 1;
    }
   
    s[t + 1] = n > 0 ? prev_s * (n - ex) / n : 0.0;
    n -= ex;
  }
  
  return s; 
}  

// Calculate marginal probability of disease progression at every time interval, given conditional probabilities.
tuple(vector, vector) calculate_marginal_dp_prob(vector cond_pf_prob, int max_all_t) {
  vector[max_all_t] log_cond_pf_prob = log(cond_pf_prob[:max_all_t]);
  vector[max_all_t] log_1m_cond_pf_prob = log(1 - cond_pf_prob[:max_all_t]);
  vector[max_all_t] marginal_dp_prob;
  vector[max_all_t] dp_cdf;

  for (m in 1:max_all_t) {
    if (m > 1) {
      // marginal_dp_prob[m] = (1 - cond_pf_prob[m]) * prod(cond_pf_prob[1:(m - 1)]);
      marginal_dp_prob[m] = log_1m_cond_pf_prob[m] + sum(log_cond_pf_prob[1:(m - 1)]);
      dp_cdf[m] = exp(marginal_dp_prob[m]) + dp_cdf[m - 1]; 
    } else {
      marginal_dp_prob[m] = log_1m_cond_pf_prob[m];
      dp_cdf[m] = exp(marginal_dp_prob[m]);
    }
  }
  
  return(exp(marginal_dp_prob), fmax(0, 1 - dp_cdf));
}  

// Create (n_patients * n_tumors) x 2 matrix of each tumor's covariates from t = 1, 2. 
matrix prepare_early_tumors_design_matrix(vector tumor_size, array[] int n_patient_tumors, array[] int n_measures, array[] int n_screening_t) {
  matrix[sum(n_patient_tumors), 2] tumor_covar;
  int n_patients = size(n_patient_tumors);
  int tumor_pos = 1;
  int tumor_size_pos = 1;
  int covar_pos = 1;
  
  for (i in 1:n_patients) {
    int n_current_tumors = n_patient_tumors[i];
    
    for (j in 1:n_current_tumors) {
      int n_current_measures = n_measures[tumor_pos];
      int tumor_size_end = tumor_size_pos + n_current_measures - 1; 
      
      if (n_current_measures < 2) {
        reject("Two measures minimum needed per tumor.");
      }
      
      tumor_covar[covar_pos, ] = tumor_size[(tumor_size_pos + n_screening_t[tumor_pos] - 1):(tumor_size_pos + n_screening_t[tumor_pos] - 1 + 1)]';
      
      covar_pos += 1;
      tumor_size_pos = tumor_size_end + 1;
      tumor_pos += 1;
    }
  }
  
  return tumor_covar;
}