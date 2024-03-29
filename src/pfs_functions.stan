// Simple regression model for the influence of tumors on surival. 
vector linear_tumor_stimulus(vector intercept, matrix coef, matrix covar) {
  return intercept + rows_dot_product(covar, coef); 
} 

// Combine influence of all tumors on survival and calculate probability of survival using a cloglog link function. 
matrix calculate_progress_linear_prob(array[] int n_patient_tumors, vector log_lambda, vector tumor_intercept, matrix tumor_coef, matrix tumor_covar) {
    vector[rows(tumor_intercept)] total_time_invar_tumor_stim = 
      sum(linear_tumor_stimulus(tumor_intercept, tumor_coef, tumor_covar)) - tumor_intercept * mean(n_patient_tumors);
    
    return inv_cloglog(rep_matrix(log_lambda, rows(tumor_intercept)) + rep_matrix(total_time_invar_tumor_stim', rows(log_lambda)));
}

// Hazard function given a base hazard and time-invariant covariates.  
matrix calculate_linear_hazard(array[] int n_patient_tumors, vector log_lambda, vector tumor_intercept, matrix tumor_coef, matrix tumor_covar) {
    vector[rows(tumor_intercept)] total_time_invar_tumor_stim = 
      sum(linear_tumor_stimulus(tumor_intercept, tumor_coef, tumor_covar)) - tumor_intercept * mean(n_patient_tumors);
    
    return exp(rep_matrix(log_lambda, rows(tumor_intercept)) + rep_matrix(total_time_invar_tumor_stim', rows(log_lambda)));
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
        
        for (t_index in t_pos:t_end) {
          int curr_t = t_measure[t_index]; 
          
          if (curr_t > pfs[i]) { 
            // Progression actually happened between pfs[i] and the next measured interval. I'm taking the min here to use closest following
            // t; some tumors might not be observed for all t, so I don't want to arbitrarily use the last one's next t.
            interval_censored[i] = interval_censored[i] > 0 ? min(curr_t - pfs[i] - 1, interval_censored[i]) : curr_t - pfs[i] - 1; 
            
            right_censored[i] = 0; // Found an observation after pfs[i] for _any_ of the tumors
            
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
tuple(vector, array[] int, array[] int, array[] int) estimate_kaplan_meier(array[] int pfs, array[] int right_censored, int max_t) {
  int n_pfs = size(pfs); // How many patients
  array[n_pfs] int sorted_pfs_idx = sort_indices_asc(pfs);
  int pfs_pos = 1;
  int n = n_pfs; // How many patients still haven't seen disease progression. 
  
  vector[max_t + 1] s = rep_vector(1.0, max_t + 1);
  array[max_t + 1] int at_risk = rep_array(n, max_t + 1);
  array[max_t + 1] int n_right_censored = rep_array(0, max_t + 1);
  array[max_t + 1] int n_exited; // = rep_array(0, max_t + 1);
  
  // For each time interval in 0..max_t see how many patiented exited and calculate proportion surviving.
  
  for (t in 0:max_t) {
    // int ex = 0; // How many saw disease progression (exited) in current interval.
    n_exited[t + 1] = 0; // How many saw disease progression (exited) in current interval.
    real prev_s = t > 0 ? s[t] : 1.0;
    
    while ((n > 0) && (pfs_pos <= n_pfs) && (right_censored[sorted_pfs_idx[pfs_pos]] || (pfs[sorted_pfs_idx[pfs_pos]] <= t))) {
      n_exited[t + 1] += !right_censored[sorted_pfs_idx[pfs_pos]];
      n_right_censored[t + 1] += right_censored[sorted_pfs_idx[pfs_pos]];
     
      pfs_pos += 1;
    }
  
    s[t + 1] = n > 0 ? prev_s * (n - n_exited[t + 1]) / n : prev_s;
    at_risk[t + 1] = n; 
    n -= n_exited[t + 1] + n_right_censored[t + 1];
  }
  
  return (s, at_risk, n_right_censored, n_exited); 
}  

// Calculate marginal probability of disease progression at every time interval, given conditional probabilities.
tuple(vector, vector) calculate_marginal_dp_prob(vector cond_pf_prob, int max_all_t) {
  vector[max_all_t] log_cond_pf_prob = log(cond_pf_prob[:max_all_t]);
  vector[max_all_t] log_1m_cond_pf_prob = log(1 - cond_pf_prob[:max_all_t]);
  vector[max_all_t] marginal_dp_log_prob;
  vector[max_all_t] dp_cdf;

  for (m in 1:max_all_t) {
    if (m > 1) {
      // marginal_dp_prob[m] = (1 - cond_pf_prob[m]) * prod(cond_pf_prob[1:(m - 1)]);
      marginal_dp_log_prob[m] = log_1m_cond_pf_prob[m] + sum(log_cond_pf_prob[1:(m - 1)]);
      dp_cdf[m] = exp(marginal_dp_log_prob[m]) + dp_cdf[m - 1]; 
    } else {
      marginal_dp_log_prob[m] = log_1m_cond_pf_prob[m];
      dp_cdf[m] = exp(marginal_dp_log_prob[m]);
    }
  }
  
  return(exp(marginal_dp_log_prob), fmax(0, 1 - dp_cdf));
}  

// Create (n_patients * n_tumors) x max_tumors matrix of each tumor's covariates from t = 1, 2, .... 
matrix prepare_early_tumors_design_matrix(vector tumor_size, array[] int n_patient_tumors, array[] int n_measures, array[] int n_screening_t, int max_measures) {
  matrix[sum(n_patient_tumors), max_measures] tumor_covar = rep_matrix(0, sum(n_patient_tumors), max_measures);
  int n_patients = size(n_patient_tumors);
  int tumor_pos = 1;
  int tumor_size_pos = 1;
  int covar_pos = 1;
  
  for (i in 1:n_patients) {
    int n_current_tumors = n_patient_tumors[i];
    
    for (j in 1:n_current_tumors) {
      int n_current_measures = n_measures[tumor_pos];
      int tumor_size_end = tumor_size_pos + n_current_measures - 1; 
      
      int measures_found = min(n_current_measures, max_measures);
      
      tumor_covar[covar_pos, :measures_found] = 
        tumor_size[(tumor_size_pos + n_screening_t[tumor_pos] - 1):(tumor_size_pos + n_screening_t[tumor_pos] - 1 + measures_found - 1)]';
      
      covar_pos += 1;
      tumor_size_pos = tumor_size_end + 1;
      tumor_pos += 1;
    }
  }
  
  return tumor_covar;
}


vector calc_pch_loglik(array[] int pfs, array[] int right_uncensored, array[] int interval_censored, int ignore_interval_censoring, vector disease_progress_prob, int max_all_t) {
  int n_patients = size(pfs);
  vector[n_patients] lp = rep_vector(0, n_patients);
  
  int pfs_interval_pos = 1;
    
  for (i in 1:n_patients) {
    int observed_pfs_interval_end = pfs_interval_pos + pfs[i] - 1; 
    
    // These are the time intervals we are sure that the patient has progression free 
    lp[i] += bernoulli_lpmf(0 | disease_progress_prob[pfs_interval_pos:observed_pfs_interval_end]);
    
    int pfs_interval_end = observed_pfs_interval_end + right_uncensored[i] + interval_censored[i]; 
    int curr_interval_censored = ignore_interval_censoring ? 0 : interval_censored[i];
    vector[curr_interval_censored + right_uncensored[i]] interval_lp = rep_vector(0, curr_interval_censored + right_uncensored[i]);
  
    for (t in 1:(curr_interval_censored + right_uncensored[i])) {
      if (t > 1) { // We need to add more possible intervals that the patient remained progression free.
        interval_lp[t] = bernoulli_lpmf(0 | disease_progress_prob[(observed_pfs_interval_end + 1):(observed_pfs_interval_end + t - 1)]);
      }
      
      // If not right censored add pdf of disease progression. 
      if (right_uncensored[i]) {
        interval_lp[t] += bernoulli_lpmf(1 | disease_progress_prob[observed_pfs_interval_end + t]);
      }
    }
   
    if (curr_interval_censored > 0) {
      // There are more than one candidate true PFS: sum of the probabilities and then log.
      lp[i] += log_sum_exp(interval_lp); 
    } else if (right_uncensored[i]) {
      lp[i] += interval_lp[1]; // PFS not observed because of right censoring.
    }
    
    // If generating PFS, jump ahead to the beginning of the next patient's probs.
    pfs_interval_pos = (max_all_t > 0 ? pfs_interval_pos + max_all_t - 1 : pfs_interval_end) + 1;
  }
  
  return lp;
}
  
real pch_lpmf(array[] int y, array[] int right_uncensored, array[] int interval_censored, int ignore_interval_censoring, vector disease_progress_prob, int max_all_t) {
  return sum(calc_pch_loglik(y, right_uncensored, interval_censored, ignore_interval_censoring, disease_progress_prob, max_all_t));
} 