/**
 * Given PFS and tumor measures data, determine interval and right censoring for each patient. 
 * 
 * @param pfs Patient-level array of the number of weeks survived without disease progression.
 * @param death_week Patient-level array of what week death was observed.
 * @param n_patient_tumors Array with the number of tumors per patient.
 * @param n_measures The number of assessments per tumor.
 * @param t_measure The week each assessment was done.
 * @return Per patient, (Number of interval censoring intervals, indicator of right censoring).
 */
tuple(array[] int, array[] int) identify_censoring(
  array[] int pfs, array[] int death_week, array[] int n_patient_tumors, array[] int n_measures, array[] int t_measure) 
{ 
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

/** Survival aggregated over all patients, S(t) = Pr[T > t], t \in {0,..., N} 
 * For each time interval in 0..max_t see how many patiented exited and calculate proportion surviving.
 *
 * @param pfs The last observed week that was progression-free
 * @param right_censored Right censoring per patient
 * @param max_t The last interval to report Kaplan-Meier results
 * @return (Proportion surviving, Number at risk, Number right censored, Number for whom disease progressed) for each week
 */
tuple(vector, array[] int, array[] int, array[] int) estimate_kaplan_meier(array[] int pfs, array[] int right_censored, int max_t) {
  int n_pfs = size(pfs); // How many patients
  array[n_pfs] int sorted_pfs_idx = sort_indices_asc(pfs);
  int pfs_pos = 1;
  int n = n_pfs; // How many patients still haven't seen disease progression. 
  
  vector[max_t + 1] s = rep_vector(1.0, max_t + 1);
  array[max_t + 1] int at_risk = rep_array(n, max_t + 1);
  array[max_t + 1] int n_right_censored = rep_array(0, max_t + 1);
  array[max_t + 1] int n_exited = rep_array(0, max_t + 1); 
  
  for (t in 1:(max_t + 1)) {
    real prev_s = t > 1 ? s[t - 1] : 1.0;
    
    while ((n > 0) && (pfs_pos <= n_pfs) && (pfs[sorted_pfs_idx[pfs_pos]] <= t)) {
      if (t <= max_t) {
        // Remember that we define "pfs" as the last interval survived not the interval of exit.
        n_exited[t + 1] += !right_censored[sorted_pfs_idx[pfs_pos]]; 
      }
      
      n_right_censored[t] += right_censored[sorted_pfs_idx[pfs_pos]];
     
      pfs_pos += 1; 
      
    }
  
    s[t] = n > 0 ? prev_s * (n - n_exited[t]) / n : prev_s;
    at_risk[t] = n; 
    n -= n_exited[t] + n_right_censored[t];
  }
  
  return (s, at_risk, n_right_censored, n_exited); 
}  

/** Convert a ragged vector tumor size measures to a matrix aligned by measurement week
 *
 * @param tumor_size Tumor sizes
 * @param n_patient_tumors Array with the number of tumors per patient.
 * @param n_measures The number of assessments per tumor.
 * @param t_measure The week each assessment was done.
 * @param n_screening_t Number of observed pre-screening assessments per tumor
 * @param max_measures Number of tumor sizes to include in output, irrespective of when they are actually observed.
 * @return (Aligned matrix, Actual t used for the columns)
 */
tuple(matrix, array[,] int) create_aligned_matrix(vector tumor_size, array[] int n_measures, array[] int t_measure, array[] int n_screening_t, int max_measures) {
  int n_tumors = size(n_measures);
  int n_total_measures = sum(n_measures);
  array[n_total_measures] int t_sort_ind = sort_indices_asc(t_measure);
  array[n_total_measures] int within_tumor_ind;
  
  int last_within_ind = 0;
  int last_t = min(t_measure) - 1; 
  
  for(n in 1:n_total_measures) {
    if (t_measure[t_sort_ind[n]] != last_t) {
      last_within_ind += 1;
    }
    
    last_t = t_measure[t_sort_ind[n]];
    within_tumor_ind[t_sort_ind[n]] = last_within_ind;
  }

  int n_screening_col = max(n_screening_t);
  int n_aligned_col = max(max_measures, last_within_ind); 
  matrix[n_tumors, n_aligned_col] aligned_matrix = rep_matrix(0, n_tumors, n_aligned_col);
  array[n_tumors, n_aligned_col] int aligned_t = rep_array(0, n_tumors, n_aligned_col);
  
  int tumor_pos = 1;
  
  for (j in 1:n_tumors) {
    int tumor_end = tumor_pos + n_measures[j] - 1;
    
    aligned_matrix[j, within_tumor_ind[tumor_pos:tumor_end]] = tumor_size[tumor_pos:tumor_end]'; 
    aligned_t[j, within_tumor_ind[tumor_pos:tumor_end]] = t_measure[tumor_pos:tumor_end]; 
    
    tumor_pos = tumor_end + 1;
  }
  
  return (aligned_matrix[, n_screening_col:(max_measures + n_screening_col - 1)], aligned_t[, n_screening_col:(max_measures + n_screening_col - 1)]);
}

/** Create (n_patients * n_tumors) x max_measures matrix of each tumor's covariates from assessment **order** (not t) = 1, 2, ...., max_measures 
 * What this function actually returns is the last pre-screening measure and the (max_measures - 1) succeeding measures. The second output in the tuple
 * provides the intervals (weeks) of these max_measures columns.
 * 
 * @param tumor_size Tumor sizes
 * @param n_patient_tumors Array with the number of tumors per patient.
 * @param n_measures The number of assessments per tumor.
 * @param t_measure The week each assessment was done.
 * @param n_screening_t Number of observed pre-screening assessments per tumor
 * @param max_measures Number of tumor sizes to include in output, irrespective of when they are actually observed.
 * @return (Design matrix, Actual t used for the columns)
 */
tuple(matrix, array[,] int) prepare_early_tumors_design_matrix(
  vector tumor_size, array[] int n_patient_tumors, array[] int n_measures, array[] int t_measure, array[] int n_screening_t, int max_measures
) {
  matrix[sum(n_patient_tumors), max_measures] tumor_covar = rep_matrix(0, sum(n_patient_tumors), max_measures);
  array[sum(n_patient_tumors), max_measures] int tumor_covar_t = rep_array(min(t_measure) - 1, sum(n_patient_tumors), max_measures);
  int n_patients = size(n_patient_tumors);
  int tumor_pos = 1;
  int tumor_size_pos = 1;
  int covar_pos = 1;
  
  for (i in 1:n_patients) {
    int n_current_tumors = n_patient_tumors[i];
    int covar_end = covar_pos + n_current_tumors - 1; 
    int tumor_end = tumor_pos + n_current_tumors - 1;
    int n_patient_measures = sum(n_measures[tumor_pos:tumor_end]);
    int tumor_size_end = tumor_size_pos + n_patient_measures - 1;
   
    (tumor_covar[covar_pos:covar_end], tumor_covar_t[covar_pos:covar_end]) = create_aligned_matrix(
      tumor_size[tumor_size_pos:tumor_size_end], n_measures[tumor_pos:tumor_end], t_measure[tumor_size_pos:tumor_size_end], n_screening_t, max_measures
    );
    
    covar_pos = covar_end + 1;
    tumor_pos = tumor_end + 1;
    tumor_size_pos = tumor_size_end + 1;
  }
  
  return (tumor_covar, tumor_covar_t);
}

tuple(matrix, matrix, vector, vector) prepare_early_tumor_sums_covar(
  vector tumor_size, 
  array[] int n_patient_tumors, array[] int n_measures, array[] int t_measure, array[] int n_screening_t, int max_measures
) {
  int n_tumors = sum(n_patient_tumors);
  int n_patients = size(n_patient_tumors);
  matrix[n_tumors, 2] tumor_covar; 
  array[n_tumors, 2] int tumor_covar_t; 
  matrix[n_patients, 2] tumor_sum_covar; 
  matrix[n_patients, 2] uncentered_tumor_sum_covar; 
  
  (tumor_covar, tumor_covar_t) = prepare_early_tumors_design_matrix(tumor_size, n_patient_tumors, n_measures, t_measure, n_screening_t, max_measures);
  
  {
    int tumor_pos = 1;
    
    for (i in 1:n_patients) {
      int tumor_end = tumor_pos + n_patient_tumors[i] - 1;
      
      tumor_sum_covar[i] = ones_row_vector(n_patient_tumors[i]) * tumor_covar[tumor_pos:tumor_end];
      
      tumor_pos = tumor_end + 1;
    }
  }
  
  vector[2] tumor_sum_covar_mean;
  vector[2] tumor_sum_covar_sd;
  
  for (c in 1:2) {
    uncentered_tumor_sum_covar[, c] = tumor_sum_covar[, c];
    
    (tumor_sum_covar_mean[c], tumor_sum_covar_sd[c], tumor_sum_covar[, c]) = standardize_tumor_sizes(tumor_sum_covar[, c]);
    
    uncentered_tumor_sum_covar[, c] /= tumor_sum_covar_sd[c];
  }
  
  return (tumor_sum_covar, uncentered_tumor_sum_covar, tumor_sum_covar_mean, tumor_sum_covar_sd);
}
