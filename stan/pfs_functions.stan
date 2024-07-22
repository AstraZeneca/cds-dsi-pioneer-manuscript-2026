/** Simple regression model for the influence of tumors on surival. 
 *
 * @param intercept Vector of tumor-level log hazard ratio model.
 * @param coef Matrix of tumor-level (rows) log hazard ratio model coefficients for the effect of tumor sizes.
 * @param covar Design matrix
 * @return Tumor-level log hazard ratios.
 */
vector linear_tumor_stimulus(vector intercept, matrix coef, matrix covar) {
  return intercept + rows_dot_product(covar, coef); 
} 

/** Combine influence of all tumors on survival and calculate probability of survival using a cloglog link function. 
 * 
 * @param log_lambda Log of baseline hazard.
 * @param itumor_ntercept log hazard ratio model intercept.
 * @param tumor_coef Log hazard ratio model coefficients for the effect of tumor sizes.
 * @param tumor_covar Tumor size covariates. 
 * @return <number of intervals> vector of probabilities of disease progress. 
 */
vector calculate_progress_linear_prob(vector log_lambda, real tumor_intercept, row_vector tumor_coef, row_vector tumor_covar) {
    real total_time_invar_tumor_stim = linear_tumor_stimulus([ tumor_intercept ]', [ tumor_coef ], [ tumor_covar ])[1];
    
    return inv_cloglog(log_lambda + rep_vector(total_time_invar_tumor_stim, rows(log_lambda))); 
}

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

/** Random PFS generator given conditional progression probability and obseration intervals. 
 *
 * @param prob Vector of conditional probability of disease progression.
 * @param t Array of weeks when tumor assessments were actually conducted.
 * @return (is interval censored, is right censored, observed PFS, actual PFS)
 */
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

/** Survival aggregated over all patients, S(t) = Pr[T > t], t \in {0,..., N} 
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
  array[max_t + 1] int n_exited; // = rep_array(0, max_t + 1);
  
  // For each time interval in 0..max_t see how many patiented exited and calculate proportion surviving.
  
  for (t in 0:max_t) {
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

/** Calculate marginal probability of disease progression at every time interval, given conditional probabilities.
 *
 * @param cond_pf_prob Conditional probability of disease progression
 * @param max_all_t The last week observed in all the data
 * @return (Marginal probabilities, Disease progression empirical CDF) for each week
 */
tuple(vector, vector) calculate_marginal_dp_prob(vector cond_pf_prob, int max_all_t) {
  vector[max_all_t] log_cond_pf_prob = log(cond_pf_prob[:max_all_t]);
  vector[max_all_t] log_1m_cond_pf_prob = log(1 - cond_pf_prob[:max_all_t]);
  vector[max_all_t] marginal_dp_log_prob;
  vector[max_all_t] dp_cdf;

  for (m in 1:max_all_t) {
    if (m > 1) {
      marginal_dp_log_prob[m] = log_1m_cond_pf_prob[m] + sum(log_cond_pf_prob[1:(m - 1)]);
      dp_cdf[m] = exp(marginal_dp_log_prob[m]) + dp_cdf[m - 1]; 
    } else {
      marginal_dp_log_prob[m] = log_1m_cond_pf_prob[m];
      dp_cdf[m] = exp(marginal_dp_log_prob[m]);
    }
  }
  
  return(exp(marginal_dp_log_prob), fmax(0, 1 - dp_cdf));
}  

/** Calculate the piecewise-constant proportional hazard log-likelihood. This returns the patients vector of log-likelihoods as opposed to the following
 * pch_lpmf() function. 
 *
 * @param pfs Observed number of weeks without disease progression.
 * @param right_censored Is right censored?
 * @param interval_censored Number of weeks over which we have interval censoring.
 * @param ignore_interval_censoring Treat `pfs` as the actual PFS.
 * @param disease_progress_prob Conditional probability of disease progress at every interval.
 * @param max_all_t The latest week assessment is done in all the data.
 * @param patient_2nd_t The week in which the first post-treatment assessment was done.
 * @return Vector of patient-level log-likelihood.
 */
vector calc_pch_loglik(
  array[] int pfs, 
  array[] int right_uncensored, array[] int interval_censored, int ignore_interval_censoring, 
  vector disease_progress_prob, int max_all_t, array[] int patient_2nd_t 
) 
{
  int n_patients = size(pfs);
  vector[n_patients] lp = rep_vector(0, n_patients);
  
  int pfs_interval_pos = 1;
    
  for (i in 1:n_patients) {
    int observed_pfs_interval_end = pfs_interval_pos + pfs[i] - 1; 
    
    // Ignoring intervals that were guaranteed for the patient to have survived because of the inclusion criteria in this meta-analysis (not the the original trials).
    pfs_interval_pos += patient_2nd_t[i] - 1;
    
    int unobs_pfs_interval_pos = max(observed_pfs_interval_end + 1, pfs_interval_pos);
   
    if (pfs_interval_pos <= observed_pfs_interval_end) { // By incrementing by patient_2nd_t we can end up outside the observed range.
      // These are the time intervals we are sure that the patient was progression free 
      lp[i] += bernoulli_lpmf(0 | disease_progress_prob[pfs_interval_pos:observed_pfs_interval_end]);
    }
    
    int pfs_interval_end = observed_pfs_interval_end + right_uncensored[i] + interval_censored[i]; 
    int curr_interval_censored = ignore_interval_censoring ? 0 : interval_censored[i];
    // vector[curr_interval_censored + right_uncensored[i]] interval_lp = rep_vector(0, curr_interval_censored + right_uncensored[i]);
    int unobs_size = pfs_interval_end - unobs_pfs_interval_pos + 1;
    vector[unobs_size] interval_lp = rep_vector(0, unobs_size);
 
    // The point of this loop is marginalize over all the potential intervals of progression, due to interval censoring. 
    // for (t in 1:(curr_interval_censored + right_uncensored[i])) {
    for (t in 1:unobs_size) {
      if (t > 1) { // We need to add more possible intervals that the patient remained progression free.
        // interval_lp[t] = bernoulli_lpmf(0 | disease_progress_prob[(observed_pfs_interval_end + 1):(observed_pfs_interval_end + t - 1)]);
        interval_lp[t] = bernoulli_lpmf(0 | disease_progress_prob[unobs_pfs_interval_pos:(unobs_pfs_interval_pos - 2 + t)]);
      }
      
      // If not right censored add pdf of disease progression. 
      if (right_uncensored[i]) {
        interval_lp[t] += bernoulli_lpmf(1 | disease_progress_prob[unobs_pfs_interval_pos + t - 1]);
      }
    }
   
    if (curr_interval_censored > 0) {
      // There are more than one candidate true PFS: sum of the probabilities and then log.
      lp[i] += log_sum_exp(interval_lp); // - log(curr_interval_censored + 1); 
    } else if (right_uncensored[i]) {
      lp[i] += interval_lp[1]; // PFS not observed because of right censoring.
    }
    
    // If generating PFS, jump ahead to the beginning of the next patient's probs.
    pfs_interval_pos = (max_all_t > 0 ? pfs_interval_pos + max_all_t - (patient_2nd_t[i] - 1) - 1 : pfs_interval_end) + 1;
  }
  
  return lp;
}
/** Calculate the piecewise-constant proportional hazard log-likelihood. This returns the patients vector of log-likelihoods as opposed to the following
 * pch_lpmf() function. 
 *
 * @param pfs Observed number of weeks without disease progression.
 * @param right_censored Is right censored?
 * @param interval_censored Number of weeks over which we have interval censoring.
 * @param ignore_interval_censoring Treat `pfs` as the actual PFS.
 * @param log_cond_prob_progress Log conditional probability of disease progress at every interval.
 * @param max_all_t The latest week assessment is done in all the data.
 * @param patient_2nd_t The week in which the first post-treatment assessment was done.
 * @return Vector of patient-level log-likelihood.
 */
vector calc_pch_loglik2(
  array[] int pfs, 
  array[] int right_uncensored, array[] int interval_censored, int ignore_interval_censoring, 
  vector log_cond_prob_progress, int max_all_t, array[] int patient_2nd_t 
) 
{
  int n_patients = size(pfs);
  vector[n_patients] lp = rep_vector(0, n_patients);
  
  int pfs_interval_pos = 1;
    
  for (i in 1:n_patients) {
    int observed_pfs_interval_end = pfs_interval_pos + pfs[i] - 1; 
    
    // Ignoring intervals that were guaranteed for the patient to have survived because of the inclusion criteria in this meta-analysis (not the the original trials).
    pfs_interval_pos += patient_2nd_t[i] - 1;
    
    int unobs_pfs_interval_pos = max(observed_pfs_interval_end + 1, pfs_interval_pos);
   
    if (pfs_interval_pos <= observed_pfs_interval_end) { // By incrementing by patient_2nd_t we can end up outside the observed range.
      // These are the time intervals we are sure that the patient was progression free 
      // lp[i] += bernoulli_lpmf(0 | disease_progress_prob[pfs_interval_pos:observed_pfs_interval_end]);
      lp[i] += sum(log_cond_prob_progress[pfs_interval_pos:observed_pfs_interval_end]);
    }
    
    int pfs_interval_end = observed_pfs_interval_end + right_uncensored[i] + interval_censored[i]; 
    int curr_interval_censored = ignore_interval_censoring ? 0 : interval_censored[i];
    int unobs_size = pfs_interval_end - unobs_pfs_interval_pos + 1;
    vector[unobs_size] interval_lp = rep_vector(0, unobs_size);
 
    // The point of this loop is marginalize over all the potential intervals of progression, due to interval censoring. 
    for (t in 1:unobs_size) {
      if (t > 1) { // We need to add more possible intervals that the patient remained progression free.
        // interval_lp[t] = bernoulli_lpmf(0 | disease_progress_prob[unobs_pfs_interval_pos:(unobs_pfs_interval_pos - 2 + t)]);
        interval_lp[t] = sum(log_cond_prob_progress[unobs_pfs_interval_pos:(unobs_pfs_interval_pos - 2 + t)]); 
      }
      
      // If not right censored add pdf of disease progression. 
      if (right_uncensored[i]) {
        // interval_lp[t] += bernoulli_lpmf(1 | disease_progress_prob[unobs_pfs_interval_pos + t - 1]);
        interval_lp[t] += log1m_exp(log_cond_prob_progress[unobs_pfs_interval_pos + t - 1]); 
      }
    }
   
    if (curr_interval_censored > 0) {
      // There are more than one candidate true PFS: sum of the probabilities and then log.
      lp[i] += log_sum_exp(interval_lp); // - log(curr_interval_censored + 1); 
    } else if (right_uncensored[i]) {
      lp[i] += interval_lp[1]; // PFS not observed because of right censoring.
    }
    
    // If generating PFS, jump ahead to the beginning of the next patient's probs.
    pfs_interval_pos = (max_all_t > 0 ? pfs_interval_pos + max_all_t - (patient_2nd_t[i] - 1) - 1 : pfs_interval_end) + 1;
  }
  
  return lp;
}

vector calc_comp_risk_pch_loglik(
  array[] int last_unclass_week,
  array[] int event_cause,
  array[] int right_censored,
  matrix log_cond_prob_surv,
  int max_confresp_week
) 
{
  int n_patients = size(last_unclass_week);
  int n_causes = cols(log_cond_prob_surv);
  vector[n_patients] lp = rep_vector(0, n_patients);
  
  int interval_pos = 1;
    
  for (i in 1:n_patients) {
    int interval_end = interval_pos + last_unclass_week[i] - right_censored[i]; 
    
    if (right_censored[i]) {
      if (last_unclass_week[i] > 0) {
        for (k in 1:n_causes) {
          lp[i] += sum(log_cond_prob_surv[interval_pos:interval_end, k]); 
        }
      }
    } else {
      lp[i] = sum(log_cond_prob_surv[interval_pos:(interval_end - 1), event_cause[i]]) + log1m_exp(log_cond_prob_surv[interval_end, event_cause[i]]);
    }
    
    interval_pos += max_confresp_week; 
  }
  
  return lp;
}
 
/** This is used to provide and easy to use Stan distribution. It just sums the log-probs. 
 *
 * @param y Observed number of weeks without disease progression.
 * @param right_censored Is right censored?
 * @param interval_censored Number of weeks over which we have interval censoring.
 * @param ignore_interval_censoring Treat `pfs` as the actual PFS.
 * @param disease_progress_prob Conditional probability of disease progress at every interval.
 * @param max_all_t The latest week assessment is done in all the data.
 * @param patient_2nd_t The week in which the first post-treatment assessment was done.
 * @return Vector of patient-level log-likelihood.
 */
real pch_lpmf(
  array[] int y, 
  array[] int right_uncensored, array[] int interval_censored, int ignore_interval_censoring, vector disease_progress_prob, int max_all_t, array[] int patient_2nd_t
) {
  return sum(calc_pch_loglik(y, right_uncensored, interval_censored, ignore_interval_censoring, disease_progress_prob, max_all_t, patient_2nd_t));
}

/** This is used to provide and easy to use Stan distribution. It just sums the log-probs. 
 *
 * @param y Observed number of weeks without disease progression.
 * @param right_censored Is right censored?
 * @param interval_censored Number of weeks over which we have interval censoring.
 * @param ignore_interval_censoring Treat `pfs` as the actual PFS.
 * @param log_cond_prob_progress Log conditional probability of disease progress at every interval.
 * @param max_all_t The latest week assessment is done in all the data.
 * @param patient_2nd_t The week in which the first post-treatment assessment was done.
 * @return Vector of patient-level log-likelihood.
 */
real pch2_lpmf(
  array[] int y, 
  array[] int right_uncensored, array[] int interval_censored, int ignore_interval_censoring, vector log_cond_prob_progress, int max_all_t, array[] int patient_2nd_t
) {
  return sum(calc_pch_loglik2(y, right_uncensored, interval_censored, ignore_interval_censoring, log_cond_prob_progress, max_all_t, patient_2nd_t));
}

real comp_risk_pch_lpmf(
  array[] int last_unclass_week,
  array[] int event_cause,
  array[] int right_censored,
  matrix log_cond_prob_surv,
  int max_confresp_week
) {
  return sum(calc_comp_risk_pch_loglik(last_unclass_week, event_cause, right_censored, log_cond_prob_surv, max_confresp_week));
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

int calc_n_tumor_covar_col(int tumor_hazard_type) {
  int n_covar_col = 0; // number of columns in covariates design matrix, excluding the intercept. 
  
  if (tumor_hazard_type == 1) {
    n_covar_col = 2;
  } else if (tumor_hazard_type == 2) {
    n_covar_col = 3;
  } else if (tumor_hazard_type == 3) {
    n_covar_col = 1;
  } else if (tumor_hazard_type == 5) {
    n_covar_col = 5;
  }
  
  return n_covar_col;
}

tuple(matrix, matrix, vector, vector) prepare_early_tumor_sums_covar(
  vector tumor_size, 
  int n_tumors, array[] int n_patient_tumors, array[] int n_measures, array[] int t_measure, array[] int n_screening_t, int max_measures
) {
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
      
      tumor_sum_covar[i] = rep_row_vector(1, n_patient_tumors[i]) * tumor_covar[tumor_pos:tumor_end];
      
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

tuple(array[,] int, array[] int, matrix, matrix, vector, vector) prepare_early_tumors_covar(
  vector tumor_size, 
  int tumor_hazard_type, 
  int n_tumors, array[] int n_patient_tumors, array[] int n_measures, array[] int t_measure, array[] int n_screening_t, int max_measures
) {
  int n_patients = size(n_patient_tumors);
  int n_covar_col = calc_n_tumor_covar_col(tumor_hazard_type);
  
  array[n_tumors, max_measures] int tumor_covar_t; // ts (weeks) of the assessments used in the covar design matrix
  array[n_patients] int patient_max_2nd_tumor_t; // What t is the second assessment in the covar design matrix 
  matrix[n_tumors, n_covar_col] tumor_covar; // This is the design matrix with the covar in the first two (or _n_) assessments.
  matrix[n_tumors, n_covar_col] uncentered_tumor_covar; // Just scaled
  vector[n_covar_col] tumor_covar_mean;
  vector[n_covar_col] tumor_covar_sd = rep_vector(0, n_covar_col);
  
  { // Preparing all the tumor covariates data 
    tuple(matrix[n_tumors, 2], array[n_tumors, 2] int) prep_res = 
      prepare_early_tumors_design_matrix(tumor_size, n_patient_tumors, n_measures, t_measure, n_screening_t, max_measures);
      
    tumor_covar_t = prep_res.2;
    
    int tumor_pos = 1;
    
    for (i in 1:n_patients) {
      int tumor_end = tumor_pos + n_patient_tumors[i] - 1;
      
      patient_max_2nd_tumor_t[i] = max(tumor_covar_t[tumor_pos:tumor_end, 2]);
      
      tumor_pos = tumor_end + 1;
    }
    
    if (tumor_hazard_type > 0) {
      if (tumor_hazard_type == 1 || tumor_hazard_type == 2 || tumor_hazard_type == 5) {  
        tumor_covar[, 1:2] = prep_res.1;
        
        if (tumor_hazard_type == 2) {
          tumor_covar[, 3] = tumor_covar[, 2] .* tumor_covar[, 1];
        }
        
        if (tumor_hazard_type == 3) {
          tumor_covar[, 3] = tumor_covar[, 2] ./ tumor_covar[, 1];
        }
        
        if (tumor_hazard_type == 5) {
          tumor_covar[, 4] = tumor_covar[, 1]^2;
          tumor_covar[, 5] = tumor_covar[, 2]^2;
        }
        
        for (c in 1:n_covar_col) {
          uncentered_tumor_covar[, c] = tumor_covar[, c];
          
          tuple(real, real, vector[n_tumors]) standardize_results = standardize_tumor_sizes(tumor_covar[, c]);
          tumor_covar_mean[c] = standardize_results.1;
          tumor_covar_sd[c] = standardize_results.2;
          tumor_covar[, c] = standardize_results.3;
          
          uncentered_tumor_covar[, c] /= tumor_covar_sd[c];
        }
       
        if (tumor_hazard_type == 2 || tumor_hazard_type == 5) {
          tumor_covar[, 3] = rep_vector(0, n_tumors);
        }
      } else if (tumor_hazard_type == 3) {
        matrix[n_tumors, 2] covar = prep_res.1;
        
        tumor_covar[, 1] = (covar[, 2] - covar[, 1]) ./ covar[, 1];
        uncentered_tumor_covar[, 1] = tumor_covar[, 1];
      } 
    } 
  }
  
  return (tumor_covar_t, patient_max_2nd_tumor_t, tumor_covar, uncentered_tumor_covar, tumor_covar_mean, tumor_covar_sd);
}

// Calculate patient-interval conditional probability of disease progression.
tuple(vector, vector) calc_tumor_stim(
  matrix tumor_covar,
  int tumor_hazard_type, array[] int n_patient_tumors,
  array[] int patient_trial,
  array[] int tumor_location,
  real tumor_stim_pop_intercept, vector tumor_stim_trial_intercept, vector tumor_stim_location_intercept, 
  row_vector tumor_stim_pop_coef, matrix tumor_stim_trial_coef, matrix tumor_stim_location_coef
) {
  int n_patients = size(n_patient_tumors);
  int n_covar_col = cols(tumor_covar);
  
  vector[n_patients] total_time_invar_tumor_stim; // This is a sum of the contribution of all a patient's tumors to their hazard 
  vector[n_patients] total_time_invar_tumor_stim_no_intercept; // This used outside the model, so don't delete it.
  
  int tumor_pos = 1;

  for (i in 1:n_patients) {
    int tumor_end = tumor_pos + n_patient_tumors[i] - 1;
    array[n_patient_tumors[i]] int patient_tumor_locations = tumor_location[tumor_pos:tumor_end];
    
    vector[n_patient_tumors[i]] patient_stim_intercept = 
      tumor_stim_pop_intercept + tumor_stim_trial_intercept[patient_trial[i]] + tumor_stim_location_intercept[patient_tumor_locations];
      
    matrix[n_patient_tumors[i], n_covar_col] patient_stim_coef =
      rep_matrix(tumor_stim_pop_coef + tumor_stim_trial_coef[patient_trial[i]], n_patient_tumors[i]) + tumor_stim_location_coef[patient_tumor_locations];
     
    vector[n_patient_tumors[i]] tumor_stim = linear_tumor_stimulus(patient_stim_intercept, patient_stim_coef, tumor_covar[tumor_pos:tumor_end]);
    total_time_invar_tumor_stim[i] = tumor_hazard_type > 0 ? sum(tumor_stim) : 0;
    total_time_invar_tumor_stim_no_intercept[i] = 
      tumor_hazard_type > 0 ? sum(linear_tumor_stimulus(rep_vector(0, n_patient_tumors[i]), patient_stim_coef, tumor_covar[tumor_pos:tumor_end])) : 0;
   
    tumor_pos = tumor_end + 1;
  }

  return (total_time_invar_tumor_stim, total_time_invar_tumor_stim_no_intercept);    
}

// Calculate patient-interval conditional probability of disease progression.
tuple(vector, vector, vector) calc_disease_progress_pred_from_early_tumors(
  array[] int pfs, matrix tumor_covar,
  int tumor_hazard_type, array[] int n_patient_tumors,
  array[] int patient_trial,
  array[] int tumor_location,
  array[] int right_uncensored, array[] int interval_censored, 
  array[] vector log_trial_lambda,
  real tumor_stim_pop_intercept, vector tumor_stim_trial_intercept, vector tumor_stim_location_intercept, 
  row_vector tumor_stim_pop_coef, matrix tumor_stim_trial_coef, matrix tumor_stim_location_coef,
  int gen_pfs, int max_all_t
) {
  int n_patients = size(n_patient_tumors);
  int n_covar_col = cols(tumor_covar);
  int n_time_periods = gen_pfs ? n_patients * max_all_t : sum(pfs) + sum(right_uncensored) + sum(interval_censored);
  
  vector[n_patients] total_time_invar_tumor_stim; // This is a sum of the contribution of all a patient's tumors to their hazard 
  vector[n_patients] total_time_invar_tumor_stim_no_intercept; // This used outside the model, so don't delete it.
  
  (total_time_invar_tumor_stim, total_time_invar_tumor_stim_no_intercept) = calc_tumor_stim(
    tumor_covar, tumor_hazard_type, n_patient_tumors, patient_trial, tumor_location, 
    tumor_stim_pop_intercept, tumor_stim_trial_intercept, tumor_stim_location_intercept, tumor_stim_pop_coef, tumor_stim_trial_coef, tumor_stim_location_coef 
  );
  
  vector[n_time_periods] disease_progress_pred;
  
  int pfs_interval_pos = 1;

  for (i in 1:n_patients) {
    int n_intervals = gen_pfs ? max_all_t : pfs[i] + right_uncensored[i] + interval_censored[i];
    int pfs_interval_end = pfs_interval_pos + n_intervals - 1; 
   
    disease_progress_pred[pfs_interval_pos:pfs_interval_end] = log_trial_lambda[patient_trial[i], 1:n_intervals] + total_time_invar_tumor_stim[i];
    
    pfs_interval_pos = pfs_interval_end + 1;
  }

  return (total_time_invar_tumor_stim, total_time_invar_tumor_stim_no_intercept, disease_progress_pred);    
}
