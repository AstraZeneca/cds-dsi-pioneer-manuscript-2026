functions {
  #include "functions.stan"
 
  tuple(array[] int, array[] int, array[] int, array[] int, array[] int) cutoff_visits(
    int cutoff_calendar_day, array[] int patient_calendar_day, array[] int t_measure, array[] int t_day_measure, array[] int patient_tumor_measure_pos
  ) {
    int n_patients = size(patient_calendar_day);
    
    // 0 is the default sentinel value if last visit is negative 
    array[n_patients] int last_visit = rep_array(0, n_patients); 
    array[n_patients] int first_testing_visit_week = rep_array(0, n_patients);
    array[n_patients] int first_testing_visit_day = rep_array(0, n_patients);
    
    // array[n_patients] int last_testing_visit_calendar_week;
    array[n_patients] int last_testing_visit_calendar_day;
    
    for (i in 1:n_patients) {
      int t_measure_pos = patient_tumor_measure_pos[i]; 
      int t_measure_end = patient_tumor_measure_pos[i + 1] - 1; 
      int n_patient_measures = t_measure_end - t_measure_pos + 1;
      int patient_cutoff_study_day = calendar_date_to_study_date(patient_calendar_day[i], cutoff_calendar_day); 
      
      array[n_patient_measures] int patient_t_measure = t_measure[t_measure_pos:t_measure_end];
      array[n_patient_measures] int patient_t_day_measure = t_day_measure[t_measure_pos:t_measure_end];
      array[n_patient_measures] int patient_measure_t_sort_idx = sort_indices_asc(patient_t_day_measure);
      int t_idx = 0;
      
      while (t_idx < n_patient_measures && patient_t_day_measure[patient_measure_t_sort_idx[t_idx + 1]] <= patient_cutoff_study_day) {
        t_idx += 1;
      }
      
      if (t_idx > 0) {  
        last_visit[i] = t_measure[t_measure_pos:t_measure_end][patient_measure_t_sort_idx[t_idx]]; 
      }
      
      if (t_idx < n_patient_measures) {
        first_testing_visit_week[i] = patient_t_measure[patient_measure_t_sort_idx[t_idx + 1]];
        first_testing_visit_day[i] = patient_t_day_measure[patient_measure_t_sort_idx[t_idx + 1]];
      }
      
      // last_testing_visit_calendar_week[i] = patient_calendar_day[i] + patient_t_measure[patient_measure_t_sort_idx[n_patient_measures]];
      last_testing_visit_calendar_day[i] = patient_calendar_day[i] + patient_t_day_measure[patient_measure_t_sort_idx[n_patient_measures]] - 1;
    }
    
    return (last_visit, first_testing_visit_week, first_testing_visit_day, last_testing_visit_calendar_day);
  } 
 
  tuple(array[] int, array[] int, array[] int) cutoff_surv_data(array[] int last_visit, array[] int event_week, array[] int right_censored, array[] int interval_censored) {
    int n_patients = size(last_visit);
    
    array[n_patients] int new_event_week = event_week;
    array[n_patients] int new_right_censored = right_censored;
    array[n_patients] int new_interval_censored = interval_censored;
    
    for (i in 1:n_patients) {
      if (last_visit[i] < event_week[i] + interval_censored[i] + (1 - right_censored[i])) {
        new_right_censored[i] = 1;
        new_interval_censored[i] = 0;
        new_event_week[i] = last_visit[i]; 
      }
    }
    
    return(new_event_week, new_right_censored, new_interval_censored);
  } 

  /** Calculate the number of cuts that represent futures.
   */
  int calc_n_oos_log_lik(array[] int sorted_last_visit_day, int cutoff_calendar_day, int cutoff_calendar_day_increment) {
    int n_oos_log_lik = 0;
    int n_patients = size(sorted_last_visit_day); 
    int n_remaining_testing_patients = n_patients; 
    int curr_cutoff_calendar_day = cutoff_calendar_day;
    int curr_patient_idx = 1;
    
    while (curr_patient_idx <= n_patients) {
      while (curr_patient_idx <= n_patients && (sorted_last_visit_day[curr_patient_idx] <= curr_cutoff_calendar_day || sorted_last_visit_day[curr_patient_idx] <= 0)) {
        curr_patient_idx += 1;
      }
      
      if (curr_patient_idx <= n_patients) {
        n_oos_log_lik += 1;
      }
      
      curr_cutoff_calendar_day += cutoff_calendar_day_increment - 1;
    }
   
    return n_oos_log_lik; 
  }
 
  /** Get the indices within the array of sorted last visit that will be used in all the LFO cuts.
   *
   * Each such index will indicate the first patient (in the sorted array) to be in each cut. In each future cut, the patients are a subset of the previous cut's
   * patients.
   */
  array[] int get_oos_patients_idx(int n_oos_log_lik, array[] int sorted_last_visit_day, int cutoff_calendar_day, int cutoff_calendar_day_increment) {
    array[n_oos_log_lik] int patient_idx;
    int n_patients = size(sorted_last_visit_day); 
    int n_remaining_testing_patients = n_patients; 
    int curr_cutoff_calendar_day = cutoff_calendar_day;
    int curr_patient_idx = 1;
    int patient_idx_pos = 1;
    
    while (curr_patient_idx <= n_patients) {
      // while (last_visit_day_sort_idx[curr_patient_idx] <= curr_cutoff_calendar_day && curr_patient_idx <= n_patients && n_remaining_patients > 0) {
      while (curr_patient_idx <= n_patients && sorted_last_visit_day[curr_patient_idx] <= curr_cutoff_calendar_day) {
        curr_patient_idx += 1;
      }
      
      if (curr_patient_idx <= n_patients) {
        patient_idx[patient_idx_pos] = curr_patient_idx;
      }
      
      curr_cutoff_calendar_day += cutoff_calendar_day_increment - 1;
    }
    
    return patient_idx; 
  } 
}

data {
  #include "data.stan"
  
  int<lower = 1> cutoff_calendar_day;
  int<lower = 1> cutoff_calendar_day_increment;
}

transformed data {
  #include "transformed_data.stan"
  
  array[n_patients] int<lower = min(t_measure), upper = max(t_measure)> cutoff_last_visit_week;
  array[n_patients] int<lower = min(t_measure), upper = max(t_measure)> after_cutoff_first_visit_week;
  array[n_patients] int<lower = min(t_day_measure), upper = max(t_day_measure)> after_cutoff_first_visit_day;
  // array[n_patients] int<lower = min(t_measure), upper = max(t_measure)> after_cutoff_last_visit_week;
  array[n_patients] int<lower = 1> after_cutoff_last_visit_calendar_day;
  
  (cutoff_last_visit_week, after_cutoff_first_visit_week, after_cutoff_first_visit_day, after_cutoff_last_visit_calendar_day) = 
    cutoff_visits(cutoff_calendar_day, calendar_day, t_measure, t_day_measure, patient_tumor_measure_pos);
  
  int<lower = 0, upper = n_patients> n_training_patients = 0;
  
  for (i in 1:n_patients) {
    if (cutoff_last_visit_week[i] > 0) {
      n_training_patients += 1;
    }
  }
  
  print("n_training_patients = ", n_training_patients);
  
  array[n_training_patients] int<lower = 1, upper = n_patients> training_patients;
  
  {
    int training_pos = 1;
    
    for (i in 1:n_patients) {
      if (cutoff_last_visit_week[i] > 0) {
        training_patients[training_pos] = i;
        training_pos += 1;
      }
    }
  }
  
  array[n_training_patients] int<lower = 0> training_last_unclassified_response_week = last_unclassified_response_week[training_patients];
  array[n_training_patients] int<lower = 0, upper = 1> training_confirmed_response_censored = confirmed_response_censored[training_patients];
  array[n_training_patients] int<lower = 0> training_confirmed_response_interval_censored = confirmed_response_interval_censored[training_patients];
  
  (training_last_unclassified_response_week, training_confirmed_response_censored, training_confirmed_response_interval_censored) =
    cutoff_surv_data(
      cutoff_last_visit_week[training_patients], training_last_unclassified_response_week, training_confirmed_response_censored, training_confirmed_response_interval_censored
    ); 
    
  array[n_training_patients] int<lower = 0> training_pfs = pfs[training_patients]; // How many periods after baseline did patient survive.
  array[n_training_patients] int<lower = 0, upper = 1> training_right_censored = right_censored[training_patients];
  array[n_training_patients] int<lower = 0> training_interval_censored = interval_censored[training_patients];

  (training_pfs, training_right_censored, training_interval_censored) =
    cutoff_surv_data(cutoff_last_visit_week[training_patients], training_pfs, training_right_censored, training_interval_censored);
    
    
  array[n_patients] int<lower = 1, upper = n_patients> after_cutoff_last_visit_calendar_day_sort_idx = sort_indices_asc(after_cutoff_last_visit_calendar_day);
    
  int<lower = 0> n_oos_log_lik = 
    calc_n_oos_log_lik(after_cutoff_last_visit_calendar_day[after_cutoff_last_visit_calendar_day_sort_idx], cutoff_calendar_day, cutoff_calendar_day_increment);
  
  array[n_oos_log_lik] int<lower = 1, upper = n_patients> n_pfs_testing_patients = 
    get_oos_patients_idx(n_oos_log_lik, after_cutoff_last_visit_calendar_day[after_cutoff_last_visit_calendar_day_sort_idx], cutoff_calendar_day, cutoff_calendar_day_increment);
}

parameters {
  #include "parameters.stan"
}

transformed parameters {
  #include "transformed_parameters.stan"
  
  matrix[n_training_patients, n_causes] training_patient_response_lp; 
  
  training_patient_response_lp[, 1] = calc_pch_loglik(
    training_pfs, training_right_censored, training_interval_censored, pfs_ignore_interval_censoring, log_cond_prob_surv[1, training_patients]
  );
  training_patient_response_lp[, 2] = calc_pch_loglik(
    training_pfs, training_right_censored, training_interval_censored, pfs_ignore_interval_censoring, log_cond_prob_surv[2, training_patients]
  );
}

model {
  #include "priors.stan"
  
  // Likelihood
  
  if (fit_data) {
    // Confirmed response model
    
    profile("crcr loglik") {
      target += reduce_sum(
        partial_sum_crcr_lupmf, training_last_unclassified_response_week, crcr_grain_size,
        confirmed_response_cause[training_patients], 
        training_confirmed_response_censored, 
        crcr_ignore_interval_censoring ? zeros_int_array(n_training_patients) : training_confirmed_response_interval_censored, 
        log_crcr_cond_prob_surv
      );
    }
    
    // PFS model
    
    for (idx in 1:n_training_patients) {
      int actual_patient_id = training_patients[idx];
      
      if (training_confirmed_response_censored[idx]) { // Unclassified
        target += log_mix(prob_non_response[actual_patient_id], training_patient_response_lp[idx, 1], training_patient_response_lp[idx, 2]);
      } else {
        target += training_patient_response_lp[idx, confirmed_response_cause[actual_patient_id]];
      }
    }
  }
}

generated quantities {
 
   
  // vector<upper = 1>[n_pfs_testing_patients] oos_log_lik = zeros_vector(n_pfs_testing_patients);
  
  /*
  vector<upper = 1>[n_pfs_testing_patients] oos_log_lik = zeros_vector(n_pfs_testing_patients);
  
  {
    array[n_confresp_testing_patients] int confresp_testing_patient_ids = training_patients[confresp_testing_patients_training_idx];
    array[n_pfs_testing_patients] int pfs_testing_patient_ids = training_patients[pfs_testing_patients_training_idx];
    array[n_confresp_testing_patients] int confresp_testing_first_visit = after_cutoff_first_visit[confresp_testing_patient_ids];
    array[n_pfs_testing_patients] int pfs_testing_first_visit = after_cutoff_first_visit[pfs_testing_patient_ids];
    
    matrix[n_pfs_testing_patients, n_causes] testing_patient_response_lp; 
    
    testing_patient_response_lp[, 1] =
      calc_pch_loglik(
        pfs[pfs_testing_patient_ids], right_censored[pfs_testing_patient_ids], interval_censored[pfs_testing_patient_ids],
        pfs_ignore_interval_censoring, log_cond_prob_surv[1, pfs_testing_patient_ids], pfs_testing_first_visit);
    testing_patient_response_lp[, 2] =
      calc_pch_loglik(
        pfs[pfs_testing_patient_ids], right_censored[pfs_testing_patient_ids], interval_censored[pfs_testing_patient_ids],
        pfs_ignore_interval_censoring, log_cond_prob_surv[2, pfs_testing_patient_ids], pfs_testing_first_visit);
      
    vector[n_training_patients] temp_log_lik = zeros_vector(n_training_patients);
    
    temp_log_lik[confresp_testing_patients_training_idx] += calc_pch_loglik(
      last_unclassified_response_week[confresp_testing_patient_ids], confirmed_response_cause[confresp_testing_patient_ids], 
      early_confirmed_response_censored[confresp_testing_patient_ids], confirmed_response_interval_censored[confresp_testing_patient_ids], 0, 
      log_crcr_cond_prob_surv[, confresp_testing_patient_ids], confresp_testing_first_visit 
    );
    
    for (idx in 1:n_pfs_testing_patients) {
      int actual_patient_id = pfs_testing_patient_ids[idx];
      
      if (confirmed_response_censored[actual_patient_id]) { // Unclassified
        temp_log_lik[pfs_testing_patients_training_idx[idx]] += 
          log_mix(prob_non_response[actual_patient_id], testing_patient_response_lp[idx, 1], testing_patient_response_lp[idx, 2]);
      } else {
        temp_log_lik[pfs_testing_patients_training_idx[idx]] += testing_patient_response_lp[idx, confirmed_response_cause[actual_patient_id]];
      }
    }
    
    oos_log_lik = temp_log_lik[pfs_testing_patients_training_idx];
  }
  */
}
