functions {
  #include "functions.stan"
 
  tuple(array[] int, array[] int, array[] int, array[] int) cutoff_visits(
    int cutoff_calendar_day, array[] int patient_calendar_day, array[] int t_measure, array[] int t_day_measure, array[] int patient_tumor_measure_pos
  ) {
    int n_patients = size(patient_calendar_day);
    
    // 0 is the default sentinel value if last visit is negative 
    array[n_patients] int last_visit = rep_array(0, n_patients); 
    array[n_patients] int first_testing_visit_week = rep_array(0, n_patients);
    array[n_patients] int first_testing_visit_day = rep_array(0, n_patients);
    
    array[n_patients] int last_visit_calendar_day;
    
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
        last_visit[i] = patient_t_measure[patient_measure_t_sort_idx[t_idx]]; 
      }
      
      if (t_idx < n_patient_measures) {
        first_testing_visit_week[i] = patient_t_measure[patient_measure_t_sort_idx[t_idx + 1]];
        first_testing_visit_day[i] = patient_t_day_measure[patient_measure_t_sort_idx[t_idx + 1]];
      }
      
      last_visit_calendar_day[i] = patient_calendar_day[i] + patient_t_day_measure[patient_measure_t_sort_idx[n_patient_measures]] - 1;
    }
    
    return (last_visit, first_testing_visit_week, first_testing_visit_day, last_visit_calendar_day);
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

  /** Get the indices within the array of sorted last visit that will be used in all the LFO cuts.
   *
   * Each such index will indicate the first patient (in the sorted array) to be in each cut. In each future cut, the patients are a subset of the previous cut's
   * patients.
   */
  array[] int get_oos_patients_idx(array[] int sorted_last_visit_calendar_day, array[] int cutoff_calendar_day) {
    int n_cutoffs = size(cutoff_calendar_day);
    int n_patients = size(sorted_last_visit_calendar_day); 
    
    array[n_cutoffs] int patient_idx;
    array[n_cutoffs, n_patients] int testing_first_visit_week = rep_array(0, n_cutoffs, n_patients);
    
    int n_remaining_testing_patients = n_patients; 
    int curr_patient_idx = 1;
    int patient_idx_pos = 1;
    
    while (curr_patient_idx <= n_patients && patient_idx_pos <= n_cutoffs) {
      while (curr_patient_idx <= n_patients && sorted_last_visit_calendar_day[curr_patient_idx] <= cutoff_calendar_day[patient_idx_pos]) {
        curr_patient_idx += 1;
      }
      
      if (curr_patient_idx <= n_patients) {
        patient_idx[patient_idx_pos] = curr_patient_idx;
        patient_idx_pos += 1;
      }
    }
    
    return patient_idx; 
  } 
  
  array[,] int get_first_testing_visit_week(
    array[] int oos_patient_idx, array[] int last_visit_calendar_day_sort_idx,
    array[] int cutoff_calendar_day, array[] int patient_calendar_day,
    array[] int t_measure, array[] int t_day_measure, array[] int patient_tumor_measure_pos
  ) {
    int n_patients = size(patient_calendar_day);
    int n_futures = size(oos_patient_idx);
    
    array[n_futures, n_patients] int first_testing_visit_week = rep_array(0, n_futures, n_patients);
    
    for (n in 1:n_futures) {
      int n_curr_patients = n_patients - oos_patient_idx[n] + 1;
      array[n_curr_patients] int curr_patients = last_visit_calendar_day_sort_idx[oos_patient_idx[n]:];
      array[n_curr_patients] int curr_cutoff_visit_days = calendar_date_to_study_date(patient_calendar_day[curr_patients], cutoff_calendar_day[n]);
      
      for (i_idx in 1:n_curr_patients) {
        int i = curr_patients[i_idx];
        int t_measure_pos = patient_tumor_measure_pos[i]; 
        int t_measure_end = patient_tumor_measure_pos[i + 1] - 1; 
        int n_patient_measures = t_measure_end - t_measure_pos + 1;
        
        array[n_patient_measures] int patient_t_measure = t_measure[t_measure_pos:t_measure_end];
        array[n_patient_measures] int patient_t_day_measure = t_day_measure[t_measure_pos:t_measure_end];
        array[n_patient_measures] int patient_measure_t_sort_idx = sort_indices_asc(patient_t_day_measure);
        int t_idx = 1;
        int curr_visit_week = patient_t_day_measure[patient_measure_t_sort_idx[t_idx]];
        
        while (t_idx <= n_patient_measures && (curr_visit_week <= max(0, curr_cutoff_visit_days[i_idx]))) {
          t_idx += 1;
          curr_visit_week = patient_t_day_measure[patient_measure_t_sort_idx[t_idx]];
        }
        
        if (t_idx <= n_patient_measures) {
          first_testing_visit_week[n, i] = patient_t_measure[patient_measure_t_sort_idx[t_idx]]; 
        } else {
          // This would only happen if we have a patient with only baseline visits.
          fatal_error("Unexpectedly could not find the first testing visit.");
        }
      }
    }
    
    return first_testing_visit_week;
  }
}

data {
  #include "data.stan"
  
  int<lower = 1> n_cutoffs;
  array[n_cutoffs] int<lower = 1> cutoff_calendar_day;
}

transformed data {
  #include "transformed_data.stan"
  
  // Training metadata: details about who to train on and with which intervals we observed for them, given the cutoff date.
  
  array[n_patients] int<lower = min(t_measure), upper = max(t_measure)> cutoff_last_visit_week; // Last visit before cutoff
  array[n_patients] int<lower = min(t_measure), upper = max(t_measure)> after_cutoff_first_visit_week; // The first visit week after cutoff
  array[n_patients] int<lower = min(t_day_measure), upper = max(t_day_measure)> after_cutoff_first_visit_day; // The first visit day after cutoff
  array[n_patients] int<lower = 1> last_visit_calendar_day; // Overall last calendar date of the last visit
  
  (cutoff_last_visit_week, after_cutoff_first_visit_week, after_cutoff_first_visit_day, last_visit_calendar_day) = 
    cutoff_visits(cutoff_calendar_day[1], calendar_day, t_measure, t_day_measure, patient_tumor_measure_pos);
    
  int<lower = 0, upper = n_patients> n_training_patients = 0;
  
  for (i in 1:n_patients) {
    if (cutoff_last_visit_week[i] > 0) {
      n_training_patients += 1;
    }
  }
  
  print("n_training_patients = ", n_training_patients);
 
  // Who are the patients we will train with 
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
 
  // We need to adjust the data we train on to account for the cutoff: right censoring when needed. We do this for both confirmed response and PFS.
  
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
    
  // Testing metadata: details needed to calculate the log likelihood for each cutoff date. Our out-of-sample observations are the weeks observed beyond the cutoff
  // dates. 
   
  // Get a list of patient IDs in the order of the calendar date of their last visit.  
  array[n_patients] int<lower = 1, upper = n_patients> last_visit_calendar_day_sort_idx = sort_indices_asc(last_visit_calendar_day);
   
  // For each cutoff we get the position in the above sort_idx array of the first patient to include for testing. All successive patients in that sort_idx list
  // would have visits after so should also be in the testing frame. This _idx array should have increasing values as we can use fewer and fewer patients for testing
  // as the cutoff date increases.
  array[n_cutoffs] int<lower = 1, upper = n_patients> pfs_testing_patient_idx = get_oos_patients_idx(last_visit_calendar_day[last_visit_calendar_day_sort_idx], cutoff_calendar_day);
   
  // The patient-specific week to start using for log likelihood calculation 
  array[n_cutoffs, n_patients] int<lower = 0, upper = n_patients> oos_patient_first_testing_visit_week =
    get_first_testing_visit_week(
      pfs_testing_patient_idx, last_visit_calendar_day_sort_idx, cutoff_calendar_day, calendar_day, t_measure, t_day_measure, patient_tumor_measure_pos
    );
    
  array[n_patients] int confirmed_response_calendar_day = study_date_to_calendar_date(calendar_day, confirmed_response_day);
}

parameters {
  #include "parameters.stan"
}

transformed parameters {
  #include "transformed_parameters.stan"
  
  // Just like we normally do, but now focused on the training patients only and their cutoff data.
 
  matrix[n_training_patients, no_prop_hazard || pfs_only ? 1 : n_causes] training_patient_response_lp; 
  
  training_patient_response_lp[, 1] = calc_pch_loglik(
    training_pfs, training_right_censored, training_interval_censored, pfs_ignore_interval_censoring, log_cond_prob_surv[1, training_patients]
  );
  
  if (!(no_prop_hazard || pfs_only)) {
    training_patient_response_lp[, 2] = calc_pch_loglik(
      training_pfs, training_right_censored, training_interval_censored, pfs_ignore_interval_censoring, log_cond_prob_surv[2, training_patients]
    );
  }
}

model {
  #include "priors.stan"
  
  // Likelihood
  
  if (fit_data) {
    // Confirmed response model
    
    // Just like we normally do, but only for training patient and their cutoff data.
    
    profile("crcr loglik") {
      target += reduce_sum(
        partial_sum_crcr_lupmf, training_last_unclassified_response_week, crcr_grain_size,
        confirmed_response_cause[training_patients], 
        training_confirmed_response_censored, 
        crcr_ignore_interval_censoring ? zeros_int_array(n_training_patients) : training_confirmed_response_interval_censored, 
        log_crcr_cond_prob_surv[, training_patients]
      );
    }
    
    // PFS model
    
    if (no_prop_hazard || pfs_only) {
      target += sum(training_patient_response_lp[, 1]);
    } else {
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
}

generated quantities {
  vector[n_cutoffs] oos_log_lik = rep_vector(0, n_cutoffs);
  
  for (n in 1:n_cutoffs) {
    int n_curr_patients = n_patients - pfs_testing_patient_idx[n] + 1; // How many patients after the current patient index
    array[n_curr_patients] int curr_patients = last_visit_calendar_day_sort_idx[pfs_testing_patient_idx[n]:]; // Who are these patients
    array[n_curr_patients] int testing_start_week = oos_patient_first_testing_visit_week[n, curr_patients]; // Which intervals do we start from
  
    vector[n_curr_patients] curr_log_lik = rep_vector(0, n_curr_patients); 
    matrix[n_curr_patients, no_prop_hazard || pfs_only ? 1 : n_causes] testing_patient_response_lp; 
   
    // Get the PFS log likelihoods for the testing intervals/weeks. 
    testing_patient_response_lp[, 1] =
      calc_pch_loglik(
        pfs[curr_patients], right_censored[curr_patients], interval_censored[curr_patients], pfs_ignore_interval_censoring, log_cond_prob_surv[1, curr_patients], testing_start_week
      );
      
    if (!(no_prop_hazard || pfs_only)) {
      testing_patient_response_lp[, 2] =
        calc_pch_loglik(
          pfs[curr_patients], right_censored[curr_patients], interval_censored[curr_patients], pfs_ignore_interval_censoring, log_cond_prob_surv[2, curr_patients], testing_start_week
        );
    }
   
    // Get the confirmed response log likelihoods for the testing frame. 
    array[n_curr_patients] int curr_confirmed_response_calendar_day = confirmed_response_calendar_day[curr_patients]; 
    array[n_curr_patients] int confirmed_response_calendar_day_sort_idx = sort_indices_asc(curr_confirmed_response_calendar_day);
    array[n_curr_patients] int curr_sorted_confirmed_response_calendar_day = curr_confirmed_response_calendar_day[confirmed_response_calendar_day_sort_idx];
    int found_conf_resp_from = 0;
    int conf_resp_from = 1;
    
    if (no_prop_hazard || pfs_only) {
      curr_log_lik += testing_patient_response_lp[, 1];
    } else {
      for (i_idx in 1:n_curr_patients) {
        int i = curr_patients[i_idx]; // This is the actual ID of the patient, i.e, their position in the full data.
  
        if (confirmed_response_censored[i]) { // Unclassified
          curr_log_lik[i_idx] += log_mix(prob_non_response[i], testing_patient_response_lp[i_idx, 1], testing_patient_response_lp[i_idx, 2]);
        } else {
          curr_log_lik[i_idx] += testing_patient_response_lp[i_idx, confirmed_response_cause[i]];
        }
       
        // We need to find which patient is the first to have their confirmed response classification after the cutoff day. All following patients
        // in the _sort_idx array should also be included
        if (!found_conf_resp_from) { 
          if (curr_sorted_confirmed_response_calendar_day[i_idx] > cutoff_calendar_day[n]) {
            found_conf_resp_from = 1;
          } else {
            conf_resp_from += 1;
          }
        }
      }
    }
   
    if (found_conf_resp_from) { // We could end up with none found if for these patients their PFS is after cutoff but their confirmed response is observed before.
      int n_curr_conf_resp_patients = n_curr_patients - conf_resp_from + 1;
      array[n_curr_conf_resp_patients] int curr_cutoff_patients_idx = confirmed_response_calendar_day_sort_idx[conf_resp_from:];
      array[n_curr_conf_resp_patients] int curr_conf_resp_patients = curr_patients[curr_cutoff_patients_idx];

      curr_log_lik[curr_cutoff_patients_idx] += calc_pch_loglik(
        last_unclassified_response_week[curr_conf_resp_patients], confirmed_response_cause[curr_conf_resp_patients],
        early_confirmed_response_censored[curr_conf_resp_patients], confirmed_response_interval_censored[curr_conf_resp_patients], 0,
        log_crcr_cond_prob_surv[, curr_conf_resp_patients], testing_start_week[curr_cutoff_patients_idx]
      );
    }
    
    oos_log_lik[n] = sum(curr_log_lik);
  }
}
