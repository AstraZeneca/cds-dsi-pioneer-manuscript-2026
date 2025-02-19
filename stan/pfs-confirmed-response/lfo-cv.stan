functions {
  #include "functions.stan"
 
  tuple(array[] int, array[] int, array[] int) cutoff_visits(
    int cutoff_calendar_day, array[] int patient_calendar_day, array[] int t_measure, array[] int t_day_measure, array[] int patient_tumor_measure_pos
  ) {
    int n_patients = size(patient_calendar_day);
    
    // 0 is the default sentinel value if last visit is negative 
    array[n_patients] int last_visit_day = rep_array(0, n_patients), last_visit_week = rep_array(0, n_patients); 

    array[n_patients] int last_visit_calendar_day;
    
    for (i in 1:n_patients) {
      int t_measure_pos = patient_tumor_measure_pos[i]; 
      int t_measure_end = patient_tumor_measure_pos[i + 1] - 1; 
      int n_patient_measures = t_measure_end - t_measure_pos + 1;
      
      // This is the study date for this patient that such a cutoff would have occured on
      int patient_cutoff_study_day = calendar_date_to_study_date(patient_calendar_day[i], cutoff_calendar_day); 
      
      array[n_patient_measures] int patient_t_measure = t_measure[t_measure_pos:t_measure_end];
      array[n_patient_measures] int patient_t_day_measure = t_day_measure[t_measure_pos:t_measure_end];
      array[n_patient_measures] int patient_measure_t_sort_idx = sort_indices_asc(patient_t_day_measure);
      int t_idx = 0;
      
      while (t_idx < n_patient_measures && patient_t_day_measure[patient_measure_t_sort_idx[t_idx + 1]] <= patient_cutoff_study_day) {
        t_idx += 1;
      }
      
      if (t_idx > 0) {  
        last_visit_day[i] = patient_t_day_measure[patient_measure_t_sort_idx[t_idx]]; 
        last_visit_week[i] = patient_t_measure[patient_measure_t_sort_idx[t_idx]]; 
      }
      
      last_visit_calendar_day[i] = patient_calendar_day[i] + patient_t_day_measure[patient_measure_t_sort_idx[n_patient_measures]] - 1;
    }
    
    return (last_visit_day, last_visit_week, last_visit_calendar_day);
  } 
 
  /** Get the indices within the array of sorted last visit that will be used in all the LFO cuts.
   *
   * Each such index will indicate the first patient (in the sorted array) to be in each cut. In each future cut, the patients are a subset of the previous cut's
   * patients.
   */
  array[] int get_oos_patients_idx(array[] int sorted_last_visit_calendar_day, array[] int cutoff_calendar_day) {
    int n_cutoffs = size(cutoff_calendar_day);
    int n_patients = size(sorted_last_visit_calendar_day);

    array[n_cutoffs] int patient_idx = rep_array(0, n_cutoffs);
    
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
      } else {
        print("Warning: no patients have any visits after cutoff ", patient_idx_pos);
      }
    }

    return patient_idx;
  }
  
  tuple(array[,] int, array[,,] int) get_testing_visit_week_bounds(
    array[] int oos_patient_idx, array[] int last_visit_calendar_day_sort_idx,
    array[] int cutoff_calendar_day, array[] int patient_calendar_day,
    array[] int t_measure, array[] int t_day_measure, array[] int patient_tumor_measure_pos
  ) {
    int n_patients = size(patient_calendar_day);
    int n_futures = size(oos_patient_idx);
    int min_all_t = min(t_measure);
    
    array[n_futures, n_patients] int first_testing_visit_week = rep_array(0, n_futures, n_patients);
    array[n_futures, n_futures, n_patients] int last_testing_visit_week = rep_array(min_all_t, n_futures, n_futures, n_patients);
    
    for (n in 1:n_futures) {
      int n_curr_patients = n_patients - oos_patient_idx[n] + 1;
      array[n_curr_patients] int curr_patients = last_visit_calendar_day_sort_idx[oos_patient_idx[n]:];
      array[n_curr_patients] int lower_cutoff_visit_days = calendar_date_to_study_date(patient_calendar_day[curr_patients], cutoff_calendar_day[n]);
      
      for (i_idx in 1:n_curr_patients) {
        int i = curr_patients[i_idx];
        int t_measure_pos = patient_tumor_measure_pos[i]; 
        int t_measure_end = patient_tumor_measure_pos[i + 1] - 1; 
        int n_patient_measures = t_measure_end - t_measure_pos + 1;
        
        array[n_patient_measures] int patient_t_measure = t_measure[t_measure_pos:t_measure_end];
        array[n_patient_measures] int patient_t_day_measure = t_day_measure[t_measure_pos:t_measure_end];
        array[n_patient_measures] int patient_measure_t_sort_idx = sort_indices_asc(patient_t_day_measure);
        int t_idx = 1;
        
        while (t_idx <= n_patient_measures && (patient_t_day_measure[patient_measure_t_sort_idx[t_idx]]<= max(0, lower_cutoff_visit_days[i_idx]))) {
          t_idx += 1;
        }
        
        if (t_idx <= n_patient_measures) {
          first_testing_visit_week[n, i] = patient_t_measure[patient_measure_t_sort_idx[t_idx]]; 
        } else {
          // This would only happen if we have a patient with only baseline visits.
          fatal_error("Unexpectedly could not find the first testing visit.");
        }
      }
      
      for (m in (n + 1):n_futures) {
        array[n_curr_patients] int upper_cutoff_visit_days = calendar_date_to_study_date(patient_calendar_day[curr_patients], cutoff_calendar_day[m]);
      
        for (i_idx in 1:n_curr_patients) {
          int i = curr_patients[i_idx];
          int t_measure_pos = patient_tumor_measure_pos[i]; 
          int t_measure_end = patient_tumor_measure_pos[i + 1] - 1; 
          int n_patient_measures = t_measure_end - t_measure_pos + 1;
          
          array[n_patient_measures] int patient_t_measure = t_measure[t_measure_pos:t_measure_end];
          array[n_patient_measures] int patient_t_day_measure = t_day_measure[t_measure_pos:t_measure_end];
          array[n_patient_measures] int patient_measure_t_sort_idx = sort_indices_desc(patient_t_day_measure);
          int t_idx = 1;
          
          while (t_idx <= n_patient_measures && (patient_t_day_measure[patient_measure_t_sort_idx[t_idx]] > max(0, upper_cutoff_visit_days[i_idx]))) {
            t_idx += 1;
          }
          
          if (t_idx <= n_patient_measures) {
            last_testing_visit_week[n, m, i] = patient_t_measure[patient_measure_t_sort_idx[t_idx]]; 
          } 
        } 
      }
    }
    
    return (first_testing_visit_week, last_testing_visit_week);
  }
}

data {
  #include "data.stan"
  
  int<lower = 1> n_cutoffs;
  array[n_cutoffs] int<lower = 1> cutoff_calendar_day;
}

transformed data {
  #include "transformed_data.stan"
  
  array[n_patients] int<lower = min(t_day_measure), upper = max(t_day_measure)> cutoff_last_visit_day;
  array[n_patients] int<lower = min(t_measure), upper = max(t_measure)> cutoff_last_visit_week;
  array[n_patients] int<lower = 1> last_visit_calendar_day; // Overall last calendar date of the last visit
  
  (cutoff_last_visit_day, cutoff_last_visit_week, last_visit_calendar_day) = cutoff_visits(cutoff_calendar_day[1], calendar_day, t_measure, t_day_measure, patient_tumor_measure_pos);
    
  // Testing metadata: details needed to calculate the log likelihood for each cutoff date. Our out-of-sample observations are the weeks observed beyond the cutoff
  // dates. 
   
  // Get a list of patient IDs in the order of the calendar date of their last visit.  
  array[n_patients] int<lower = 1, upper = n_patients> last_visit_calendar_day_sort_idx = sort_indices_asc(last_visit_calendar_day);
   
  // For each cutoff we get the position in the above sort_idx array of the first patient to include for testing. All successive patients in that sort_idx list
  // would have visits after so should also be in the testing frame. This _idx array should have increasing values as we can use fewer and fewer patients for testing
  // as the cutoff date increases.
  array[n_cutoffs] int<lower = 1, upper = n_patients> pfs_testing_patient_idx = get_oos_patients_idx(last_visit_calendar_day[last_visit_calendar_day_sort_idx], cutoff_calendar_day);
   
  // The patient-specific week to start using for log likelihood calculation 
  array[n_cutoffs, n_patients] int<lower = 0> oos_patient_first_testing_visit_week;
  array[n_cutoffs, n_cutoffs, n_patients] int oos_patient_last_testing_visit_week;
 
  (oos_patient_first_testing_visit_week, oos_patient_last_testing_visit_week) = get_testing_visit_week_bounds(
    pfs_testing_patient_idx, last_visit_calendar_day_sort_idx, cutoff_calendar_day, calendar_day, t_measure, t_day_measure, patient_tumor_measure_pos
  );
  
  array[n_patients] int confirmed_response_calendar_day = study_date_to_calendar_date(calendar_day, confirmed_response_day);
}

parameters {
  #include "parameters.stan"
}

transformed parameters {
  #include "transformed_parameters.stan"
}

model {
  #include "priors.stan"
  
  // Likelihood
  
  if (fit_data) {
    last_unclassified_response_week ~ pch(
      confirmed_response_cause,
      confirmed_response_censored,
      confirmed_response_interval_censored,
      crcr_ignore_interval_censoring,
      log_crcr_cond_prob_surv,
      ones_int_array(n_patients), cutoff_last_visit_week
    );
    
    matrix[n_patients, no_prop_hazard || pfs_only ? 1 : n_causes] training_patient_response_lp;

    training_patient_response_lp[, 1] = calc_pch_loglik(
      pfs, right_censored, interval_censored, pfs_ignore_interval_censoring, log_cond_prob_surv[1],
      ones_int_array(n_patients), cutoff_last_visit_week
    );

    if (no_prop_hazard || pfs_only) {
      target += sum(training_patient_response_lp[, 1]);
    } else {
      training_patient_response_lp[, 2] = calc_pch_loglik(
        pfs, right_censored, interval_censored, pfs_ignore_interval_censoring, log_cond_prob_surv[2],
        ones_int_array(n_patients), cutoff_last_visit_week
      );
      
      for (i in 1:n_patients) {
        if (confirmed_response_censored[i] || confirmed_response_day[i] > cutoff_last_visit_day[i]) { // Unclassified
          // target += log_mix(prob_non_response[i], training_patient_response_lp[i, 1], training_patient_response_lp[i, 2]);
          
          target += log_sum_exp(log_cif[1, i, max_confresp_week] + training_patient_response_lp[i, 1], log_cif[2, i, max_confresp_week] + training_patient_response_lp[i, 2]) -
            log_sum_exp(log_cif[1, i, max_confresp_week], log_cif[2, i, max_confresp_week]);
        } else {
          target += training_patient_response_lp[i, confirmed_response_cause[i]];
        }
      }
    }
  }
}

generated quantities {
  matrix[n_cutoffs, n_cutoffs] oos_log_lik = rep_matrix(0, n_cutoffs, n_cutoffs); 
  matrix[n_cutoffs, n_cutoffs] oos_pfs_log_lik = rep_matrix(0, n_cutoffs, n_cutoffs); 
  matrix[n_cutoffs, n_cutoffs] oos_crcr_log_lik = rep_matrix(0, n_cutoffs, n_cutoffs); 
  
  for (n in 1:n_cutoffs) {
    int n_curr_patients = n_patients - pfs_testing_patient_idx[n] + 1; // How many patients after the current patient index
    array[n_curr_patients] int curr_patients = last_visit_calendar_day_sort_idx[pfs_testing_patient_idx[n]:]; // Who are these patients
    array[n_curr_patients] int testing_start_week = oos_patient_first_testing_visit_week[n, curr_patients]; // Which intervals do we start from
    
    for (m in n:n_cutoffs) {
      array[n_curr_patients] int testing_end_week = m < n_cutoffs ? oos_patient_last_testing_visit_week[n, m + 1, curr_patients] : rep_array(max_all_t, n_curr_patients);
    
      vector[n_curr_patients] curr_log_lik = rep_vector(0, n_curr_patients); 
      matrix[n_curr_patients, no_prop_hazard || pfs_only ? 1 : n_causes] testing_patient_response_lp; 
      
      // Get the PFS log likelihoods for the testing intervals/weeks. 
      testing_patient_response_lp[, 1] =
        calc_pch_loglik(
          pfs[curr_patients], right_censored[curr_patients], interval_censored[curr_patients], 
          pfs_ignore_interval_censoring, 
          log_cond_prob_surv[1, curr_patients], 
          testing_start_week, testing_end_week
        );
        
      if (no_prop_hazard || pfs_only) {
        curr_log_lik += testing_patient_response_lp[, 1];
      } else {
        testing_patient_response_lp[, 2] =
          calc_pch_loglik(
            pfs[curr_patients], right_censored[curr_patients], interval_censored[curr_patients], 
            pfs_ignore_interval_censoring, 
            log_cond_prob_surv[2, curr_patients], 
            testing_start_week, testing_end_week
          );
        
        for (i_idx in 1:n_curr_patients) {
          int i = curr_patients[i_idx]; // This is the actual ID of the patient, i.e, their position in the full data.
    
          if (confirmed_response_censored[i]) { // Unclassified
            // curr_log_lik[i_idx] += log_mix(prob_non_response[i], testing_patient_response_lp[i_idx, 1], testing_patient_response_lp[i_idx, 2]);
            curr_log_lik[i_idx] += 
              log_sum_exp(log_cif[1, i, max_confresp_week] + testing_patient_response_lp[i_idx, 1], log_cif[2, i, max_confresp_week] + testing_patient_response_lp[i_idx, 2]) -
              log_sum_exp(log_cif[1, i, max_confresp_week], log_cif[2, i, max_confresp_week]);
          } else {
            curr_log_lik[i_idx] += testing_patient_response_lp[i_idx, confirmed_response_cause[i]];
          }
        }
      }
      
      // Get the confirmed response log likelihoods for the testing frame. 
      array[n_curr_patients] int curr_confirmed_response_calendar_day = confirmed_response_calendar_day[curr_patients]; 
      array[n_curr_patients] int confirmed_response_calendar_day_sort_idx = sort_indices_asc(curr_confirmed_response_calendar_day);
      array[n_curr_patients] int curr_sorted_confirmed_response_calendar_day = curr_confirmed_response_calendar_day[confirmed_response_calendar_day_sort_idx];
      
      int found_conf_resp_from = 0, conf_resp_from = 1;
      
      while (!found_conf_resp_from) {
        // We need to find which patient is the first to have their confirmed response classification after the cutoff day. All following patients
        // in the _sort_idx array should also be included
        if (curr_sorted_confirmed_response_calendar_day[conf_resp_from] > cutoff_calendar_day[n]) {
          found_conf_resp_from = 1;
        } else {
          conf_resp_from += 1;
        }
      }
      
      oos_pfs_log_lik[n, m] = sum(curr_log_lik);
     
      if (found_conf_resp_from) { // We could end up with none found if for these patients their PFS is after cutoff but their confirmed response is observed before.
        int n_curr_conf_resp_patients = n_curr_patients - conf_resp_from + 1;
        array[n_curr_conf_resp_patients] int curr_cutoff_patients_idx = confirmed_response_calendar_day_sort_idx[conf_resp_from:];
        array[n_curr_conf_resp_patients] int curr_conf_resp_patients = curr_patients[curr_cutoff_patients_idx];
        vector[n_curr_conf_resp_patients] curr_crcr_log_lik = calc_pch_loglik(
          last_unclassified_response_week[curr_conf_resp_patients], confirmed_response_cause[curr_conf_resp_patients],
          early_confirmed_response_censored[curr_conf_resp_patients], confirmed_response_interval_censored[curr_conf_resp_patients], 0,
          log_crcr_cond_prob_surv[, curr_conf_resp_patients], testing_start_week[curr_cutoff_patients_idx], testing_end_week[curr_cutoff_patients_idx]
        );
  
        curr_log_lik[curr_cutoff_patients_idx] += curr_crcr_log_lik;
        oos_crcr_log_lik[n, m] = sum(curr_crcr_log_lik);
      }
      
      oos_log_lik[n, m] = sum(curr_log_lik);
    }
  }
}
