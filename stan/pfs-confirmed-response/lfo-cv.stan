functions {
  #include "functions.stan"
  #include "../pos.stan"
  #include "../lfo.stan"
  #include "../gp.stan"
}

data {
  #include "data.stan"
  #include "../tumor/fine_tumor_data.stan"
  
  int<lower = 0, upper = 1> train_beyond_cutoff;
  
  int<lower = 1> n_cutoffs;
  array[n_cutoffs] int<lower = 1> cutoff_calendar_day;
}

transformed data {
  #include "transformed_data.stan"
  
  array[n_patients] int<lower = min(t_day_measure), upper = max(t_day_measure)> cutoff_last_visit_day;
  array[n_patients] int<lower = min(t_measure), upper = max(t_measure)> cutoff_last_visit_week;
  array[n_patients] int<lower = 1> last_visit_calendar_day; // Overall last calendar date of the last visit
  
  (cutoff_last_visit_day, cutoff_last_visit_week, last_visit_calendar_day) = fine_cutoff_visits(
    cutoff_calendar_day[1], calendar_day, t_measure, t_day_measure, patient_tumor_measure_pos
  );
    
  // Testing metadata: details needed to calculate the log likelihood for each cutoff date. Our out-of-sample observations are the weeks observed beyond the cutoff
  // dates. 
   
  // Get a list of patient IDs in the order of the calendar date of their last visit.  
  array[n_patients] int<lower = 1, upper = n_patients> last_visit_calendar_day_sort_idx = sort_indices_asc(last_visit_calendar_day);
   
  // For each cutoff we get the position in the above sort_idx array of the first patient to include for testing. All successive patients in that sort_idx list
  // would have visits after so should also be in the testing frame. This _idx array should have increasing values as we can use fewer and fewer patients for testing
  // as the cutoff date increases.
  array[n_cutoffs] int<lower = 1, upper = n_patients> pfs_testing_patient_idx = get_oos_patients_idx(last_visit_calendar_day[last_visit_calendar_day_sort_idx], cutoff_calendar_day);
  
  int<lower = 0, upper = n_patients> n_all_testing_patients = n_patients - pfs_testing_patient_idx[1] + 1;
  array[n_all_testing_patients] int<lower = 1, upper = n_patients> all_testing_patients = last_visit_calendar_day_sort_idx[pfs_testing_patient_idx[1]:];
   
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
      ones_int_array(n_patients), train_beyond_cutoff ? rep_array(max_confresp_week, n_patients) : cutoff_last_visit_week
    );
    
    matrix[n_patients, no_prop_hazard || pfs_only ? 1 : n_causes] training_patient_response_lp;

    training_patient_response_lp[, 1] = calc_pch_loglik(
      pfs, right_censored, interval_censored, pfs_ignore_interval_censoring, log_cond_prob_surv[1],
      ones_int_array(n_patients), train_beyond_cutoff ? rep_array(max_all_t, n_patients) : cutoff_last_visit_week
    );

    if (no_prop_hazard || pfs_only) {
      target += sum(training_patient_response_lp[, 1]);
    } else {
      training_patient_response_lp[, 2] = calc_pch_loglik(
        pfs, right_censored, interval_censored, pfs_ignore_interval_censoring, log_cond_prob_surv[2],
        ones_int_array(n_patients), train_beyond_cutoff ? rep_array(max_all_t, n_patients) : cutoff_last_visit_week 
      );
      
      for (i in 1:n_patients) {
        if (confirmed_response_censored[i] || confirmed_response_day[i] > cutoff_last_visit_day[i]) { // Unclassified
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
  array[n_cutoffs, n_cutoffs] vector[n_all_testing_patients] patient_log_lik, patient_pfs_log_lik, patient_crcr_log_lik;
  
  // matrix[n_cutoffs, n_cutoffs] oos_log_lik = rep_matrix(0, n_cutoffs, n_cutoffs); 
  // matrix[n_cutoffs, n_cutoffs] oos_pfs_log_lik = rep_matrix(0, n_cutoffs, n_cutoffs); 
  // matrix[n_cutoffs, n_cutoffs] oos_crcr_log_lik = rep_matrix(0, n_cutoffs, n_cutoffs); 
  
  for (n in 1:n_cutoffs) {
    int n_curr_patients = n_patients - pfs_testing_patient_idx[n] + 1; // How many patients after the current patient index
    int curr_first_testing_patient_idx = n_all_testing_patients - n_curr_patients + 1; 
    array[n_curr_patients] int curr_patients = last_visit_calendar_day_sort_idx[pfs_testing_patient_idx[n]:]; // Who are these patients
    array[n_curr_patients] int testing_start_week = oos_patient_first_testing_visit_week[n, curr_patients]; // Which intervals do we start from
    
    for (m in 1:n_cutoffs) {
      patient_log_lik[n, m] = zeros_vector(n_all_testing_patients);
      patient_pfs_log_lik[n, m] = zeros_vector(n_all_testing_patients);
      patient_crcr_log_lik[n, m] = zeros_vector(n_all_testing_patients);
    }
    
    for (m in n:n_cutoffs) {
      array[n_curr_patients] int testing_end_week = m < n_cutoffs ? oos_patient_last_testing_visit_week[n, m + 1, curr_patients] : rep_array(max_all_t, n_curr_patients);
    
      vector[n_curr_patients] curr_log_lik = zeros_vector(n_curr_patients); 
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
            curr_log_lik[i_idx] += 
              log_sum_exp(log_cif[1, i, max_confresp_week] + testing_patient_response_lp[i_idx, 1], log_cif[2, i, max_confresp_week] + testing_patient_response_lp[i_idx, 2]) -
              log_sum_exp(log_cif[1, i, max_confresp_week], log_cif[2, i, max_confresp_week]);
          } else {
            curr_log_lik[i_idx] += testing_patient_response_lp[i_idx, confirmed_response_cause[i]];
          }
        }
      }
      
      patient_pfs_log_lik[n, m, curr_first_testing_patient_idx:] = curr_log_lik;
      // oos_pfs_log_lik[n, m] = sum(curr_log_lik);
      
      // Get the confirmed response log likelihoods for the testing frame. 
      array[n_curr_patients] int curr_confirmed_response_calendar_day = confirmed_response_calendar_day[curr_patients]; 
      array[n_curr_patients] int confirmed_response_calendar_day_sort_idx = sort_indices_asc(curr_confirmed_response_calendar_day);
      array[n_curr_patients] int curr_sorted_confirmed_response_calendar_day = curr_confirmed_response_calendar_day[confirmed_response_calendar_day_sort_idx];
      
      int found_conf_resp_from = 0, conf_resp_from = 1;
      
      while (!found_conf_resp_from && conf_resp_from <= n_curr_patients) {
        // We need to find which patient is the first to have their confirmed response classification after the cutoff day. All following patients
        // in the _sort_idx array should also be included
        if (curr_sorted_confirmed_response_calendar_day[conf_resp_from] > cutoff_calendar_day[n]) {
          found_conf_resp_from = 1;
        } else {
          conf_resp_from += 1;
        }
      }
     
      if (found_conf_resp_from) { // We could end up with none found if for these patients their PFS is after cutoff but their confirmed response is observed before.
        int n_curr_conf_resp_patients = n_curr_patients - conf_resp_from + 1;
        array[n_curr_conf_resp_patients] int curr_cutoff_patients_idx = confirmed_response_calendar_day_sort_idx[conf_resp_from:];
        array[n_curr_conf_resp_patients] int curr_conf_resp_patients = curr_patients[curr_cutoff_patients_idx];
        vector[n_curr_conf_resp_patients] curr_crcr_log_lik = calc_pch_loglik(
          last_unclassified_response_week[curr_conf_resp_patients], confirmed_response_cause[curr_conf_resp_patients],
          confirmed_response_censored[curr_conf_resp_patients], confirmed_response_interval_censored[curr_conf_resp_patients], 0,
          log_crcr_cond_prob_surv[, curr_conf_resp_patients], testing_start_week[curr_cutoff_patients_idx], testing_end_week[curr_cutoff_patients_idx]
        );
 
        // patient_crcr_log_lik[n, m, curr_first_testing_patient_idx:][curr_cutoff_patients_idx] = curr_crcr_log_lik; 
        patient_crcr_log_lik[n, m, (curr_first_testing_patient_idx + conf_resp_from - 1):] = curr_crcr_log_lik; 
        curr_log_lik[curr_cutoff_patients_idx] += curr_crcr_log_lik;
        // oos_crcr_log_lik[n, m] = sum(curr_crcr_log_lik);
      }
      
      patient_log_lik[n, m, curr_first_testing_patient_idx:] = curr_log_lik; 
      // oos_log_lik[n, m] = sum(curr_log_lik);
    }
  }
}
