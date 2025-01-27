functions {
  #include "functions.stan"
 
  tuple(array[] int, array[] int) cutoff_visits(int cutoff_calendar_week, array[] int patient_calendar_week, array[] int t_measure, array[] int patient_tumor_measure_pos) {
    int n_patients = size(patient_calendar_week);
    
    // 0 is the default sentinel value if last visit is negative 
    array[n_patients] int last_visit = rep_array(0, n_patients); 
    array[n_patients] int first_testing_visit = rep_array(0, n_patients);
    
    for (i in 1:n_patients) {
      int t_measure_pos = patient_tumor_measure_pos[i]; 
      int t_measure_end = patient_tumor_measure_pos[i + 1] - 1; 
      int n_patient_measures = t_measure_end - t_measure_pos + 1;
      int patient_cutoff_study_week = calendar_week_to_study_week(patient_calendar_week[i], cutoff_calendar_week); 
      
      array[n_patient_measures] int sorted_patient_measure_t = sort_asc(t_measure[t_measure_pos:t_measure_end]);
      int t_idx = 0;
      
      while (t_idx < n_patient_measures && sorted_patient_measure_t[t_idx + 1] <= patient_cutoff_study_week) {
        t_idx += 1;
      }
      
      if (t_idx > 0) {  
        last_visit[i] = sorted_patient_measure_t[t_idx]; 
      }
      
      if (t_idx < n_patient_measures) {
        first_testing_visit[i] = sorted_patient_measure_t[t_idx + 1];
      }
    }
    
    return (last_visit, first_testing_visit);
  } 
  
  tuple(array[] int, array[] int, array[] int) identify_cutoff_training_and_testing_patients(
     array[] int calendar_week, int cutoff_calendar_week, array[] int last_visit, int n_training_patients, 
    int n_confresp_testing_patients, array[] int last_unclass_week, array[] int confresp_right_censored, array[] int confresp_interval_censored, 
    int n_pfs_testing_patients, array[] int pfs, array[] int pfs_right_censored, array[] int pfs_interval_censored
  ) {
    int n_patients = size(calendar_week);
    
    array[n_training_patients] int training_patients;
    array[n_confresp_testing_patients] int confresp_testing_patients_training_idx;
    array[n_pfs_testing_patients] int pfs_testing_patients_training_idx;
    
    int training_pos = 1, pfs_testing_pos = 1, confresp_testing_pos = 1;
    
    for (i in 1:n_patients) {
      if (last_visit[i] > 0) {
        training_patients[training_pos] = i;
        
        if (calendar_week[i] + pfs[i] - pfs_right_censored[i] + pfs_interval_censored[i] > cutoff_calendar_week) {
          pfs_testing_patients_training_idx[pfs_testing_pos] = training_pos;
          pfs_testing_pos += 1;
          
          if (calendar_week[i] + last_unclass_week[i] - confresp_right_censored[i] + confresp_interval_censored[i] > cutoff_calendar_week) {
            confresp_testing_patients_training_idx[confresp_testing_pos] = training_pos;
            confresp_testing_pos += 1;
          }
        }
        
        training_pos += 1;
      }
    }
    
    return (training_patients, confresp_testing_patients_training_idx, pfs_testing_patients_training_idx);
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
}

data {
  #include "data.stan"
  
  int<lower = 1, upper = max(calendar_week)> cutoff_calendar_week;
}

transformed data {
  #include "transformed_data.stan"
  
  array[n_patients] int<lower = min(t_measure), upper = max(t_measure)> cutoff_last_visit;
  array[n_patients] int<lower = min(t_measure), upper = max(t_measure)> after_cutoff_first_visit;
  
  (cutoff_last_visit, after_cutoff_first_visit) = cutoff_visits(cutoff_calendar_week, calendar_week, t_measure, patient_tumor_measure_pos);
  
  int<lower = 1, upper = n_patients> n_pfs_testing_patients = 0;
  int<lower = 0, upper = n_pfs_testing_patients> n_confresp_testing_patients = 0;
  int<lower = n_pfs_testing_patients, upper = n_patients> n_training_patients = 0;
  
  for (i in 1:n_patients) {
    if (cutoff_last_visit[i] > 0) {
      n_training_patients += 1;
      
      if (calendar_week[i] + pfs[i] - right_censored[i] + interval_censored[i] > cutoff_calendar_week) {
        n_pfs_testing_patients += 1;
        
        if (calendar_week[i] + last_unclassified_response_week[i] - confirmed_response_censored[i] + confirmed_response_interval_censored[i] > cutoff_calendar_week) {
          n_confresp_testing_patients += 1;
        }
      }
    }
  }
  
  print("n_training_patients = ", n_training_patients, ", n_pfs_testing_patients = ", n_pfs_testing_patients, ", n_confresp_testing_patients = ", n_confresp_testing_patients);
  
  array[n_training_patients] int<lower = 1, upper = n_patients> training_patients;
  array[n_confresp_testing_patients] int<lower = 1, upper = n_training_patients> confresp_testing_patients_training_idx; // Index in training_patients 
  array[n_pfs_testing_patients] int<lower = 1, upper = n_training_patients> pfs_testing_patients_training_idx; // Index in training_patients 
  
  (training_patients, confresp_testing_patients_training_idx, pfs_testing_patients_training_idx) = 
    identify_cutoff_training_and_testing_patients(
      calendar_week, cutoff_calendar_week, cutoff_last_visit, n_training_patients, 
      n_confresp_testing_patients, last_unclassified_response_week, confirmed_response_censored, confirmed_response_interval_censored,
      n_pfs_testing_patients, pfs, right_censored, interval_censored
    );
  
  array[n_training_patients] int<lower = 0> training_last_unclassified_response_week = last_unclassified_response_week[training_patients];
  array[n_training_patients] int<lower = 0, upper = 1> training_confirmed_response_censored = confirmed_response_censored[training_patients];
  array[n_training_patients] int<lower = 0> training_confirmed_response_interval_censored = confirmed_response_interval_censored[training_patients];
  
  (training_last_unclassified_response_week, training_confirmed_response_censored, training_confirmed_response_interval_censored) =
    cutoff_surv_data(
      cutoff_last_visit[training_patients], training_last_unclassified_response_week, training_confirmed_response_censored, training_confirmed_response_interval_censored
    ); 
    
  array[n_training_patients] int<lower = 0> training_pfs = pfs[training_patients]; // How many periods after baseline did patient survive.
  array[n_training_patients] int<lower = 0, upper = 1> training_right_censored = right_censored[training_patients];
  array[n_training_patients] int<lower = 0> training_interval_censored = interval_censored[training_patients];

  (training_pfs, training_right_censored, training_interval_censored) =
    cutoff_surv_data(cutoff_last_visit[training_patients], training_pfs, training_right_censored, training_interval_censored);
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
}
