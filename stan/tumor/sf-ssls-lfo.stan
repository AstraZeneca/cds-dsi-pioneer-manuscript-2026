functions {
  #include "../util.stan"
  #include "../pos.stan"
  #include "../gp.stan"
  #include "../pfs_functions.stan"
  #include "../lfo.stan"
  #include "sf-ssls_functions.stan"
  #include "recist.stanfunctions"
}  

data {
  #include "../base_data.stan"
  #include "base_data.stan"

  #include "sf-ssls-hyperparam.stan"

  // --- LFO CV specific ---
  int<lower = 0, upper = 1> train_beyond_cutoff;
  int<lower = 1> n_cutoffs;
  array[n_cutoffs] int<lower = 1> cutoff_calendar_day;
} 

transformed data {
  #include "../base_transformed_data.stan"
  #include "tumor_transformed_data.stan"
  #include "sf-transformed_data.stan"
//   #include "other_events_transformed_data.stan"

  // --- LFO CV specific (visit-based) ---
  array[n_patients] int cutoff_last_visit_day;
  array[n_patients] int cutoff_last_visit_week;
  array[n_patients] int last_visit_calendar_day;
  array[n_patients] int last_visit_week_observed;
  
  // Precompute for each patient the last visit index before or at cutoff_last_visit_week
  array[n_patients] int cutoff_last_visit_idx;

  (cutoff_last_visit_day, cutoff_last_visit_week, last_visit_calendar_day, last_visit_week_observed, cutoff_last_visit_idx) = cutoff_visits(
    cutoff_calendar_day[1], calendar_day, t_patient_visits, t_patient_visits_day, patient_visit_pos
  );

  // This is an array of patient IDs (sorted by last visit calendar day)
  array[n_patients] int<lower = 1, upper = n_patients> last_visit_calendar_day_sort_idx = sort_indices_asc(last_visit_calendar_day);
  
  // Per cutoff, which index in the above sorted list of patient IDs, identifying the first patient to be included in the out-of-sample testing set
  array[n_cutoffs] int<lower = 1, upper = n_patients> testing_patient_idx = 
    get_oos_patients_idx(last_visit_calendar_day[last_visit_calendar_day_sort_idx], cutoff_calendar_day);

  print("testing_patient_idx = ", testing_patient_idx);

  assert_ascending(testing_patient_idx);

  int<lower = 0, upper = n_patients> n_all_testing_patients = n_patients - testing_patient_idx[1] + 1;

  // These are the patient IDs of all patients that are included in the out-of-sample testing set, sorted by last visit calendar day
  array[n_all_testing_patients] int<lower = 1, upper = n_patients> all_testing_patients = last_visit_calendar_day_sort_idx[testing_patient_idx[1]:];

  print("n_all_testing_patients = ", n_all_testing_patients);

  array[n_cutoffs, n_patients] int<lower = 0> oos_patient_first_testing_visit_week;
  array[n_cutoffs, n_cutoffs, n_patients] int oos_patient_last_testing_visit_week;
  array[n_cutoffs, n_patients] int testing_start_idx;
  array[n_cutoffs, n_cutoffs, n_patients] int testing_end_idx;

  (oos_patient_first_testing_visit_week, oos_patient_last_testing_visit_week, testing_start_idx, testing_end_idx) = get_testing_visit_week_bounds(
    testing_patient_idx, last_visit_calendar_day_sort_idx, cutoff_calendar_day, calendar_day, t_patient_visits, t_patient_visits_day, patient_visit_pos
  );

  for (n in 1:n_cutoffs) {
    int n_curr_patients = n_patients - testing_patient_idx[n] + 1; // How many patients after the current patient index
    array[n_curr_patients] int curr_patients = last_visit_calendar_day_sort_idx[testing_patient_idx[n]:]; // Who are these patients

    array[n_cutoffs] int m_size = rep_array(-1, n_cutoffs);
    
    for (m in 1:n_cutoffs) {
      if (m >= n) {
        m_size[m] = 0;

        for (i_idx in 1:n_curr_patients) {
          // Note: i is the original patient ID (1-based index from input data), not a sort position.
          // curr_patients contains original patient IDs that were reordered by sorting on last_visit_calendar_day
          int i = curr_patients[i_idx];

          int visit_start, visit_end;
          (visit_start, visit_end) = get_pos(patient_visit_pos, i); 

          int start_idx = testing_start_idx[n, i];
          int end_idx = m < n_cutoffs ? testing_end_idx[n, m + 1, i] : visit_end;

          if (start_idx > 0 && end_idx >= start_idx) {
            m_size[m] += 1;    
          }
        }
      }
    }

    print("Number of future visits:");
    print(n,": ", m_size);
  }

  array[n_patients + 1] int<lower = 1> testing_visit_pos = zeros_int_array(n_patients + 1); 
  array[n_patients] int n_patient_testing_visits;

  for (i in 1:n_patients) {
    int visit_start, visit_end;
    (visit_start, visit_end) = get_pos(patient_visit_pos, i); 

    int start_idx = testing_start_idx[1, i];

    n_patient_testing_visits[i] = visit_end - start_idx + 1;
  }

  testing_visit_pos = create_pos(n_patient_testing_visits);
}

parameters {
//   #include "other_events_parameters.stan"
  #include "sf-ssls-parameters.stan"
}

transformed parameters {
//   #include "other_events_transformed_parameters.stan"
  #include "sf-ssls-transformed_parameters.stan"
}

model {
//   #include "other_events_priors.stan"
  #include "sf-ssls-priors.stan"

  if (fit_tumor_data) {
    // --- LFO CV specific ---
    for (i in 1:n_patients) {
      if (cutoff_last_visit_idx[i] > 0) {
        int visit_start, visit_end;
        (visit_start, visit_end) = get_pos(patient_visit_pos, i);

        int cutoff_idx = cutoff_last_visit_idx[i];

        normalized_sld[visit_start:cutoff_idx] ~ sf_log_space_obs(states[visit_start:cutoff_idx], measure_sd, log_lod - log(sum_tumor_size[visit_start]));
      }
    }

    //   matrix[n_patients, n_causes] patient_response_lp = rep_matrix(0, n_patients, n_causes);

    //   // Non-target progression
    //   patient_response_lp[, 1] = calc_pch_loglik(
    //     non_target_pfs, 
    //     non_target_right_censored, 
    //     zeros_int_array(n_train_patients),
    //     0, 
    //     log_cond_prob_surv[1],
    //     ones_int_array(n_patients), train_beyond_cutoff ? rep_array(max_all_t, n_patients) : cutoff_last_visit_week
    //   );

    //   for (k in 1:n_causes) { 
    //     target += sum(patient_response_lp[, k]);
    //   }
  }
}

generated quantities {
  // Reminder to self: log_lik can be positive; probability densities aren't restricted to [0, 1]
  array[n_cutoffs, n_cutoffs] vector[n_all_testing_patients] patient_log_lik;

  array[sum(n_patient_testing_visits)] int<lower = CR, upper = PD + 1> oos_recist = rep_array(PD + 1, sum(n_patient_testing_visits));
  array[n_cutoffs, n_cutoffs] matrix<lower = 0>[PD, PD] oos_recist_confusion_matrix; // rows = observed, cols = predicted

  for (i in last_visit_calendar_day_sort_idx[testing_patient_idx[1]:]) {
    int visit_start, visit_screening_end, visit_treat_pos, visit_end;
    (visit_start, visit_screening_end, visit_treat_pos, visit_end) = get_visit_pos(patient_visit_pos, i, n_patient_screening_visits[i]);

    int cutoff_idx = cutoff_last_visit_idx[i];
    cutoff_idx = cutoff_idx > 0 ? cutoff_idx : visit_start;

    int visit_size = cutoff_idx - visit_start + 1;
    int treat_visit_size = max(0, cutoff_idx - visit_treat_pos + 1);
    int start_idx = testing_start_idx[1, i];
    int n_oos_visits = visit_end - start_idx + 1; 
    
    array[n_oos_visits + 1] int forecast_time = get_int_sub_array(t_patient_visits, patient_visit_pos, i)[visit_size:];      
    matrix[n_oos_visits, 2] forecast_patient_states;
    vector[visit_size] rep_patient_log_sld;
    vector[n_oos_visits] forecast_patient_log_sld;

    (forecast_patient_states, rep_patient_log_sld, forecast_patient_log_sld) = 
      generate_patient_states_rng(
        states[visit_start:cutoff_idx],
        forecast_time,
        patient_log_decrease_rate[i], patient_log_growth_rate[i],
        sum_tumor_size[cutoff_idx], 
        0.0001, 0.0001, // exp(patient_log_growth_lag[train_idx]), exp(pop_log_growth_transition_rate),
        rep_matrix(0.0, n_oos_visits, 2), // Hardcode zeros for forecast process noise
        measure_sd
      );

    array[treat_visit_size + n_oos_visits] int full_predict_recist = calculate_target_recist(
      exp(append_row(rep_patient_log_sld, forecast_patient_log_sld)) * 10,
      n_patient_screening_visits[i]
    );

    int oos_recist_start, oos_recist_end;
    (oos_recist_start, oos_recist_end) = get_pos(testing_visit_pos, i);

    oos_recist[oos_recist_start:oos_recist_end] = full_predict_recist[(treat_visit_size + 1):];
  }

  for (n in 1:n_cutoffs) {
    int n_curr_patients = n_patients - testing_patient_idx[n] + 1; // How many patients after the current patient index
    int curr_first_testing_patient_idx = n_all_testing_patients - n_curr_patients + 1; 
    array[n_curr_patients] int curr_patients = last_visit_calendar_day_sort_idx[testing_patient_idx[n]:]; // Who are these patients
    
    for (m in 1:n_cutoffs) {
      patient_log_lik[n, m] = zeros_vector(n_all_testing_patients);
      oos_recist_confusion_matrix[n, m] = rep_matrix(0, PD, PD);
      oos_recist_confusion_matrix[n, m] = rep_matrix(0, PD, PD);
    
      if (m >= n) {
        patient_log_lik[n, m] = zeros_vector(n_all_testing_patients);

        for (i_idx in 1:n_curr_patients) {
          // Note: i is the original patient ID (1-based index from input data), not a sort position.
          // curr_patients contains original patient IDs that were reordered by sorting on last_visit_calendar_day
          int i = curr_patients[i_idx];

          int visit_start, visit_end, first_start_idx = testing_start_idx[1, i];
          (visit_start, visit_end) = get_pos(patient_visit_pos, i); 

          int start_idx = testing_start_idx[n, i];
          int end_idx = m < n_cutoffs ? testing_end_idx[n, m + 1, i] : visit_end;
         
          // This does not exclude patients with post cutoff visits but no training visits (patients who aren't even in the study at the cutoff).
          if (start_idx > 0 && end_idx >= start_idx) {
            patient_log_lik[n, m, curr_first_testing_patient_idx + i_idx - 1] += sf_log_space_obs_lpdf(
                normalized_sld[start_idx:end_idx] | states[start_idx:end_idx], measure_sd, log_lod - log(sum_tumor_size[visit_start]));

            int oos_recist_start, oos_recist_end;
            (oos_recist_start, oos_recist_end) = get_pos(testing_visit_pos, i);

            int n_curr_testing_visits = end_idx - start_idx + 1;
            int test_start_offset = start_idx - first_start_idx;

            // print(
            //   "n: ", n, " | m: ", m, " | i: ", i,
            //   " | oos_recist_start: ", oos_recist_start,
            //   " | oos_recist_end: ", oos_recist_end,
            //   " | n_curr_testing_visits: ", n_curr_testing_visits,
            //   " | start_idx: ", start_idx,
            //   " | end_idx: ", end_idx
            // );

            // assert_equal(oos_recist_start + n_curr_testing_visits - 1, oos_recist_end);

            for (t_idx in 1:n_curr_testing_visits) {
              oos_recist_confusion_matrix[recist[start_idx + t_idx - 1], oos_recist[oos_recist_start + test_start_offset + t_idx - 1]] += 1;
            }
          }
        }
      } else {
        patient_log_lik[n, m] = rep_vector(negative_infinity(), n_all_testing_patients);
      }
    }
  } 
}
