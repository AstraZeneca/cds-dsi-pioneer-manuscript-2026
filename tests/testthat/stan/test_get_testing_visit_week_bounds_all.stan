functions {
#include "util.stanfunctions"
#include "pos.stanfunctions"
#include "lfo.stanfunctions"
}

data {
  int<lower=1> N_cases;
  int<lower=1> max_n_patients;
  int<lower=1> max_n_cutoffs;
  int<lower=1> max_n_visits;
  
  // Test case dimensions
  array[N_cases] int<lower=1> n_patients;
  array[N_cases] int<lower=1> n_cutoffs;
  array[N_cases] int<lower=1> n_visits;
  
  // Input data for get_testing_visit_week_bounds
  array[N_cases, max_n_cutoffs] int oos_patient_idx;
  array[N_cases, max_n_patients] int last_visit_calendar_day_sort_idx;
  array[N_cases, max_n_cutoffs] int cutoff_calendar_day;
  array[N_cases, max_n_patients] int patient_calendar_day;
  array[N_cases, max_n_visits] int t_patient_visits;
  array[N_cases, max_n_visits] int t_patient_visits_day;
  array[N_cases, max_n_patients+1] int patient_visit_pos;
  
  // Error handling flags
  array[N_cases] int expect_error;
}

generated quantities {
  // Output arrays for all test cases
  array[N_cases, max_n_cutoffs, max_n_patients] int first_testing_visit_week;
  array[N_cases, max_n_cutoffs, max_n_cutoffs, max_n_patients] int last_testing_visit_week;
  array[N_cases, max_n_cutoffs, max_n_patients] int testing_start_idx;
  array[N_cases, max_n_cutoffs, max_n_cutoffs, max_n_patients] int testing_end_idx;
  array[N_cases] int case_status;
  
  // Initialize arrays with loops (needed for 4D arrays)
  for (case in 1:N_cases) {
    case_status[case] = 0; // 0 = success, 1 = error
    for (n in 1:max_n_cutoffs) {
      for (i in 1:max_n_patients) {
        first_testing_visit_week[case, n, i] = -999;
        testing_start_idx[case, n, i] = -999;
        for (m in 1:max_n_cutoffs) {
          last_testing_visit_week[case, n, m, i] = -999;
          testing_end_idx[case, n, m, i] = -999;
        }
      }
    }
  }
  
  for (case in 1:N_cases) {
    int np = n_patients[case];
    int nc = n_cutoffs[case];
    int nv = n_visits[case];
    
    // Handle potential errors by wrapping in conditional logic
    // Note: Stan doesn't have true try-catch, so we rely on fatal_error for expected errors
    if (expect_error[case] == 0) {
      // Expected to succeed - call the function
      array[nc, np] int first_week;
      array[nc, nc, np] int last_week;
      array[nc, np] int start_idx;
      array[nc, nc, np] int end_idx;
      
      (first_week, last_week, start_idx, end_idx) = get_testing_visit_week_bounds(
        oos_patient_idx[case, 1:nc],
        last_visit_calendar_day_sort_idx[case, 1:np],
        cutoff_calendar_day[case, 1:nc],
        patient_calendar_day[case, 1:np],
        t_patient_visits[case, 1:nv],
        t_patient_visits_day[case, 1:nv],
        patient_visit_pos[case, 1:(np+1)]
      );
      
      // Store results in padded output arrays
      for (n in 1:nc) {
        for (i in 1:np) {
          first_testing_visit_week[case, n, i] = first_week[n, i];
          testing_start_idx[case, n, i] = start_idx[n, i];
          
          for (m in 1:nc) {
            last_testing_visit_week[case, n, m, i] = last_week[n, m, i];
            testing_end_idx[case, n, m, i] = end_idx[n, m, i];
          }
        }
      }
      case_status[case] = 0; // Success
      
    } else {
      // Expected to error - mark as error case for validation
      case_status[case] = 1; // Expected error
      // Leave outputs as -999 (error sentinel values)
    }
  }
}
