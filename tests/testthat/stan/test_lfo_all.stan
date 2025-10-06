functions {
#include "util.stan"
#include "pos.stan"
#include "lfo.stan"
  // ...include any other dependencies if needed...
}


data {
  int<lower=1> N_cases;
  int<lower=1> max_n_patients;
  int<lower=1> max_n_visits;
  int<lower=1> max_n_measures;
  array[N_cases] int cutoff_calendar_day;
  array[N_cases, max_n_patients] int patient_calendar_day;
  array[N_cases, max_n_visits] int t_patient_visits;
  array[N_cases, max_n_visits] int t_patient_visits_day;
  array[N_cases, max_n_patients+1] int patient_visit_pos;
  array[N_cases] int n_patients;
  array[N_cases] int n_visits;
  array[N_cases] int n_measures;
  array[N_cases, max_n_measures] int t_measure;
  array[N_cases, max_n_measures] int t_day_measure;
  array[N_cases, max_n_patients+1] int patient_tumor_measure_pos;
  // OOS test cases
  int<lower=1> n_oos_cases;
  int<lower=1> max_n_patients_oos;
  int<lower=1> max_n_cutoffs_oos;
  array[n_oos_cases] int<lower=1> n_patients_oos;
  array[n_oos_cases] int<lower=1> n_cutoffs_oos;
  array[n_oos_cases, max_n_patients_oos] int sorted_last_visit_calendar_day_oos;
  array[n_oos_cases, max_n_cutoffs_oos] int cutoff_calendar_day_oos;
  
  // Testing visit week bounds test cases
  int<lower=0> n_testing_bounds_cases;
  int<lower=1> max_n_patients_testing;
  int<lower=1> max_n_cutoffs_testing;
  int<lower=1> max_n_visits_testing;
  array[n_testing_bounds_cases] int<lower=1> n_patients_testing;
  array[n_testing_bounds_cases] int<lower=1> n_cutoffs_testing;
  array[n_testing_bounds_cases] int<lower=1> n_visits_testing;
  array[n_testing_bounds_cases, max_n_cutoffs_testing] int oos_patient_idx_testing;
  array[n_testing_bounds_cases, max_n_patients_testing] int last_visit_calendar_day_sort_idx_testing;
  array[n_testing_bounds_cases, max_n_cutoffs_testing] int cutoff_calendar_day_testing;
  array[n_testing_bounds_cases, max_n_patients_testing] int patient_calendar_day_testing;
  array[n_testing_bounds_cases, max_n_visits_testing] int t_patient_visits_testing;
  array[n_testing_bounds_cases, max_n_visits_testing] int t_patient_visits_day_testing;
  array[n_testing_bounds_cases, max_n_patients_testing+1] int patient_visit_pos_testing;
}

generated quantities {
  array[N_cases, max_n_patients] int last_visit_day = rep_array(0, N_cases, max_n_patients);
  array[N_cases, max_n_patients] int last_visit_week = rep_array(0, N_cases, max_n_patients);
  array[N_cases, max_n_patients] int last_visit_calendar_day = rep_array(0, N_cases, max_n_patients);
  array[N_cases, max_n_patients] int cutoff_last_visit_idx = rep_array(0, N_cases, max_n_patients);

  for (case in 1:N_cases) {
    int np = n_patients[case];
    int nv = n_visits[case];
    array[np] int lcd;
    array[np] int lwk;
    array[np] int lcal;
    array[np] int lwk_obs;
    array[np] int clvi;
    (lcd, lwk, lcal, lwk_obs, clvi) = cutoff_visits(
      cutoff_calendar_day[case],
      patient_calendar_day[case, 1:np],
      t_patient_visits[case, 1:nv],
      t_patient_visits_day[case, 1:nv],
      patient_visit_pos[case, 1:(np+1)]
    );
    for (i in 1:np) {
      last_visit_day[case, i] = lcd[i];
      last_visit_week[case, i] = lwk[i];
      last_visit_calendar_day[case, i] = lcal[i];
      cutoff_last_visit_idx[case, i] = clvi[i];
    }
  }

  // --- fine_cutoff_visits block ---
  array[N_cases, max_n_patients] int fine_last_visit_day = rep_array(0, N_cases, max_n_patients);
  array[N_cases, max_n_patients] int fine_last_visit_week = rep_array(0, N_cases, max_n_patients);
  array[N_cases, max_n_patients] int fine_last_visit_calendar_day = rep_array(0, N_cases, max_n_patients);

  for (case in 1:N_cases) {
    int np = n_patients[case];
    int nm = n_measures[case];
    array[np] int lcd;
    array[np] int lwk;
    array[np] int lcal;
    (lcd, lwk, lcal) = fine_cutoff_visits(
      cutoff_calendar_day[case],
      patient_calendar_day[case, 1:np],
      t_measure[case, 1:nm],
      t_day_measure[case, 1:nm],
      patient_tumor_measure_pos[case, 1:(np+1)]
    );
    for (i in 1:np) {
      fine_last_visit_day[case, i] = lcd[i];
      fine_last_visit_week[case, i] = lwk[i];
      fine_last_visit_calendar_day[case, i] = lcal[i];
    }
  }
  // --- get_testing_visit_week_bounds block ---
  array[n_testing_bounds_cases, max_n_cutoffs_testing, max_n_patients_testing] int first_testing_visit_week_out;
  array[n_testing_bounds_cases, max_n_cutoffs_testing, max_n_cutoffs_testing, max_n_patients_testing] int last_testing_visit_week_out;
  array[n_testing_bounds_cases, max_n_cutoffs_testing, max_n_patients_testing] int testing_start_idx_out;
  array[n_testing_bounds_cases, max_n_cutoffs_testing, max_n_cutoffs_testing, max_n_patients_testing] int testing_end_idx_out;
  
  // Initialize arrays
  for (case in 1:n_testing_bounds_cases) {
    for (n in 1:max_n_cutoffs_testing) {
      for (i in 1:max_n_patients_testing) {
        first_testing_visit_week_out[case, n, i] = 0;
        testing_start_idx_out[case, n, i] = 0;
        for (m in 1:max_n_cutoffs_testing) {
          last_testing_visit_week_out[case, n, m, i] = 0;
          testing_end_idx_out[case, n, m, i] = 0;
        }
      }
    }
  }
  
  for (case in 1:n_testing_bounds_cases) {
    int np = n_patients_testing[case];
    int nc = n_cutoffs_testing[case];
    int nv = n_visits_testing[case];
    
    array[nc, np] int first_week;
    array[nc, nc, np] int last_week;
    array[nc, np] int start_idx;
    array[nc, nc, np] int end_idx;
    
    (first_week, last_week, start_idx, end_idx) = get_testing_visit_week_bounds(
      oos_patient_idx_testing[case, 1:nc],
      last_visit_calendar_day_sort_idx_testing[case, 1:np],
      cutoff_calendar_day_testing[case, 1:nc],
      patient_calendar_day_testing[case, 1:np],
      t_patient_visits_testing[case, 1:nv],
      t_patient_visits_day_testing[case, 1:nv],
      patient_visit_pos_testing[case, 1:(np+1)]
    );
    
    for (n in 1:nc) {
      for (i in 1:np) {
        first_testing_visit_week_out[case, n, i] = first_week[n, i];
        testing_start_idx_out[case, n, i] = start_idx[n, i];
        
        for (m in 1:nc) {
          last_testing_visit_week_out[case, n, m, i] = last_week[n, m, i];
          testing_end_idx_out[case, n, m, i] = end_idx[n, m, i];
        }
      }
    }
  }
}
