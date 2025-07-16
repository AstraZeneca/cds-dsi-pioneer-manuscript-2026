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
}

generated quantities {
  array[N_cases, max_n_patients] int last_visit_day = rep_array(0, N_cases, max_n_patients);
  array[N_cases, max_n_patients] int last_visit_week = rep_array(0, N_cases, max_n_patients);
  array[N_cases, max_n_patients] int last_visit_calendar_day = rep_array(0, N_cases, max_n_patients);

  for (case in 1:N_cases) {
    int np = n_patients[case];
    int nv = n_visits[case];
    array[np] int lcd;
    array[np] int lwk;
    array[np] int lcal;
    (lcd, lwk, lcal) = cutoff_visits(
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
  // --- get_oos_patients_idx block ---
  array[n_oos_cases, max_n_cutoffs_oos] int oos_idx = rep_array(0, n_oos_cases, max_n_cutoffs_oos);
  for (case in 1:n_oos_cases) {
    array[n_cutoffs_oos[case]] int idx = get_oos_patients_idx(
      sorted_last_visit_calendar_day_oos[case, 1:n_patients_oos[case]],
      cutoff_calendar_day_oos[case, 1:n_cutoffs_oos[case]]
    );
    for (j in 1:n_cutoffs_oos[case]) {
      oos_idx[case, j] = idx[j];
    }
    for (j in (n_cutoffs_oos[case]+1):max_n_cutoffs_oos) {
      oos_idx[case, j] = 0;
    }
  }
}
