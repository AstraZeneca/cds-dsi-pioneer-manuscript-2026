functions {
#include "util.stanfunctions"
#include "pos.stanfunctions"
#include "lfo.stanfunctions"
}

data {
  int n_patients;
  int n_cutoffs;
  int n_visits;
  
  array[n_cutoffs] int oos_patient_idx;
  array[n_patients] int last_visit_calendar_day_sort_idx;
  array[n_cutoffs] int cutoff_calendar_day;
  array[n_patients] int patient_calendar_day;
  array[n_visits] int t_patient_visits;
  array[n_visits] int t_patient_visits_day;
  array[n_patients+1] int patient_visit_pos;
}

generated quantities {
  // Test the function with simple inputs
  array[n_cutoffs, n_patients] int first_testing_visit_week;
  array[n_cutoffs, n_cutoffs, n_patients] int last_testing_visit_week;
  array[n_cutoffs, n_patients] int testing_start_idx;
  array[n_cutoffs, n_cutoffs, n_patients] int testing_end_idx;
  
  (first_testing_visit_week, last_testing_visit_week, testing_start_idx, testing_end_idx) = get_testing_visit_week_bounds(
    oos_patient_idx,
    last_visit_calendar_day_sort_idx,
    cutoff_calendar_day,
    patient_calendar_day,
    t_patient_visits,
    t_patient_visits_day,
    patient_visit_pos
  );
}
