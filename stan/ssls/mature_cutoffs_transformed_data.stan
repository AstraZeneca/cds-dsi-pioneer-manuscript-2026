// Relocated mature_cutoffs_transformed_data.stan
array[n_patients, n_mature_cutoffs_calendar_days] int patient_relative_day_at_cutoff;
array[n_patients, n_mature_cutoffs_calendar_days] int patient_relative_week_at_cutoff;
array[n_patients, n_mature_cutoffs_calendar_days] int cutoff_observed_patients_mask;
array[n_mature_cutoffs_calendar_days] int<lower = 0> n_cutoff_observed_patients;
int truncated_max_all_t = 0; 
for (j in 1:n_mature_cutoffs_calendar_days) {
  array[n_patients] int last_visit_day;
  array[n_patients] int last_visit_week;
  array[n_patients] int last_visit_calendar_day_dummy;
  array[n_patients] int cutoff_last_visit_idx_dummy;
  (last_visit_day, last_visit_week, last_visit_calendar_day_dummy, cutoff_observed_patients_mask[, j], cutoff_last_visit_idx_dummy) = cutoff_visits(
    mature_cutoffs_calendar_days[j], calendar_day, t_patient_visits, t_patient_visits_day, patient_visit_pos);
  n_cutoff_observed_patients[j] = sum(cutoff_observed_patients_mask[, j]);
  patient_relative_day_at_cutoff[, j] = last_visit_day; patient_relative_week_at_cutoff[, j] = last_visit_week;
  truncated_max_all_t = max(truncated_max_all_t, max(last_visit_week));
}
array[sum(n_cutoff_observed_patients)] int<lower = 1, upper = n_patients> cutoff_observed_patients;
array[n_mature_cutoffs_calendar_days + 1] int<lower = 1> cutoff_observed_patients_pos = create_pos(n_cutoff_observed_patients);
for (j in 1:n_mature_cutoffs_calendar_days) {
  int cutoff_start, cutoff_end; (cutoff_start, cutoff_end) = get_pos(cutoff_observed_patients_pos, j);
  cutoff_observed_patients[cutoff_start:cutoff_end] = which(cutoff_observed_patients_mask[, j]);
}
