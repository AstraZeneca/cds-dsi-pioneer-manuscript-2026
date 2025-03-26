vector[n_patients] patient_tumor_gp_intercept_effect;

for (s in 1:n_trials) {
  int patient_pos, patient_end;
  (patient_pos, patient_end) = get_pos(trial_patient_pos, s);
  
  patient_tumor_gp_intercept_effect[patient_pos:patient_end] = 
    raw_patient_tumor_gp_intercept_effect[patient_pos:patient_end] * patient_tumor_gp_intercept_sd[min(s, n_tumor_separate_trials)];
}