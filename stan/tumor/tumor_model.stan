if (fit_tumor_data) {
  for (s in 1:n_trials) {
    int actual_s = min(s, n_tumor_separate_trials);
   
    int curr_patient_pos, curr_patient_end;
    (curr_patient_pos, curr_patient_end) = get_pos(trial_patient_pos, s);
    
    for (i in curr_patient_pos:curr_patient_end) {
      int n_non_measured_visits = n_patient_non_measured_tumor_visits[i];
      int n_measured_visits = n_patient_post_treat_visits[i] - n_non_measured_visits;
      
      vector[n_patient_unique_visits[i]] curr_patient_tumor_gp = get_sub_vector(patient_obs_tumor_gp, patient_unique_visits_pos, i); 
      vector[n_patient_unique_visits[i]] curr_patient_sld = get_sub_vector(post_treat_sld, post_treat_visits_pos, i);
      
      if (n_measured_visits > 0) {
        array[n_measured_visits] int measured2patient_idx = get_int_sub_array(measured2patient_visits_idx, patient_measured_tumor_visits_pos, i); 
        vector[n_measured_visits] measured_sld = curr_patient_sld[measured2patient_idx]; 
        
        if (prod(measured_sld) <= 0) {
          fatal_error("Cannot have non-positive measured SLD values.");
        }
      
        target += normal_lpdf(measured_sld | curr_patient_tumor_gp[measured2patient_idx], tumor_measure_error_sd[actual_s]);
      }
      
      if (n_non_measured_visits > 0)  {
        array[n_non_measured_visits] int non_measured2patient_idx = get_int_sub_array(non_measured2patient_visits_idx, patient_non_measured_tumor_visits_pos, i);

        if (sum(curr_patient_sld[non_measured2patient_idx]) > 0) {
          fatal_error("All SLD values should be zero for non-measured visits.");
        }

        target += normal_lcdf(zeros_vector(n_non_measured_visits) | curr_patient_tumor_gp[non_measured2patient_idx], tumor_measure_error_sd[actual_s]);
      }
    }
  }
}