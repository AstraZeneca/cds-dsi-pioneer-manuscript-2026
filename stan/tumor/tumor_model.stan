if (fit_tumor_data) {
  for (s in 1:n_trials) {
    int curr_patient_pos = trial_patient_pos[s];
    int curr_patient_end = trial_patient_pos[s + 1] - 1;
    int actual_s = min(s, n_tumor_separate_trials);
    
    for (i in curr_patient_pos:curr_patient_end) {
      int n_non_measured_visits = n_patient_non_measured_tumors[i], n_measured_visits = n_patient_visits[i] - n_non_measured_visits;
      
      array[n_non_measured_visits] int non_measured_visit_idx = get_int_sub_array(global_non_measured_tumors_idx, patient_non_measured_tumors_pos, i);
      array[n_measured_visits] int measured_visit_idx = get_int_sub_array(global_measured_tumors_idx, patient_measured_tumors_pos, i);
      
      if (!model_all_measures) {
        non_measured_visit_idx = pop_unique_visits_idx_dict[non_measured_visit_idx];
        measured_visit_idx = pop_unique_visits_idx_dict[measured_visit_idx];
      }
      
      if (n_non_measured_visits > 0)  {
        target += normal_lcdf(0 | tumor_mean[actual_s] + latent_pop_tumor_gp[actual_s, non_measured_visit_idx], tumor_sd[actual_s]);
      }
      
      if (n_measured_visits > 0) {
        vector[n_measured_visits] measured_sld = sum_tumor_size[get_int_sub_array(measured_tumors, patient_measured_tumors_pos, i)];
       
        if (prod(measured_sld) <= 0) {
          fatal_error("Cannot have non-positive SLD values.");
        }
         
        target += normal_lpdf(measured_sld | tumor_mean[actual_s] + latent_pop_tumor_gp[actual_s, measured_visit_idx], tumor_sd[actual_s]);
      }
    }
  }
}