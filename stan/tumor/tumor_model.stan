if (fit_tumor_data) {
  for (s in 1:n_trials) {
    int actual_s = min(s, n_tumor_separate_trials);
    
    int trial_gp_pos, trial_gp_end;
    (trial_gp_pos, trial_gp_end) = get_pos(trial_unique_visits_pos, s);
   
    int curr_patient_pos, curr_patient_end;
    (curr_patient_pos, curr_patient_end) = get_pos(trial_patient_pos, s);
    
    for (i in curr_patient_pos:curr_patient_end) {
      int n_non_measured_visits = n_patient_non_measured_tumor_visits[i];
      int n_measured_visits = n_patient_post_treat_visits[i] - n_non_measured_visits;
      
      vector[n_patient_unique_visits[i]] curr_patient_sld = get_sub_vector(post_treat_sld, post_treat_visits_pos, i);
      vector[n_patient_unique_visits[i]] curr_patient_trial_tumor_gp = 
        trial_obs_tumor_gp[trial_gp_pos:trial_gp_end][get_int_sub_array(patient2trial_unique_visit_idx, patient_unique_visits_pos, i)];
      array[n_patient_unique_visits[i]] int curr_patient_idx = get_int_sub_array(patient_unique_visits_idx, patient_unique_visits_pos, i);
      matrix[n_patient_unique_visits[i], n_patient_unique_visits[i]] patient_K = gp_exp_quad_cov(
        all_measure_t[curr_patient_idx],
        patient_tumor_gp_alpha[actual_s], patient_tumor_gp_rho[actual_s], tumor_sd[actual_s]^2
      );
      
      array[n_measured_visits] int measured2patient_idx = get_int_sub_array(measured2patient_visits_idx, patient_measured_tumor_visits_pos, i); 
      matrix[n_measured_visits, n_measured_visits] measured_K = patient_K[measured2patient_idx, measured2patient_idx];
      vector[n_measured_visits] measured_sld = curr_patient_sld[measured2patient_idx]; 
      vector[n_measured_visits] mu_measured;
      
      if (n_measured_visits > 0) {
        mu_measured = curr_patient_trial_tumor_gp[measured2patient_idx] + patient_tumor_gp_intercept_effect[i];
        matrix[n_measured_visits, n_measured_visits] L_measured_K = cholesky_decompose(measured_K);

        if (prod(measured_sld) <= 0) {
          fatal_error("Cannot have non-positive measured SLD values.");
        }
      
        target += multi_normal_cholesky_lpdf(measured_sld | mu_measured, L_measured_K);
      }
      
      if (n_non_measured_visits > 0)  {
        array[n_non_measured_visits] int non_measured2patient_idx = get_int_sub_array(non_measured2patient_visits_idx, patient_non_measured_tumor_visits_pos, i);

        matrix[n_non_measured_visits, n_non_measured_visits] non_measured_K = patient_K[non_measured2patient_idx, non_measured2patient_idx];
        vector[n_non_measured_visits] non_measured_sld = curr_patient_sld[non_measured2patient_idx];

        if (sum(non_measured_sld) > 0) {
          fatal_error("All SLD values should be zero for non-measured visits.");
        }

        vector[n_non_measured_visits] mu_non_measured = curr_patient_trial_tumor_gp[non_measured2patient_idx] + patient_tumor_gp_intercept_effect[i];

        if (n_measured_visits > 0) {
          matrix[n_non_measured_visits, n_measured_visits] cross_K = gp_exp_quad_cov(
            all_measure_t[curr_patient_idx[non_measured2patient_idx]], all_measure_t[curr_patient_idx[measured2patient_idx]], 
            patient_tumor_gp_alpha[actual_s], patient_tumor_gp_rho[actual_s]
          ); 
          
          target += multi_normal_lcdf(zeros_vector(n_non_measured_visits) | mu_measured, mu_non_measured, measured_sld, measured_K, cross_K, non_measured_K, delta);
        } else {
          matrix[n_non_measured_visits, n_non_measured_visits] L_non_measured_K = cholesky_decompose(non_measured_K);

          target += multi_normal_cholesky_lcdf(zeros_vector(n_non_measured_visits) | mu_non_measured, L_non_measured_K);
        }
      }
    }
  }
}