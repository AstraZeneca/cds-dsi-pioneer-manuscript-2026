if (fit_tumor_data) {
  int n_curr_pop_unique_visits = separate_trial_tumor_gp ? 0 : n_pop_unique_visits;
  vector[n_curr_pop_unique_visits] pop_tumor_gp;
 
  if (!separate_trial_tumor_gp) { 
    pop_tumor_gp = calc_gp_pred(all_measure_t[pop_unique_visits_idx], tumor_mean[1], pop_tumor_gp_alpha[1], pop_tumor_gp_rho[1], delta, pop_tumor_gp_eta);
  }
  
  for (s in 1:n_trials) {
    int actual_s = min(s, n_tumor_separate_trials);
   
    vector[n_trial_unique_visits[s]] trial_tumor_gp;
    
    int curr_patient_pos, curr_patient_end;
    (curr_patient_pos, curr_patient_end) = get_pos(trial_patient_pos, s);
    
    if (separate_trial_tumor_gp) {
      trial_tumor_gp = calc_gp_pred(
        all_measure_t[get_int_sub_array(trial_unique_visits_idx, trial_unique_visits_pos, s)], 
        0, pop_tumor_gp_alpha[s], pop_tumor_gp_rho[s], delta, get_sub_vector(pop_tumor_gp_eta, trial_unique_visits_pos, s) 
      );
    } else {
      array[n_trial_unique_visits[s]] int curr_trial2pop_idx = get_int_sub_array(trial2pop_unique_visit_idx, trial_unique_visits_pos, s);
      trial_tumor_gp = pop_tumor_gp[curr_trial2pop_idx];
    }
    
    for (i in curr_patient_pos:curr_patient_end) {
      int n_non_measured_visits = n_patient_non_measured_tumor_visits[i];
      int n_measured_visits = n_patient_post_treat_visits[i] - n_non_measured_visits;
      
      vector[n_patient_unique_visits[i]] curr_patient_sld = get_sub_vector(post_treat_sld, post_treat_visits_pos, i);
      vector[n_patient_unique_visits[i]] curr_patient_trial_tumor_gp = trial_tumor_gp[get_int_sub_array(patient2trial_unique_visit_idx, patient_unique_visits_pos, i)];
      
      matrix[n_patient_unique_visits[i], n_patient_unique_visits[i]] patient_K = calc_gp_vcov(
        all_measure_t[get_int_sub_array(patient_unique_visits_idx, patient_unique_visits_pos, i)],
        patient_tumor_gp_alpha[actual_s], patient_tumor_gp_rho[actual_s], tumor_sd[actual_s]
      );
      
      if (n_non_measured_visits > 0)  {
        array[n_non_measured_visits] int non_measured_idx = get_int_sub_array(non_measured2patient_visits_idx, patient_non_measured_tumor_visits_pos, i);  

        matrix[n_non_measured_visits, n_non_measured_visits] non_measured_K = patient_K[non_measured_idx, non_measured_idx];
        matrix[n_non_measured_visits, n_non_measured_visits] L_non_measured_K = cholesky_decompose(non_measured_K);
        vector[n_non_measured_visits] non_measured_sld = curr_patient_sld[non_measured_idx]; 
        
        if (sum(non_measured_sld) > 0) {
          fatal_error("All SLD values should be zero for non-measured visits.");
        }

        target += multi_normal_cholesky_lcdf(
          zeros_vector(n_non_measured_visits) | curr_patient_trial_tumor_gp[non_measured_idx] + patient_tumor_gp_intercept_effect[i], L_non_measured_K
        );
      }
      
      if (n_measured_visits > 0) {
        array[n_measured_visits] int measured_idx = get_int_sub_array(measured2patient_visits_idx, patient_measured_tumor_visits_pos, i); 
        
        matrix[n_measured_visits, n_measured_visits] measured_K = patient_K[measured_idx, measured_idx];
        matrix[n_measured_visits, n_measured_visits] L_measured_K = cholesky_decompose(measured_K); 
        vector[n_measured_visits] measured_sld = curr_patient_sld[measured_idx]; 
       
        if (prod(measured_sld) <= 0) {
          fatal_error("Cannot have non-positive SLD values.");
        }
         
        target += multi_normal_cholesky_lpdf(measured_sld | curr_patient_trial_tumor_gp[measured_idx] + patient_tumor_gp_intercept_effect[i], L_measured_K);
      }
    }
  }
}