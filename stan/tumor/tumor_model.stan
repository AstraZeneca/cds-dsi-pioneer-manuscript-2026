tumor_mean ~ normal(2.8, 0.1);
// tumor_sd ~ normal(0, 1.25);
tumor_sd ~ normal(0, 2);

patient_tumor_gp_intercept_sd ~ normal(0, 0.25);

if (use_tumor_model && multilevel_patient) { 
  patient_tumor_gp_intercept_effect ~ normal(0, patient_tumor_gp_intercept_sd);
}

if (use_tumor_model && multilevel_tumor) {
  tumor_gp_intercept_sd ~ normal(0, 0.1);
  
  int tumor_pos = 1;
  
  for (i in 1:n_patients) {
    for (j in 1:n_patient_tumors[i]) {
      tumor_gp_intercept_effect ~ normal(0, tumor_gp_intercept_sd[i]);
      
      tumor_pos += 1;
    }
  }
}

pop_tumor_gp_rho ~ inv_gamma(pop_tumor_gp_rho_alpha, pop_tumor_gp_rho_beta);

if (use_tumor_model && fit_tumor_data) {
  int tumor_pos = 1;
  int t_measure_pos = 1;
  int t_missing_measure_pos = 1;

  for (i in 1:n_patients) {
    for (j in 1:n_patient_tumors[i]) {
      int t_measure_end = t_measure_pos + n_measures[tumor_pos] - 1;
      
      matrix[n_measures[tumor_pos], n_measures[tumor_pos]] L_current_tumor_vcov = 
        calc_gp_cholesky_vcov(all_measure_t[patient_t_measure_idx[t_measure_pos:t_measure_end]], pop_tumor_gp_alpha, pop_tumor_gp_rho, delta);
      
      log(tumor_size[t_measure_pos:t_measure_end]) ~ multi_normal_cholesky(rep_vector(tumor_gp_intercept[tumor_pos], n_measures[tumor_pos]), L_current_tumor_vcov);

      tumor_pos += 1;
      t_measure_pos = t_measure_end + 1;
    }
  }
}  