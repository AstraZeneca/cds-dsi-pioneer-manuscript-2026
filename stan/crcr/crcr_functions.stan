tuple(array[] matrix, matrix) calc_cif(int n_patients, matrix log_crcr_cond_prob_surv, int max_confresp_week) {
  int n_causes = cols(log_crcr_cond_prob_surv);
  array[n_patients] matrix[max_confresp_week, n_causes] cif; // cumulative incidence function
  matrix[n_patients, n_causes] prob_cause; 
  
  for (i in 1:n_patients) { 
    int patient_prob_pos = 1 + (i - 1) * max_confresp_week; 
    
    for (t in 1:max_confresp_week) {
      cif[i, t] = 
        exp(sum(log_crcr_cond_prob_surv[patient_prob_pos:(patient_prob_pos + t - 2)]) + log1m_exp(log_crcr_cond_prob_surv[patient_prob_pos + t - 1])); 
    }
  
    for (k in 1:n_causes) {
      cif[i, , k] = cumulative_sum(cif[i, , k]);
    }
    
    prob_cause[i] = cif[i, max_confresp_week];
    prob_cause[i] /= sum(prob_cause[i]); 
  }
  
  return(cif, prob_cause);
}