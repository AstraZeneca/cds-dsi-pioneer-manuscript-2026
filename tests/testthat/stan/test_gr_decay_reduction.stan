data {
  int<lower=1> n_forecast_patients;
  int<lower=0> n_covar;
  int<lower=0,upper=1> enable_gr_decay;
  int<lower=0,upper=1> enable_pop_cov_gr_decay;
  array[n_forecast_patients] int forecast_patient_idx;
  matrix[n_forecast_patients, n_covar] Q_covar_design_matrix;
  array[enable_gr_decay ? 1 : 0] real gr_decay_log_loc_pop;
  vector[(enable_gr_decay && enable_pop_cov_gr_decay) ? n_covar : 0] gr_decay_coef_qr_pop;
}
generated quantities {
  vector[n_forecast_patients] gr_decay_log_loc_patient = zeros_vector(n_forecast_patients);
  if (enable_gr_decay) {
    vector[n_forecast_patients] linpred = enable_pop_cov_gr_decay
      ? (Q_covar_design_matrix[forecast_patient_idx, :] * gr_decay_coef_qr_pop)
      : zeros_vector(n_forecast_patients);
    gr_decay_log_loc_patient = rep_vector(gr_decay_log_loc_pop[1], n_forecast_patients) + linpred;
  }
}
