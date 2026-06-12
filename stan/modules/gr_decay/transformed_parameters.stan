// gr_decay/transformed_parameters.stan
// Per-patient decay rate kappa_i = exp(pop intercept + optional QR covariate slopes).
// Sizing: ALWAYS vector[n_forecast_patients], filled inside the flag guard and defaulting
// to zeros when off (the value is never consumed when off — warp sites take the t branch).
// Indexed by the forecast-local patient index j (matches frac_log_growth_patient, init_*).

vector[n_forecast_patients] gr_decay_log_loc_patient = zeros_vector(n_forecast_patients);
vector[n_forecast_patients] gr_decay_kappa = zeros_vector(n_forecast_patients);

if (enable_gr_decay) {
  vector[n_forecast_patients] gr_decay_linpred_pop = enable_pop_cov_gr_decay
    ? (Q_covar_design_matrix[forecast_patient_idx, :] * gr_decay_coef_qr_pop)
    : zeros_vector(n_forecast_patients);
  gr_decay_log_loc_patient = rep_vector(gr_decay_log_loc_pop[1], n_forecast_patients)
    + gr_decay_linpred_pop;          // pop-only reduction: == gr_decay_log_loc_pop[1] when covariates off
  gr_decay_kappa = exp(gr_decay_log_loc_patient);
}
