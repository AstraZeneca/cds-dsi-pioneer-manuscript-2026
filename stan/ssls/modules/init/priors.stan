// init/priors.stan — active priors for initial proportion module
// Use new naming via transformed data aliases
// init/priors.stan
// initial proportion module priors
init_logit_loc_pop ~ normal(init_logit_loc_pop_mean, init_logit_loc_pop_sd);

if (enable_pop_cov_init) {
  init_coef_qr_pop ~ normal(init_coef_qr_pop_mean, init_coef_qr_pop_sd);
}

// Always put priors on SD parameters so they aren't left priorless when raw effects disabled
init_sd_trial_intercept   ~ normal(0, init_sd_trial_intercept_sd);
init_sd_patient_intercept ~ normal(0, init_sd_patient_intercept_sd);

init_raw_trial_intercept ~ std_normal();
init_raw_patient_intercept ~ std_normal();

if (enable_trial_cov_init) {
  init_sd_trial_slope ~ normal(0, init_sd_trial_slope_sd);
}
to_vector(init_raw_trial_slope) ~ std_normal();

if (enable_patient_cov_init) {
  init_sd_patient_slope ~ normal(0, init_sd_patient_slope_sd);
}
to_vector(init_raw_patient_slope) ~ std_normal();
