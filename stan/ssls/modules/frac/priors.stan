// frac/priors.stan — active priors for fraction (decrease share) module
// Use new naming via transformed data aliases
pop_decrease_frac_logit ~ normal(frac_logit_loc_pop_mean, frac_logit_loc_pop_sd);

if (enable_pop_cov_frac) {
  frac_coef_qr_pop ~ normal(0, 1);
}

if (enable_trial_intercept_frac) {
  frac_sd_trial_intercept ~ normal(0, frac_sd_trial_intercept_sd);
  frac_raw_trial_intercept ~ std_normal();
}

if (enable_patient_intercept_frac) {
  frac_sd_patient_intercept ~ normal(0, frac_sd_patient_intercept_sd);
  frac_raw_patient_intercept ~ std_normal();
}

if (enable_trial_cov_frac) {
  frac_sd_trial_slope ~ normal(0, 1);
  to_vector(frac_raw_trial_slope) ~ std_normal();
}

if (enable_patient_cov_frac) {
  frac_sd_patient_slope ~ normal(0, 1);
  to_vector(frac_raw_patient_slope) ~ std_normal();
}
