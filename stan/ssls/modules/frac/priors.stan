// frac/priors.stan — active priors for fraction (decrease share) module
// Use new naming via transformed data aliases
// frac/priors.stan
// fraction module priors
frac_logit_loc_pop ~ normal(frac_logit_loc_pop_mean, frac_logit_loc_pop_sd);

// Population covariate effects (QR) only if enabled
if (enable_pop_cov_frac) {
  frac_coef_qr_pop ~ normal(frac_coef_qr_pop_mean, frac_coef_qr_pop_sd);
}

// Always give the SD parameters a prior (otherwise when effect disabled they become prior-less)
frac_sd_trial_intercept   ~ normal(0, frac_sd_trial_intercept_sd);
frac_sd_patient_intercept ~ normal(0, frac_sd_patient_intercept_sd);

// Raw random effects only sampled if enabled (length 0 otherwise)
frac_raw_trial_intercept ~ std_normal();
frac_raw_patient_intercept ~ std_normal();

if (n_covar > 0 && enable_trial_cov_frac) {
  frac_sd_trial_slope ~ normal(0, frac_sd_trial_slope_sd);
  to_vector(frac_raw_trial_slope) ~ std_normal();
}

if (n_covar > 0 && enable_patient_cov_frac) {
  frac_sd_patient_slope ~ normal(0, frac_sd_patient_slope_sd);
  to_vector(frac_raw_patient_slope) ~ std_normal();
}
