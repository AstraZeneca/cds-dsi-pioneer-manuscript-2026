// tr/priors.stan — active priors for total rate module
// Using new naming convention:
tr_loc_pop ~ normal(tr_loc_pop_mean, tr_loc_pop_sd);

if (enable_pop_cov_tr) {
	tr_coef_qr_pop ~ normal(tr_coef_qr_pop_mean, tr_coef_qr_pop_sd);
}

// Always put priors on SDs (they exist even when raw vectors are length 0)
tr_sd_trial_intercept   ~ normal(0, tr_sd_trial_intercept_sd);
tr_sd_patient_intercept ~ normal(0, tr_sd_patient_intercept_sd);

tr_raw_trial_intercept ~ std_normal();
tr_raw_patient_intercept ~ std_normal();

if (n_covar > 0 && enable_trial_cov_tr) {
	tr_sd_trial_slope ~ normal(0, tr_sd_trial_slope_sd);
}

to_vector(tr_raw_trial_slope) ~ std_normal();

if (n_covar > 0 && enable_patient_cov_tr) {
	tr_sd_patient_slope ~ normal(0, tr_sd_patient_slope_sd);
}

to_vector(tr_raw_patient_slope) ~ std_normal();
