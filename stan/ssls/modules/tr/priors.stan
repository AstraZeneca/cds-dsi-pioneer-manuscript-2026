// tr/priors.stan — active priors for total rate module
// Using new naming convention (legacy values mapped via aliases in transformed data):
pop_log_total_rate ~ normal(tr_loc_pop_mean, tr_loc_pop_sd);

if (enable_pop_cov_tr) {
	tr_coef_qr_pop ~ normal(0, 1);
}

tr_sd_trial_intercept ~ normal(0, tr_sd_trial_intercept_sd);
if (enable_trial_intercept_tr) {
	tr_raw_trial_intercept ~ std_normal();
}

tr_sd_patient_intercept ~ normal(0, tr_sd_patient_intercept_sd);
if (enable_patient_intercept_tr) {
	tr_raw_patient_intercept ~ std_normal();
}

if (enable_trial_cov_tr) {
	tr_sd_trial_slope ~ normal(0, tr_sd_trial_slope_sd);
	to_vector(tr_raw_trial_slope) ~ std_normal();
}

if (enable_patient_cov_tr) {
	tr_sd_patient_slope ~ normal(0, tr_sd_patient_slope_sd);
	to_vector(tr_raw_patient_slope) ~ std_normal();
}
