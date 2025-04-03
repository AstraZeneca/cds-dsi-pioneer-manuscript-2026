tumor_mean ~ normal(tumor_mean_mean, tumor_mean_sd);
tumor_measure_error_sd ~ normal(0, tumor_measure_error_sd_sd);

trial_tumor_gp_intercept_sd ~ normal(0, trial_tumor_gp_intercept_sd_sd);
raw_trial_tumor_gp_intercept_effect ~ std_normal();

patient_tumor_gp_intercept_sd ~ normal(0, patient_tumor_gp_intercept_sd_sd);
raw_patient_tumor_gp_intercept_effect ~ std_normal();

// pop_lod ~ normal(0, 0.1);

pop_tumor_gp_alpha ~ normal(0, pop_tumor_gp_alpha_sd);
pop_tumor_gp_rho ~ inv_gamma(pop_tumor_gp_rho_alpha, pop_tumor_gp_rho_beta);

trial_tumor_gp_alpha ~ normal(0, trial_tumor_gp_alpha_sd);
trial_tumor_gp_rho ~ inv_gamma(trial_tumor_gp_rho_alpha, trial_tumor_gp_rho_beta);

patient_tumor_gp_alpha ~ normal(0, patient_tumor_gp_alpha_sd);
patient_tumor_gp_rho ~ inv_gamma(patient_tumor_gp_rho_alpha, patient_tumor_gp_rho_beta);

pop_tumor_gp_eta ~ std_normal();
trial_tumor_gp_eta ~ std_normal();
patient_tumor_gp_eta ~ std_normal();

