pop_tumor_intercept ~ normal(tumor_mean_mean, tumor_mean_sd);

trial_tumor_intercept_sd ~ normal(0, trial_tumor_intercept_sd_sd);
raw_trial_tumor_intercept_effect ~ std_normal();

patient_tumor_intercept_sd ~ normal(0, patient_tumor_intercept_sd_sd);
raw_patient_tumor_intercept_effect ~ std_normal();

// pop_lod ~ normal(0, 0.1);

pop_tumor_gp_alpha ~ normal(0, pop_tumor_gp_alpha_sd);
log_pop_tumor_gp_rho ~ normal(pop_tumor_gp_rho_meanlog, pop_tumor_gp_rho_sdlog);

log_trial_tumor_gp_rho_sd ~ normal(0, log_trial_tumor_gp_rho_sd_sd);
raw_log_trial_tumor_gp_rho_effect ~ std_normal(); 

log_patient_tumor_gp_rho_sd ~ normal(0, log_patient_tumor_gp_rho_sd_sd);
raw_log_patient_tumor_gp_rho_effect ~ std_normal(); 

patient_tumor_gp_eta ~ std_normal();

tumor_measure_error_sd ~ normal(0, tumor_measure_error_sd_sd);