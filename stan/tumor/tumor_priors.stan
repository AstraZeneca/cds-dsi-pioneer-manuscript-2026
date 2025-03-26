tumor_mean ~ normal(3, 5);
tumor_sd ~ normal(0, 10);

patient_tumor_gp_intercept_sd ~ normal(0, patient_tumor_gp_sd_sd);
raw_patient_tumor_gp_intercept_effect ~ std_normal();

// pop_lod ~ normal(0, 0.1);

pop_tumor_gp_alpha ~ normal(0, pop_tumor_gp_alpha_sd);
pop_tumor_gp_rho ~ inv_gamma(pop_tumor_gp_rho_alpha, pop_tumor_gp_rho_beta);

patient_tumor_gp_alpha ~ normal(0, patient_tumor_gp_alpha_sd);
patient_tumor_gp_rho ~ inv_gamma(patient_tumor_gp_rho_alpha, patient_tumor_gp_rho_beta);

pop_tumor_gp_eta ~ std_normal();
patient_tumor_gp_eta ~ std_normal();

