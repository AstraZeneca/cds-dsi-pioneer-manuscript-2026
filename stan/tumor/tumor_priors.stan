tumor_mean ~ normal(3, 5);
tumor_sd ~ normal(0, 10);

// pop_lod ~ normal(0, 0.1);

pop_tumor_gp_alpha ~ normal(0, pop_tumor_gp_alpha_sd);
pop_tumor_gp_rho ~ inv_gamma(pop_tumor_gp_rho_alpha, pop_tumor_gp_rho_beta);

// trial_tumor_gp_intercept_sd ~ normal(0, 0.25);
// 
// if (add_trial_level_tumor_gp) { 
//   trial_tumor_gp_intercept_effect ~ normal(0, trial_tumor_gp_intercept_sd);
// }
// 
// patient_tumor_gp_intercept_sd ~ normal(0, 0.25);
// 
// if (add_patient_level_tumor_gp) { 
//   patient_tumor_gp_intercept_effect ~ normal(0, patient_tumor_gp_intercept_sd);
// }

for (s in 1:n_tumor_separate_trials) {
  latent_pop_tumor_gp_eta[s] ~ std_normal();
}