// MIGRATED total rate priors -> modules/tr/priors.stan

// Growth lag and transition rate priors (missing previously)
pop_log_growth_lag ~ normal(growth_lag_mean, growth_lag_sd);
pop_log_growth_transition_rate ~ normal(0, log_growth_transition_rate_sd); // Center at 0 (log scale); adjust hyperparam name if mean provided later
patient_log_growth_lag_sd ~ normal(0, patient_log_growth_lag_sd_sd);

raw_patient_log_growth_lag ~ std_normal();

// GP rho (length-scale) prior (population) and hierarchical patient variation
log_pop_tumor_gp_rho ~ normal(pop_tumor_gp_rho_meanlog, pop_tumor_gp_rho_sdlog);
log_patient_tumor_gp_rho_sd ~ normal(0, log_patient_tumor_gp_rho_sd_sd);

raw_log_patient_tumor_gp_rho_effect ~ std_normal();

// If patient-level:
// patient_log_total_rate_effect ~ normal(0, sigma_patient_total);

pop_process_sd[1] ~ normal(0, pop_decrease_process_sd_sd);
pop_process_sd[2] ~ normal(0, pop_growth_process_sd_sd);

measure_sd ~ normal(0, measure_sd_sd);

if (!independ_cross_process_noise) {
  L_process_corr ~ lkj_corr_cholesky(process_corr_param);
}

// Process noise innovations (unused earlier but declared); treat as standard normal if model uses non-centered noise draws later
to_vector(raw_patient_process_noise) ~ std_normal();


// log_lod ~ normal(log(lod), log_lod_sd);

// Priors for covariate effects using hyperparameters from data block (original space)
