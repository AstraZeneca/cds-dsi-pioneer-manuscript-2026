pop_log_net_rate ~ normal(pop_log_net_rate_mean, pop_log_net_rate_sd);
pop_log_rate_ratio ~ normal(pop_log_rate_ratio_mean, pop_log_rate_ratio_sd);
patient_log_net_rate_sd ~ normal(0, patient_log_net_rate_sd_sd);
// patient_log_net_rate ~ normal(pop_log_net_rate, patient_log_net_rate_sd);
raw_patient_log_net_rate ~ std_normal(); 
// patient_log_rate_ratio_sd ~ normal(0, patient_log_rate_ratio_sd_sd);
// raw_patient_log_rate_ratio ~ std_normal();
trial_log_net_rate_sd ~ normal(0, trial_log_net_rate_sd_sd);
raw_trial_log_net_rate ~ std_normal(); 

pop_log_growth_lag ~ normal(growth_lag_mean, growth_lag_sd);
pop_log_growth_transition_rate ~ normal(0, log_growth_transition_rate_sd);
patient_log_growth_lag_sd ~ normal(0, patient_log_growth_lag_sd_sd);
// patient_log_growth_lag ~ normal(pop_log_growth_lag, patient_log_growth_lag_sd);
raw_patient_log_growth_lag ~ std_normal(); 

// pop_tumor_gp_alpha ~ normal(0, pop_tumor_gp_alpha_sd);
log_pop_tumor_gp_rho ~ normal(pop_tumor_gp_rho_meanlog, pop_tumor_gp_rho_sdlog);
log_patient_tumor_gp_rho_sd ~ normal(0, log_patient_tumor_gp_rho_sd_sd);
to_vector(raw_log_patient_tumor_gp_rho_effect) ~ std_normal(); 

to_vector(raw_patient_process_noise) ~ std_normal();

pop_process_sd[1] ~ normal(0, pop_decrease_process_sd_sd);
pop_process_sd[2] ~ normal(0, pop_growth_process_sd_sd);

measure_sd ~ normal(0, measure_sd_sd);

if (!independ_cross_process_noise) {
  L_process_corr ~ lkj_corr_cholesky(process_corr_param);
}

pop_decrease_prop_logis ~ normal(pop_decrease_prop_logis_mean, pop_decrease_prop_logis_sd);
trial_decrease_prop_logis_sd ~ normal(0, trial_decrease_prop_logis_sd_sd);
raw_trial_decrease_prop_logis ~ std_normal();
patient_decrease_prop_logis_sd ~ normal(0, patient_decrease_prop_logis_sd_sd);
// patient_decrease_prop_logis ~ normal(pop_decrease_prop_logis, patient_decrease_prop_logis_sd);
raw_patient_decrease_prop_logis ~ std_normal();

// log_lod ~ normal(log(lod), log_lod_sd);

// Priors for covariate effects using hyperparameters from data block
pop_log_net_rate_coef ~ normal(pop_log_net_rate_coef_mean, pop_log_net_rate_coef_sd);
pop_decrease_prop_logis_coef ~ normal(pop_decrease_prop_logis_coef_mean, pop_decrease_prop_logis_coef_sd);

// Jacobian adjustment for QR reparam
target += log_abs_det_R_covar_design_matrix;

if (!pop_rates_param_only && add_trial_level_net_rate) {
  trial_log_net_rate_coef_sd ~ normal(0, trial_log_net_rate_coef_sd_sd);
  to_vector(raw_trial_log_net_rate_coef) ~ std_normal();
}