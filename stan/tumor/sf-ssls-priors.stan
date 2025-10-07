pop_log_total_rate ~ normal(pop_log_total_rate_mean, pop_log_total_rate_sd);
pop_decrease_frac_logit ~ normal(pop_decrease_frac_logit_mean, pop_decrease_frac_logit_sd);

trial_log_total_rate_sd ~ normal(0, trial_log_total_rate_sd_sd);
patient_log_total_rate_sd ~ normal(0, patient_log_total_rate_sd_sd);

raw_trial_log_total_rate ~ std_normal();
raw_patient_log_total_rate ~ std_normal();

// Priors on ORIGINAL-scale (beta) coefficients recovered from QR space.
// No Jacobian needed: sampling occurs in QR space (theta parameters named *_qr) but priors are placed on
// recovered original-scale coefficients (data * backsolve). R is data, so log|det R^{-1}| is constant and omitted.
pop_decrease_frac_logit_coef ~ normal(pop_decrease_frac_logit_coef_mean, pop_decrease_frac_logit_coef_sd);
to_vector(raw_trial_decrease_frac_logit_coef) ~ std_normal();
trial_decrease_frac_logit_coef_sd ~ normal(0, trial_decrease_frac_logit_coef_sd_sd);

if (add_patient_level_frac) {
  to_vector(raw_patient_decrease_frac_logit_coef) ~ std_normal();
  patient_decrease_frac_logit_coef_sd ~ normal(0, patient_decrease_frac_logit_coef_sd_sd);
  // Patient-level fraction intercept (hierarchical)
  patient_decrease_frac_logit_sd ~ normal(0, patient_decrease_frac_logit_sd_sd);
  raw_patient_decrease_frac_logit ~ std_normal();
}

// Growth lag and transition rate priors (missing previously)
pop_log_growth_lag ~ normal(growth_lag_mean, growth_lag_sd);
pop_log_growth_transition_rate ~ normal(0, log_growth_transition_rate_sd); // Center at 0 (log scale); adjust hyperparam name if mean provided later
patient_log_growth_lag_sd ~ normal(0, patient_log_growth_lag_sd_sd);
if (!pop_growth_lag_param_only) {
  raw_patient_log_growth_lag ~ std_normal();
}

// GP rho (length-scale) prior (population) and hierarchical patient variation
log_pop_tumor_gp_rho ~ lognormal(pop_tumor_gp_rho_meanlog, pop_tumor_gp_rho_sdlog);
log_patient_tumor_gp_rho_sd ~ normal(0, log_patient_tumor_gp_rho_sd_sd);
if (!independ_long_process_noise && !pop_rho_param_only) {
  raw_log_patient_tumor_gp_rho_effect ~ std_normal();
}

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

pop_decrease_prop_logis ~ normal(pop_decrease_prop_logis_mean, pop_decrease_prop_logis_sd);
trial_decrease_prop_logis_sd ~ normal(0, trial_decrease_prop_logis_sd_sd);
raw_trial_decrease_prop_logis ~ std_normal();
patient_decrease_prop_logis_sd ~ normal(0, patient_decrease_prop_logis_sd_sd);
// patient_decrease_prop_logis ~ normal(pop_decrease_prop_logis, patient_decrease_prop_logis_sd);
raw_patient_decrease_prop_logis ~ std_normal();

// log_lod ~ normal(log(lod), log_lod_sd);

// Priors for covariate effects using hyperparameters from data block (original space)
pop_decrease_prop_logis_coef ~ normal(pop_decrease_prop_logis_coef_mean, pop_decrease_prop_logis_coef_sd);


