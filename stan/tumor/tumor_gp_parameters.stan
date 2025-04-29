// Patient-level GP
real<lower = 0> pop_tumor_gp_alpha;
real log_pop_tumor_gp_rho;

// Hierarchical length-scale parameter
real<lower = 0> log_trial_tumor_gp_rho_sd;
vector[add_trial_level_tumor_gp_param ? n_trials : 0] raw_log_trial_tumor_gp_rho_effect; 

real<lower = 0> log_patient_tumor_gp_rho_sd;
vector[n_patients] raw_log_patient_tumor_gp_rho_effect; 

vector[sum(n_patient_unique_visits)] patient_tumor_gp_eta;
