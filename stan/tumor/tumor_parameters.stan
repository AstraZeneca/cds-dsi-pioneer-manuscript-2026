real<lower = 0> tumor_mean; // Actual tumor mean, not the lognormal mean
real<lower = 0> tumor_sd; // Actual tumor SD, not the lognormal one

real<lower = 0> pop_tumor_gp_rho;

// Multilevel intercepts

real<lower = 0> patient_tumor_gp_intercept_sd; 
vector[use_tumor_model && multilevel_patient ? n_patients : 0] patient_tumor_gp_intercept_effect;

vector<lower = 0>[use_tumor_model && multilevel_tumor ? n_patients : 0] tumor_gp_intercept_sd; 
vector[use_tumor_model && multilevel_tumor ? sum(n_patient_tumors) : 0] tumor_gp_intercept_effect;
