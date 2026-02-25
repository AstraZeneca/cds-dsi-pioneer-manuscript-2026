// tr/flags.stan
// Module: Total Rate (tr)
// Defines gating flags for total rate module. Population intercept always on.
// Multi-level hierarchy: flags are arrays indexed by level (1..n_levels)

int<lower=0,upper=1> enable_pop_cov_tr;  // population covariate linear model

// Per-level flags for hierarchical intercepts and slopes
// Index 1 = first grouping level (e.g., trial)
// Index n_levels = patient level
array[n_levels] int<lower=0,upper=1> enable_level_intercept_tr;
array[n_levels] int<lower=0,upper=1> enable_level_cov_tr;

// Process noise flags (patient-specific features, unchanged)
int<lower=0,upper=1> enable_patient_process_noise_tr;
int<lower=0,upper=1> enable_patient_process_noise_sd_tr;
int<lower=0,upper=1> enable_patient_process_noise_phi_tr;
int<lower=0,upper=1> enable_pop_process_noise_tr; 
