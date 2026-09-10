// tr/flags.stan
// Module: Total Rate (tr)
// Defines gating flags for total rate module. Population intercept always on.
// Multi-level hierarchy: flags are arrays indexed by level (1..n_levels)

int<lower=0,upper=1> enable_pop_cov_tr;  // population covariate linear model

// Per-level flags for hierarchical intercepts and slopes
// Index 1 = first grouping level (e.g., trial)
// Index n_levels = patient level
array[n_levels] int<lower=0,upper=4> enable_level_intercept_tr;
array[n_levels] int<lower=0,upper=1> enable_level_cov_tr;

// Process noise flags (patient-specific features, unchanged)
int<lower=0,upper=1> enable_patient_process_noise_tr;
int<lower=0,upper=1> enable_patient_process_noise_sd_tr;
int<lower=0,upper=1> enable_patient_process_noise_phi_tr;
int<lower=0,upper=1> enable_pop_process_noise_tr;

// ===== SD sub-hierarchy (issue #110) =====
// Mode matrix: [location-level L, sub-level ℓ]
// Entries with ℓ >= L must be 0 (validated in transformed_data via split_sd_cp_ncp_pos).
// 0=NONE, 1=FE, 2=RE, 4=RE_CP. (RE_GP=3 is rejected.)
array[n_levels, n_levels] int<lower=0, upper=4> enable_sd_level_intercept_mode_tr;
