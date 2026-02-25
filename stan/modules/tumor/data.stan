// ============================================================================
// TUMOR MEASUREMENT DATA
// ============================================================================
// Tumor/SLD-specific longitudinal measurements and outcomes.
// This module contains data aligned with the visit schedule defined in base_data.stan.
//
// Array lengths: All measurement arrays have length sum(n_patient_visits) and are
// aligned with t_patient_visits from base_data.stan.

// Sum of target lesion diameters (SLD) at each visit
vector<lower = 0>[sum(n_patient_visits)] sum_tumor_size; // cm

// RECIST response category at each visit
// 1 = CR (Complete Response)
// 2 = PR (Partial Response)
// 3 = SD (Stable Disease)
// 4 = PD (Progressive Disease)
// 5 = Not Evaluable
array[sum(n_patient_visits)] int<lower = 1, upper = 5> recist;

// ============================================================================
// TUMOR-BASED PROGRESSION OUTCOMES
// ============================================================================

// Progression-free survival based on tumor measurements
array[n_patients] int<lower = 0> pfs; // How many weeks after baseline did patient survive without progression
array[n_patients] int<lower = 0, upper = 1> right_censored;
array[n_patients] int<lower = 0> interval_censored; // Weeks after `pfs` that actual progression could have occurred

// Progression-free survival based on target lesions only (excluding non-target and new lesions)
array[n_patients] int<lower = 0> target_pfs;
array[n_patients] int<lower = 0, upper = 1> target_right_censored;

// Death timing
array[n_patients] int<lower = 0> death_week;

// ============================================================================
// MODEL CONFIGURATION
// ============================================================================

// Whether to fit the tumor dynamics model (0 = skip, 1 = fit)
int<lower = 0, upper = 1> fit_tumor_data;

// Number of time points for state-space trajectory replication
int<lower = 0> sf_rep_T;

// Debug mode flag
int<lower = 0, upper = 1> debug;

// Number of shards for parallel processing (map_rect)
int<lower = 1, upper = n_patients> n_shards;
