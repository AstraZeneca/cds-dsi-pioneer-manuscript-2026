// ============================================================================
// Base Hierarchy Data
// Shared across all models (full joint model, standalone multistate, etc.)
// ============================================================================

// Trial and Patient Level Data
int<lower=1> n_trials;
int<lower=0> n_patients;

// Number of patients handled by full HMC (dimensions per-patient parameters).
// When Laplace is disabled, n_hmc_patients == n_patients.
int<lower=0> n_hmc_patients;

// Laplace routing key: identifies which patients receive full HMC treatment.
// laplace_split_level: hierarchy level to split on (must be non-patient level, e.g. 1 = trial).
// laplace_target_group: group ID at that level whose patients get full HMC.
// Both are 0 when enable_laplace_nontarget = 0.
int<lower=0> laplace_split_level;
int<lower=0> laplace_target_group;

// Patient to Trial Mapping (kept for backward compatibility)
array[n_patients] int<lower=1, upper=n_trials> patient_trial;

// Multi-level hierarchy configuration
// Supports arbitrary N-level hierarchies (e.g., trial → region → site → patient)
int<lower=1> n_levels;
array[n_levels] int<lower=1> n_groups_per_level;
array[n_patients, n_levels] int<lower=1> patient_level_groups;

// ============================================================================
// PATIENT-LEVEL COVARIATES
// ============================================================================

int<lower=0> n_covar;
matrix[n_patients, n_covar] covar_design_matrix;
