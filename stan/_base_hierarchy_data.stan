// ============================================================================
// Base Hierarchy Data
// Shared across all models (full joint model, standalone multistate, etc.)
// ============================================================================

// Trial and Patient Level Data
int<lower=1> n_trials;
int<lower=0> n_patients;

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
