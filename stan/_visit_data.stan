// ============================================================================
// Visit Schedule Data
// Shared by all models (full biomarker models and standalone multistate).
//
// Requires in scope:
//   n_patients — from _hierarchy_data.stan
// ============================================================================

array[n_patients] int<lower=1> n_patient_visits;
array[sum(n_patient_visits)] int t_patient_visits;
