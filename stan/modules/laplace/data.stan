// laplace/data.stan
//
// All patient data lives in unified arrays:
//   _base_hierarchy_data.stan  — patient_level_groups, covar_design_matrix
//   _base_data.stan            — n_patient_visits, t_patient_visits, ...
//   modules/psa/data.stan      — psa_values, psa_measured, ...
//   modules/multistate/data.stan — ms_final_state, ms_time_01, ...
//
// Routing: laplace_split_level + laplace_target_group (in _base_hierarchy_data.stan)
// determine which patients are HMC vs Laplace-marginalized.
// Index arrays hmc_patient_idx and laplace_patient_idx are computed in
// laplace/transformed_data.stan.
//
// No separate laplace_* data arrays are needed.
