// Include fragment: transformed parameters block content for SSLS models
// This file contains state-space model transformed parameters (states computation)
// Previously in legacy/sf-ssls-transformed_parameters.stan

// ============================================================================
// Combined patient-level rates from tumor regression and growth fraction modules
// ============================================================================

vector[n_patients] patient_log_decrease_rate = tr_loc_patient + frac_log_decrease_patient;
vector[n_patients] patient_log_growth_rate   = tr_loc_patient + frac_log_growth_patient;

// ============================================================================
// Growth lag parameters (to be refactored into module)
// ============================================================================

vector[n_patients] patient_log_growth_lag_effect = zeros_vector(n_patients);
vector[n_patients] patient_log_growth_lag = rep_vector(pop_log_growth_lag, n_patients);

if (!pop_growth_lag_param_only) {
  patient_log_growth_lag_effect = patient_log_growth_lag_sd * raw_patient_log_growth_lag;
  patient_log_growth_lag += patient_log_growth_lag_effect;
}

// ============================================================================
// Tumor GP correlation parameters (to be refactored into module)
// ============================================================================

vector[n_patients] patient_tumor_gp_rho = independ_long_process_noise ? zeros_vector(n_patients) : rep_vector(exp(log_pop_tumor_gp_rho), n_patients);  

// ============================================================================
// State matrices
// ============================================================================

// Backward compatibility: states at actual visit times for observation model
matrix[n_total_visits, 2] states; 

vector[n_patients] patient_decrease_rate = exp(patient_log_decrease_rate);
vector[n_patients] patient_growth_rate = exp(patient_log_growth_rate);

// Full grid: [n_patients × max_t_width]
// Position 1 = each patient's first visit (different absolute weeks)
// Position t = t weeks after first visit for that patient
array[2] matrix[n_patients, max_t_width] states_full_grid;

profile("states") {
  vector[independ_long_process_noise || pop_rho_param_only ? 0 : n_patients] log_patient_tumor_gp_rho_effect; 
  vector[independ_long_process_noise || pop_rho_param_only ? 0 : n_trials] log_trial_tumor_gp_rho_effect;
  if (!independ_long_process_noise && !pop_rho_param_only) {
    log_patient_tumor_gp_rho_effect = log_patient_tumor_gp_rho_sd * raw_log_patient_tumor_gp_rho_effect;
    log_trial_tumor_gp_rho_effect = zeros_vector(n_trials);
    patient_tumor_gp_rho = exp(log_pop_tumor_gp_rho + log_trial_tumor_gp_rho_effect[patient_trial] + log_patient_tumor_gp_rho_effect);
  }

  // ============================================================================
  // Compute states on FULL time grid for all patients (VECTORIZED)
  // Each patient gets states at every week from their min to max visit
  // Grid is RELATIVE to each patient's first visit, size = max_t_width
  // This enables time-varying covariates in the hazard model
  // ============================================================================
  
  profile("compute full states") {
    // Matrix multiplication: states = init + rate * time
    // Each row is one patient, each column is one time point
    // This computes all states for all patients in two matrix operations
    // time_since_first_visit defined in base_transformed_data.stan
    // Note: (time_since_first_visit - 1) maps index 1 to 0 weeks elapsed (baseline)
    states_full_grid[1] = init_log_decrease_patient * ones_row_vector(max_t_width) + (-patient_decrease_rate) * (time_since_first_visit - 1);
    states_full_grid[2] = init_log_growth_patient * ones_row_vector(max_t_width) + patient_growth_rate * (time_since_first_visit - 1);
  }
  
  // ============================================================================
  // Extract visit-time states for observation model (backward compatibility)
  // ============================================================================
  
  profile("extract states") {
    for (i in 1:n_patients) {
      int visit_start, visit_end;
      (visit_start, visit_end) = get_pos(patient_visit_pos, i);
      
      // t_patient_visit_idx gives indices relative to patient's first visit
      // These map directly to columns in states_full_grid[*][i, :]
      array[get_pos_size(patient_visit_pos, i)] int visit_indices = t_patient_visit_idx[visit_start:visit_end];
      
      // Extract states at actual visit times from full grid (vectorized row extraction)
      states[visit_start:visit_end, 1] = to_vector(states_full_grid[1][i, visit_indices]);
      states[visit_start:visit_end, 2] = to_vector(states_full_grid[2][i, visit_indices]);
    }
  }
}
