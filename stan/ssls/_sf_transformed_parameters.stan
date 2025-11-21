// Include fragment: transformed parameters block content for SSLS models
// This file contains state-space model transformed parameters (states computation)
// Previously in legacy/sf-ssls-transformed_parameters.stan

// ============================================================================
// Combined patient-level rates from tumor regression and growth fraction modules
// ============================================================================

matrix[n_patients, max_t_width] patient_log_decrease_rate = (tr_loc_patient + frac_log_decrease_patient) * ones_row_vector(max_t_width);
matrix[n_patients, max_t_width] patient_log_growth_rate   = (tr_loc_patient + frac_log_growth_patient) * ones_row_vector(max_t_width);

if (enable_patient_process_noise_tr) {
  // Add AR(1) process noise (mean-reverting to baseline rate)
  // log_rate[i,t] = log_baseline_rate[i] + deviation[i,t]
  patient_log_decrease_rate += tr_patient_process_noise;
  patient_log_growth_rate   += tr_patient_process_noise;
}

// ============================================================================
// State matrices
// ============================================================================

// Backward compatibility: states at actual visit times for observation model
matrix[n_total_visits, 2] states; 

matrix[n_patients, max_t_width] patient_decrease_rate = exp(patient_log_decrease_rate);
matrix[n_patients, max_t_width] patient_growth_rate = exp(patient_log_growth_rate);

// Full grid: [n_patients × max_t_width]
// Position 1 = each patient's first visit (different absolute weeks)
// Position t = t weeks after first visit for that patient
array[2] matrix[n_patients, max_t_width] states_full_grid;

profile("states") {
  // ============================================================================
  // Compute states on FULL time grid for all patients
  // Each patient gets states at every week from their min to max visit
  // Grid is RELATIVE to each patient's first visit, size = max_t_width
  // This enables time-varying covariates in the hazard model
  // ============================================================================
  
  profile("compute full states") {
    if (enable_patient_process_noise_tr) {
      // Time-varying rates: use cumulative sum for numerical integration
      // Computational cost: O(n_patients × max_t_width)
      for (i in 1:n_patients) {
        // Initial states (t=1): no rates applied yet
        states_full_grid[1][i, 1] = init_log_decrease_patient[i];
        states_full_grid[2][i, 1] = init_log_growth_patient[i];

        if (max_t_width > 1) {
          // Subsequent states (t>1): initial state + cumulative sum of rates
          // state[t] = init + sum(rates[1:(t-1)])
          states_full_grid[1][i, 2:] = to_row_vector(init_log_decrease_patient[i] + cumulative_sum(-patient_decrease_rate[i, :(max_t_width-1)]));
          states_full_grid[2][i, 2:] = to_row_vector(init_log_growth_patient[i] + cumulative_sum(patient_growth_rate[i, :(max_t_width-1)]));
        }
      }
    } else {
      // Constant rates: simple vectorized computation (optimal for no process noise)
      // Computational cost: O(n_patients) - just vector-scalar multiplications
      // states[i,t] = init[i] + rate[i] * (t - 1)
      // Note: time_since_first_visit starts at 1, so subtract 1 to get [0, 1, 2, ...]
      states_full_grid[1] = init_log_decrease_patient * ones_row_vector(max_t_width) + -patient_decrease_rate[, 1] * (time_since_first_visit - 1);
      states_full_grid[2] = init_log_growth_patient * ones_row_vector(max_t_width) + patient_growth_rate[, 1] * (time_since_first_visit - 1);
    }
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

      // if (is_nan(states_full_grid[1][i, 1]) || is_nan(states_full_grid[2][i, 1])) {
      //   print("NaN detected in states computation for patient ", i, ".");
      //   print(i, ": patient_decrease_rate[i] = ", patient_decrease_rate[i]);
      //   print(i, ": patient_growth_rate[i] = ", patient_growth_rate[i]);
      //   print(i, ": tr_phi_patient_process_noise[i] = ", tr_phi_patient_process_noise[i]);
      //   print(i, ": tr_logit_phi_pop_process_noise = ", tr_logit_phi_pop_process_noise);
      //   print(i, ": tr_raw_patient_phi_process_noise[i] * tr_sd_patient_phi_process_noise = ", 
      //     tr_raw_patient_phi_process_noise[i] * tr_sd_patient_phi_process_noise);
      // }
    }
  }
}
