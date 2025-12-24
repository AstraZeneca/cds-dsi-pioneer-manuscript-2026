// Include fragment: transformed parameters block content for SSLS models
// This file contains state-space model transformed parameters (states computation)
// Previously in legacy/sf-ssls-transformed_parameters.stan

// ============================================================================
// Combined patient-level rates from tumor regression and growth fraction modules
// ============================================================================

matrix[n_patients, enable_patient_process_noise_tr ? max_t_width : 1] patient_log_decrease_rate; 
matrix[n_patients, enable_patient_process_noise_tr ? max_t_width : 1] patient_log_growth_rate;

if (enable_patient_process_noise_tr) {
  // Add AR(1) process noise (mean-reverting to baseline rate)
  // log_rate[i,t] = log_baseline_rate[i] + deviation[i,t]
  patient_log_decrease_rate = (tr_loc_patient + frac_log_decrease_patient) * ones_row_vector(max_t_width) + tr_patient_process_noise;
  patient_log_growth_rate   = (tr_loc_patient + frac_log_growth_patient) * ones_row_vector(max_t_width) + tr_patient_process_noise;
} else {
  // No process noise: constant rates (stored in single column)
  patient_log_decrease_rate[, 1] = tr_loc_patient + frac_log_decrease_patient;
  patient_log_growth_rate[, 1]   = tr_loc_patient + frac_log_growth_patient;
}

// ============================================================================
// State matrices
// ============================================================================

// Backward compatibility: states at actual visit times for observation model
matrix[n_total_visits, 2] states; 

// Conditionally sized matrices for exponentiated rates
matrix[n_patients, enable_patient_process_noise_tr ? max_t_width : 1] patient_decrease_rate;
matrix[n_patients, enable_patient_process_noise_tr ? max_t_width : 1] patient_growth_rate;

patient_decrease_rate = exp(patient_log_decrease_rate);
patient_growth_rate = exp(patient_log_growth_rate);

// Full grid: [n_patients × max_t_width]
// Position 1 = each patient's first visit (different absolute weeks)
// Position t = t weeks after first visit for that patient
// Needed when: process noise is ON, or other_events uses time-varying tumor covariates
array[2] matrix[n_patients, (enable_patient_process_noise_tr || oe_enable_pop_tumor_cov) ? max_t_width : 0] states_full_grid;

profile("states") {
  if (enable_patient_process_noise_tr) {
    // ============================================================================
    // PROCESS NOISE ON: Compute states using cumulative sum (time-varying rates)
    // ============================================================================

    profile("compute full states") {
      for (i in 1:n_patients) {
        states_full_grid[1][i, 1] = init_log_decrease_patient[i];
        states_full_grid[2][i, 1] = init_log_growth_patient[i];

        if (max_t_width > 1) {
          states_full_grid[1][i, 2:] = to_row_vector(init_log_decrease_patient[i] + cumulative_sum(-patient_decrease_rate[i, :(max_t_width-1)]));
          states_full_grid[2][i, 2:] = to_row_vector(init_log_growth_patient[i] + cumulative_sum(patient_growth_rate[i, :(max_t_width-1)]));
        }
      }
    }

    profile("extract states") {
      for (i in 1:n_patients) {
        int visit_start, visit_end;
        (visit_start, visit_end) = get_pos(patient_visit_pos, i);
        array[get_pos_size(patient_visit_pos, i)] int visit_indices = t_patient_visit_idx[visit_start:visit_end];
        states[visit_start:visit_end, 1] = to_vector(states_full_grid[1][i, visit_indices]);
        states[visit_start:visit_end, 2] = to_vector(states_full_grid[2][i, visit_indices]);
      }
    }
  } else if (oe_enable_pop_tumor_cov) {
    // ============================================================================
    // PROCESS NOISE OFF, but other_events needs full grid: Use vectorized approach
    // ============================================================================

    profile("compute full states") {
      // Vectorized outer product: states = init + rate * time
      states_full_grid[1] = init_log_decrease_patient * ones_row_vector(max_t_width) + (-patient_decrease_rate[, 1]) * (time_since_first_visit - 1);
      states_full_grid[2] = init_log_growth_patient * ones_row_vector(max_t_width) + patient_growth_rate[, 1] * (time_since_first_visit - 1);
    }

    profile("extract states") {
      for (i in 1:n_patients) {
        int visit_start, visit_end;
        (visit_start, visit_end) = get_pos(patient_visit_pos, i);
        array[get_pos_size(patient_visit_pos, i)] int visit_indices = t_patient_visit_idx[visit_start:visit_end];
        states[visit_start:visit_end, 1] = to_vector(states_full_grid[1][i, visit_indices]);
        states[visit_start:visit_end, 2] = to_vector(states_full_grid[2][i, visit_indices]);
      }
    }
  } else {
    // ============================================================================
    // PROCESS NOISE OFF, no full grid needed: Use sparse matrix multiplication
    // ============================================================================

    profile("compute sparse states") {
      matrix[n_patients, max_unique_visit] D_increments = rep_matrix(-patient_decrease_rate[, 1], max_unique_visit);
      matrix[n_patients, max_unique_visit] G_increments = rep_matrix(patient_growth_rate[, 1], max_unique_visit);

      matrix[n_pop_unique_visits, n_patients] cumsum_d = visit_cumsum_mat * D_increments' + rep_matrix(init_log_decrease_patient', n_pop_unique_visits);
      matrix[n_pop_unique_visits, n_patients] cumsum_g = visit_cumsum_mat * G_increments' + rep_matrix(init_log_growth_patient', n_pop_unique_visits);

      for (i in 1:n_patients) {
        int visit_start, visit_end;
        (visit_start, visit_end) = get_pos(patient_visit_pos, i);
        states[visit_start:visit_end, 1] = cumsum_d[patient2pop_unique_visit_idx[visit_start:visit_end], i];
        states[visit_start:visit_end, 2] = cumsum_g[patient2pop_unique_visit_idx[visit_start:visit_end], i];
      }
    }
  }
}
