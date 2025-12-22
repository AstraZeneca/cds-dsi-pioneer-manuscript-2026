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

// Full grid: [n_patients × max_t_width] - only needed when process noise is ON
// Position 1 = each patient's first visit (different absolute weeks)
// Position t = t weeks after first visit for that patient
// When process noise is OFF, this is sized to 0 columns to save computation
array[2] matrix[n_patients, enable_patient_process_noise_tr ? max_t_width : 0] states_full_grid;

profile("states") {
  if (enable_patient_process_noise_tr) {
    // ============================================================================
    // PROCESS NOISE ON: Compute states on FULL time grid for all patients
    // Each patient gets states at every week from their min to max visit
    // Grid is RELATIVE to each patient's first visit, size = max_t_width
    // This enables time-varying covariates in the hazard model
    // ============================================================================

    profile("compute full states") {
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
    }

    // Extract visit-time states from full grid
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
  } else {
    // ============================================================================
    // PROCESS NOISE OFF: Use efficient matrix multiplication approach (original)
    // Only compute states at unique visit times, not full dense grid
    // This is the pre-process-noise implementation for optimal performance
    // ============================================================================

    profile("compute sparse states") {
      // Build increment matrices [n_patients × max_unique_visit]
      // Using unit time steps (delta_t = 1 week)
      matrix[n_patients, max_unique_visit] D_increments = rep_matrix(-patient_decrease_rate[, 1], max_unique_visit);
      matrix[n_patients, max_unique_visit] G_increments = rep_matrix(patient_growth_rate[, 1], max_unique_visit);

      // Single batched matrix multiply for ALL patients
      // visit_cumsum_mat is [n_pop_unique_visits × max_unique_visit] precomputed in transformed_data
      // Result is [n_pop_unique_visits × n_patients] - states at unique visit times only
      matrix[n_pop_unique_visits, n_patients] cumsum_d = visit_cumsum_mat * D_increments' + rep_matrix(init_log_decrease_patient', n_pop_unique_visits);
      matrix[n_pop_unique_visits, n_patients] cumsum_g = visit_cumsum_mat * G_increments' + rep_matrix(init_log_growth_patient', n_pop_unique_visits);

      // Extract states at actual patient visit times using precomputed index mapping
      for (i in 1:n_patients) {
        int visit_start, visit_end;
        (visit_start, visit_end) = get_pos(patient_visit_pos, i);

        // patient2pop_unique_visit_idx maps each patient's visits to population unique visit indices
        states[visit_start:visit_end, 1] = cumsum_d[patient2pop_unique_visit_idx[visit_start:visit_end], i];
        states[visit_start:visit_end, 2] = cumsum_g[patient2pop_unique_visit_idx[visit_start:visit_end], i];
      }
    }
  }
}
