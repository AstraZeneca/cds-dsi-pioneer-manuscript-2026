// Include fragment: transformed parameters block content for SSLS models
// This file contains state-space model transformed parameters (states computation)
// Previously in legacy/sf-ssls-transformed_parameters.stan

// ============================================================================
// Combined patient-level rates from tumor regression and growth fraction modules
// ============================================================================

// Time-varying rates if either pop or patient process noise is enabled
// Note: enable_any_process_noise_tr is computed in transformed data
matrix[n_patients, enable_any_process_noise_tr ? max_t_width : 1] patient_log_decrease_rate;
matrix[n_patients, enable_any_process_noise_tr ? max_t_width : 1] patient_log_growth_rate;

if (enable_pop_process_noise_tr && enable_patient_process_noise_tr) {
  // Both levels: pop + patient (additive)
  // log_rate[i,t] = log_baseline_rate[i] + pop_deviation[t] + patient_deviation[i,t]
  patient_log_decrease_rate = (tr_loc_patient + frac_log_decrease_patient) * ones_row_vector(max_t_width)
    + rep_matrix(tr_pop_process_noise, n_patients)
    + tr_patient_process_noise;
  patient_log_growth_rate = (tr_loc_patient + frac_log_growth_patient) * ones_row_vector(max_t_width)
    + rep_matrix(tr_pop_process_noise, n_patients)
    + tr_patient_process_noise;
} else if (enable_pop_process_noise_tr) {
  // Pop-level only: shared temporal trend
  // log_rate[i,t] = log_baseline_rate[i] + pop_deviation[t]
  patient_log_decrease_rate = (tr_loc_patient + frac_log_decrease_patient) * ones_row_vector(max_t_width)
    + rep_matrix(tr_pop_process_noise, n_patients);
  patient_log_growth_rate = (tr_loc_patient + frac_log_growth_patient) * ones_row_vector(max_t_width)
    + rep_matrix(tr_pop_process_noise, n_patients);
} else if (enable_patient_process_noise_tr) {
  // Patient-level only (existing behavior)
  // log_rate[i,t] = log_baseline_rate[i] + patient_deviation[i,t]
  patient_log_decrease_rate = (tr_loc_patient + frac_log_decrease_patient) * ones_row_vector(max_t_width) + tr_patient_process_noise;
  patient_log_growth_rate = (tr_loc_patient + frac_log_growth_patient) * ones_row_vector(max_t_width) + tr_patient_process_noise;
} else {
  // No process noise: constant rates (stored in single column)
  patient_log_decrease_rate[, 1] = tr_loc_patient + frac_log_decrease_patient;
  patient_log_growth_rate[, 1] = tr_loc_patient + frac_log_growth_patient;
}

// ============================================================================
// State matrices
// ============================================================================

// Backward compatibility: states at actual visit times for observation model
matrix[n_total_visits, 2] states; 

// Conditionally sized matrices for exponentiated rates
matrix[n_patients, enable_any_process_noise_tr ? max_t_width : 1] patient_decrease_rate;
matrix[n_patients, enable_any_process_noise_tr ? max_t_width : 1] patient_growth_rate;

patient_decrease_rate = exp(patient_log_decrease_rate);
patient_growth_rate = exp(patient_log_growth_rate);

// Full grid: [n_patients × max_t_width]
// Position 1 = each patient's first visit (different absolute weeks)
// Position t = t weeks after first visit for that patient
// Needed when: process noise is ON, or multistate uses time-varying covariates
array[2] matrix[n_patients, (enable_any_process_noise_tr || enable_ms_pop_time_varying_cov) ? max_t_width : 0] states_full_grid;

profile("states") {
  if (enable_patient_process_noise_tr) {
    // ============================================================================
    // PATIENT-LEVEL PROCESS NOISE: Per-patient loop required (each patient has unique temporal pattern)
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
  } else if (enable_pop_process_noise_tr) {
    // ============================================================================
    // POP-ONLY PROCESS NOISE: Vectorized computation exploiting shared temporal structure
    // ============================================================================
    // Key insight: rate[i,t] = baseline_rate[i] * exp(pop_noise[t])
    // Therefore: cumsum(rate[i,1:t-1]) = baseline_rate[i] * cumsum(exp(pop_noise[1:t-1]))
    // The cumsum of exp(pop_noise) is identical for all patients, so compute once and use outer product

    profile("compute full states") {
      // Baseline rates without process noise (patient-specific, time-invariant)
      vector[n_patients] baseline_decrease_rate = exp(tr_loc_patient + frac_log_decrease_patient);
      vector[n_patients] baseline_growth_rate = exp(tr_loc_patient + frac_log_growth_patient);

      // Shared temporal cumulative sum: same for all patients
      // pop_cumsum_exp[t] = sum_{k=1}^{t} exp(pop_noise[k])
      row_vector[max_t_width] pop_exp = exp(tr_pop_process_noise);
      row_vector[max_t_width] pop_cumsum_exp = cumulative_sum(pop_exp);

      // First timepoint: just the initial state
      states_full_grid[1][, 1] = init_log_decrease_patient;
      states_full_grid[2][, 1] = init_log_growth_patient;

      if (max_t_width > 1) {
        // Vectorized outer product: states[i,t] = init[i] +/- baseline_rate[i] * pop_cumsum_exp[t-1]
        // For decrease: state decreases, so subtract
        // For growth: state increases, so add
        states_full_grid[1][, 2:] = rep_matrix(init_log_decrease_patient, max_t_width - 1)
          - baseline_decrease_rate * pop_cumsum_exp[:(max_t_width - 1)];
        states_full_grid[2][, 2:] = rep_matrix(init_log_growth_patient, max_t_width - 1)
          + baseline_growth_rate * pop_cumsum_exp[:(max_t_width - 1)];
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
  } else if (enable_ms_pop_time_varying_cov) {
    // ============================================================================
    // PROCESS NOISE OFF, but multistate needs full grid for time-varying covariates: Use vectorized approach
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
