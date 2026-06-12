// Include fragment: transformed parameters block content for SSLS models
// This file contains state-space model transformed parameters (states computation)
// Previously in legacy/sf-ssls-transformed_parameters.stan

// ============================================================================
// Combined patient-level rates from tumor regression and growth fraction modules
// ============================================================================

// Time-varying rates if either pop or patient process noise is enabled
// Note: enable_any_process_noise_tr is computed in transformed data
matrix[n_forecast_patients, enable_any_process_noise_tr ? max_t_width : 1] patient_log_decrease_rate;
matrix[n_forecast_patients, enable_any_process_noise_tr ? max_t_width : 1] patient_log_growth_rate;

if (enable_pop_process_noise_tr && enable_patient_process_noise_tr) {
  // Both levels: pop + patient (additive)
  // log_rate[i,t] = log_baseline_rate[i] + pop_deviation[t] + patient_deviation[i,t]
  patient_log_decrease_rate = (tr_loc_patient + frac_log_decrease_patient) * ones_row_vector(max_t_width)
    + rep_matrix(tr_pop_process_noise, n_forecast_patients)
    + tr_patient_process_noise;
  patient_log_growth_rate = (tr_loc_patient + frac_log_growth_patient) * ones_row_vector(max_t_width)
    + rep_matrix(tr_pop_process_noise, n_forecast_patients)
    + tr_patient_process_noise;
} else if (enable_pop_process_noise_tr) {
  // Pop-level only: shared temporal trend
  // log_rate[i,t] = log_baseline_rate[i] + pop_deviation[t]
  patient_log_decrease_rate = (tr_loc_patient + frac_log_decrease_patient) * ones_row_vector(max_t_width)
    + rep_matrix(tr_pop_process_noise, n_forecast_patients);
  patient_log_growth_rate = (tr_loc_patient + frac_log_growth_patient) * ones_row_vector(max_t_width)
    + rep_matrix(tr_pop_process_noise, n_forecast_patients);
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
matrix[n_forecast_visits, 2] states;

// Conditionally sized matrices for exponentiated rates
matrix[n_forecast_patients, enable_any_process_noise_tr ? max_t_width : 1] patient_decrease_rate;
matrix[n_forecast_patients, enable_any_process_noise_tr ? max_t_width : 1] patient_growth_rate;

patient_decrease_rate = exp(patient_log_decrease_rate);
patient_growth_rate = exp(patient_log_growth_rate);

// Full grid: [n_patients × max_t_width]
// Position 1 = each patient's first visit (different absolute weeks)
// Position t = t weeks after first visit for that patient
// Needed when: process noise is ON, explicit flag is set, or MS time-varying covariates are active
array[2] matrix[n_forecast_patients, need_states_full_grid ? max_t_width : 0] states_full_grid;

profile("states") {
  if (enable_patient_process_noise_tr) {
    // ============================================================================
    // PATIENT-LEVEL PROCESS NOISE: Per-patient loop required (each patient has unique temporal pattern)
    // ============================================================================

    profile("compute full states") {
      for (j in 1:n_forecast_patients) {
        states_full_grid[1][j, 1] = init_log_decrease_patient[j];
        states_full_grid[2][j, 1] = init_log_growth_patient[j];

        if (max_t_width > 1) {
          states_full_grid[1][j, 2:] = to_row_vector(init_log_decrease_patient[j] + cumulative_sum(-patient_decrease_rate[j, :(max_t_width-1)]));
          if (enable_gr_decay) {
            // phi-differences over elapsed time so cumsum telescopes to growth_rate*phi(t).
            row_vector[max_t_width] warp = growth_warp(time_since_first_visit - 1, gr_decay_kappa[j]);
            row_vector[max_t_width-1] warp_diff = warp[2:] - warp[:(max_t_width-1)];
            states_full_grid[2][j, 2:] = to_row_vector(init_log_growth_patient[j]
              + cumulative_sum(patient_growth_rate[j, :(max_t_width-1)] .* warp_diff));
          } else {
            states_full_grid[2][j, 2:] = to_row_vector(init_log_growth_patient[j] + cumulative_sum(patient_growth_rate[j, :(max_t_width-1)]));
          }
        }
      }
    }

    profile("extract states") {
      for (j in 1:n_forecast_patients) {
        int p = forecast_patient_idx[j];
        int data_start, data_end;
        (data_start, data_end) = get_pos(patient_visit_pos, p);
        int state_start, state_end;
        (state_start, state_end) = get_pos(forecast_visit_pos, j);
        array[data_end - data_start + 1] int visit_indices = t_patient_visit_idx[data_start:data_end];
        states[state_start:state_end, 1] = to_vector(states_full_grid[1][j, visit_indices]);
        states[state_start:state_end, 2] = to_vector(states_full_grid[2][j, visit_indices]);
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
      vector[n_forecast_patients] baseline_decrease_rate = exp(tr_loc_patient + frac_log_decrease_patient);
      vector[n_forecast_patients] baseline_growth_rate = exp(tr_loc_patient + frac_log_growth_patient);

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
        if (enable_gr_decay) {
          // Outer-product factorization breaks: the decay weight depends on BOTH patient
          // (via kappa_i) and time. Per-patient phi-difference cumsum (mirrors B1).
          for (i in 1:n_forecast_patients) {
            row_vector[max_t_width] warp = growth_warp(time_since_first_visit - 1, gr_decay_kappa[i]);
            row_vector[max_t_width-1] warp_diff = warp[2:] - warp[:(max_t_width-1)];
            states_full_grid[2][i, 2:] = init_log_growth_patient[i]
              + cumulative_sum(baseline_growth_rate[i] * (pop_exp[:(max_t_width-1)] .* warp_diff));
          }
        } else {
          states_full_grid[2][, 2:] = rep_matrix(init_log_growth_patient, max_t_width - 1)
            + baseline_growth_rate * pop_cumsum_exp[:(max_t_width - 1)];
        }
      }
    }

    profile("extract states") {
      for (j in 1:n_forecast_patients) {
        int p = forecast_patient_idx[j];
        int data_start, data_end;
        (data_start, data_end) = get_pos(patient_visit_pos, p);
        int state_start, state_end;
        (state_start, state_end) = get_pos(forecast_visit_pos, j);
        array[data_end - data_start + 1] int visit_indices = t_patient_visit_idx[data_start:data_end];
        states[state_start:state_end, 1] = to_vector(states_full_grid[1][j, visit_indices]);
        states[state_start:state_end, 2] = to_vector(states_full_grid[2][j, visit_indices]);
      }
    }
  } else if (need_states_full_grid) {
    // ============================================================================
    // PROCESS NOISE OFF, but downstream modules need full grid: Use vectorized approach
    // (Triggered by enable_states_full_grid flag OR MS time-varying covariates)
    // ============================================================================

    profile("compute full states") {
      // Vectorized outer product: states = init + rate * time
      states_full_grid[1] = init_log_decrease_patient * ones_row_vector(max_t_width) + (-patient_decrease_rate[, 1]) * (time_since_first_visit - 1);
      if (enable_gr_decay) {
        // Per-patient warp of elapsed time (time_since_first_visit - 1), since kappa_i differs.
        for (i in 1:n_forecast_patients) {
          states_full_grid[2][i] = init_log_growth_patient[i]
            + patient_growth_rate[i, 1] * growth_warp(time_since_first_visit - 1, gr_decay_kappa[i]);
        }
      } else {
        states_full_grid[2] = init_log_growth_patient * ones_row_vector(max_t_width) + patient_growth_rate[, 1] * (time_since_first_visit - 1);
      }
    }

    profile("extract states") {
      for (j in 1:n_forecast_patients) {
        int p = forecast_patient_idx[j];
        int data_start, data_end;
        (data_start, data_end) = get_pos(patient_visit_pos, p);
        int state_start, state_end;
        (state_start, state_end) = get_pos(forecast_visit_pos, j);
        array[data_end - data_start + 1] int visit_indices = t_patient_visit_idx[data_start:data_end];
        states[state_start:state_end, 1] = to_vector(states_full_grid[1][j, visit_indices]);
        states[state_start:state_end, 2] = to_vector(states_full_grid[2][j, visit_indices]);
      }
    }
  } else {
    // ============================================================================
    // PROCESS NOISE OFF, no full grid needed: Direct visit-time computation
    // ============================================================================
    // With constant rates, state is a linear function of time from baseline:
    //   state_d[v] = init_d - decrease_rate * (visit_week - baseline_week)
    //   state_g[v] = init_g + growth_rate * (visit_week - baseline_week)

    profile("compute visit states") {
      for (j in 1:n_forecast_patients) {
        int p = forecast_patient_idx[j];
        int data_start, data_end;
        (data_start, data_end) = get_pos(patient_visit_pos, p);
        int state_start, state_end;
        (state_start, state_end) = get_pos(forecast_visit_pos, j);

        int baseline_visit_idx = data_start + n_patient_screening_visits[p] - 1;
        int patient_baseline_week = t_patient_visits[baseline_visit_idx];

        // Time from baseline to each visit (weeks)
        vector[data_end - data_start + 1] dt =
          to_vector(t_patient_visits[data_start:data_end]) - patient_baseline_week;

        states[state_start:state_end, 1] = init_log_decrease_patient[j] - patient_decrease_rate[j, 1] * dt;
        if (enable_gr_decay) {
          vector[data_end - data_start + 1] dt_warp;
          for (v in 1:rows(dt)) dt_warp[v] = growth_warp(dt[v], gr_decay_kappa[j]);
          states[state_start:state_end, 2] = init_log_growth_patient[j] + patient_growth_rate[j, 1] * dt_warp;
        } else {
          states[state_start:state_end, 2] = init_log_growth_patient[j] + patient_growth_rate[j, 1] * dt;
        }
      }
    }
  }
}
