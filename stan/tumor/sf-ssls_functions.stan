row_vector sf_log_space_transition(row_vector current_x, 
                                   real time_next, real time_current,
                                   real decrease_rate, real growth_rate) {
  // Calculate time interval
  real delta_t = time_next - time_current;
  
  if (delta_t <= 0) {
      reject("Time must move forward; delta_t = ", delta_t);
  }
  
  // Expected means after time interval
  row_vector[2] expected_x = current_x + [-decrease_rate, growth_rate] * delta_t;
  
  return expected_x;
}

tuple(row_vector, vector) sf_log_space_transition(row_vector current_x, 
                                          real time_next, real time_current,
                                          real decrease_rate, real growth_rate,
                                          vector process_sd) {
  // Calculate time interval
  real delta_t = time_next - time_current;
  
  if (delta_t <= 0) {
      reject("Time must move forward; delta_t = ", delta_t);
  }
  
  // Expected means after time interval
  row_vector[2] expected_x = sf_log_space_transition(current_x, time_next, time_current, decrease_rate, growth_rate); 
  
  // Scale process noise standard deviations with sqrt(delta_t)
  vector[2] log_scaled_sd = log(process_sd) + 0.5 * log(delta_t);
  
  return (expected_x, log_scaled_sd);
}
                                         

/**
 * Log probability density function for state transitions in the log-space Stein-Fojo model
 *
 * This function calculates the log probability density of observing the next state given 
 * the current state and parameters in the log-transformed Stein-Fojo tumor growth model.
 * It properly handles irregular time points by scaling process noise with sqrt(delta_t).
 *
 * @param next_x Array of next state values (log-space)
 * @param current_x Array of current state values (log-space)
 * @param time_next Next observation time
 * @param time_current Current observation time
 * @param decrease_rate Tumor regression rate constant (d)
 * @param growth_rate Tumor growth rate constant (g)
 * @param process_sd_reg Process noise SD for regression component
 * @param process_sd_growth Process noise SD for growth component
 *
 * @return Log probability density
 */
real sf_log_space_transition_lpdf(row_vector next_x, row_vector current_x, 
                                 real time_next, real time_current,
                                 real decrease_rate, real growth_rate,
                                 vector process_sd, matrix L_process_corr) {
    // Calculate time interval
    real delta_t = time_next - time_current;
    
    if (delta_t <= 0) {
        reject("Time must move forward; delta_t = ", delta_t);
    }
    
    // Expected means after time interval
    row_vector[2] expected_x; 
    vector[2] log_scaled_sd; 
    (expected_x, log_scaled_sd) = sf_log_space_transition(current_x, time_next, time_current, decrease_rate, growth_rate, process_sd);
    
    matrix[2, 2] scaled_L_process_cov = diag_pre_multiply(exp(log_scaled_sd), L_process_corr);
    
    // Log probability for each component
    return multi_normal_cholesky_lpdf(next_x | expected_x, scaled_L_process_cov);
}

tuple(row_vector, row_vector) sf_log_space_transition_ncp(row_vector raw_next_x, row_vector current_x, 
                                                          real time_next, real time_current,
                                                          real decrease_rate, real growth_rate,
                                                          vector process_sd, matrix L_process_corr) {
                                                          // real process_sd_reg, real process_sd_growth) {
    row_vector[2] expected_x;
    vector[2] log_scaled_sd;
    (expected_x, log_scaled_sd) = sf_log_space_transition(current_x, time_next, time_current, decrease_rate, growth_rate, process_sd);
    
    matrix[2, 2] scaled_L_process_cov = diag_pre_multiply(exp(log_scaled_sd), L_process_corr);
    
    return (expected_x, expected_x + raw_next_x * scaled_L_process_cov'); 
}

tuple(row_vector, row_vector) sf_log_space_transition_ncp(row_vector current_x, 
                                                          real time_next, real time_current,
                                                          real decrease_rate, real growth_rate,
                                                          row_vector process_noise) {
    row_vector[2] expected_x = sf_log_space_transition(current_x, time_next, time_current, decrease_rate, growth_rate);
    
    return (expected_x, expected_x + process_noise); 
}

/**
 * Log probability density function for observations in the log-space Stein-Fojo model
 *
 * This function calculates the log probability density of the observation given
 * the current state, and handles left-censoring for measurements below LOD.
 *
 * @param y Observation (normalized)
 * @param x State values (log-space)
 * @param measure_sd Measurement noise standard deviation
 * @param log_lod Log of limit of detection (for censoring)
 *
 * @return Log probability density
 */
real sf_log_space_obs_lpdf(vector normalized_y, matrix x, real measure_sd, real log_normalized_lod) {
  assert_equal(rows(normalized_y), rows(x));
  int T = rows(normalized_y);
  real lp = 0;
  
  for (t in 1:T) {
    // Expected log observation using log-sum-exp
    real log_pred = log_sum_exp(x[t]);
    
    // Handle censoring
    if (normalized_y[t] > 0) {
      lp += normal_lpdf(log(normalized_y[t]) | log_pred, measure_sd);
    } else {
      lp += normal_lcdf(log_normalized_lod | log_pred, measure_sd);
    }
  }
  
  return lp;
}

/**
 * Calculate log likelihood for a complete patient trajectory
 *
 * @param y Array of observations
 * @param x Array of state vectors
 * @param times Array of observation times
 * @param decrease_rate Tumor regression rate
 * @param growth_rate Tumor growth rate
 * @param process_sd_reg Process noise SD for regression
 * @param process_sd_growth Process noise SD for growth
 * @param measure_sd Measurement noise SD
 * @param log_lod Log of limit of detection
 * @param alpha Initial state parameter
 * @param state_sd Initial state uncertainty
 *
 * @return Log likelihood
 */
real sf_log_space_trajectory_lpdf(matrix x, row_vector x0, array[] int times, real decrease_rate, real growth_rate, vector process_sd, matrix L_process_corr) {
  int T = rows(x) + 1;
  real log_prob = 0;
  
  // State transitions
  for (t in 1:(T - 1)) {
    log_prob += sf_log_space_transition_lpdf(x[t] | t > 1 ? x[t - 1] : x0, times[t + 1], times[t],
                                            decrease_rate, growth_rate,
                                            process_sd, L_process_corr);
  }
  
  return log_prob;
}

vector get_growth_lag_factor(vector time_points, real growth_lag, real transition_rate) {
  return inv_logit((time_points - growth_lag) / transition_rate);
}

vector get_growth_lag_factor(array[] real time_points, real growth_lag, real transition_rate) {
  return inv_logit((to_vector(time_points) - growth_lag) / transition_rate);
}

tuple(matrix, matrix) sf_log_space_trajectory_ncp(matrix raw_x, row_vector x0, array[] real times, real decrease_rate, real growth_rate, vector process_sd, matrix L_process_corr) {
  return sf_log_space_trajectory_ncp(raw_x, x0, times, decrease_rate, growth_rate, negative_infinity(), 1, process_sd, L_process_corr); 
}

tuple(matrix, matrix) sf_log_space_trajectory_ncp(
  matrix raw_x, row_vector x0, array[] real times,
  real decrease_rate, real growth_rate, real growth_lag, real transition_rate,
  vector process_sd, matrix L_process_corr
) {
  int T = rows(raw_x) + 1;
  matrix[T, 2] x;
  matrix[T, 2] expected_x;
  x[1] = x0;
  expected_x[1] = x0;
  
  vector[T] time_varying_factor = get_growth_lag_factor(times, growth_lag, transition_rate); 
  
  // State transitions
  for (t in 2:T) {
    (expected_x[t], x[t]) = sf_log_space_transition_ncp(raw_x[t - 1], x[t - 1], times[t], times[t - 1],
                                                        decrease_rate, time_varying_factor[t] * growth_rate,
                                                        process_sd, L_process_corr);
  }
  
  return (expected_x, x);
}

tuple(matrix, matrix) sf_log_space_trajectory_ncp(
  row_vector x0, array[] real times,
  real decrease_rate, real growth_rate, real growth_lag, real transition_rate,
  matrix process_noise
) {
  return sf_log_space_trajectory_ncp(x0, times, decrease_rate, growth_rate, growth_lag, transition_rate, process_noise, 0); 
}

tuple(matrix, matrix) sf_log_space_trajectory_ncp(
  row_vector x0, array[] real times,
  real decrease_rate, real growth_rate, real growth_lag, real transition_rate,
  matrix process_noise, int debug
) {
  int T = size(times);
  matrix[T, 2] x;
  matrix[T, 2] expected_x;
  x[1] = x0;
  expected_x[1] = x0;
  
  vector[T] time_varying_factor = get_growth_lag_factor(times, growth_lag, transition_rate); 
  
  if (debug) {
    print("times = ", times, ", time_varying_factor = ", time_varying_factor, ", decrease_rate = ", decrease_rate, ", growth_rate = ", growth_rate);
    print("noise = ", process_noise);
    print("log x[1] = ", x0);
  }
  
  // State transitions
  for (t in 2:T) {
    (expected_x[t], x[t]) = sf_log_space_transition_ncp(x[t - 1], times[t], times[t - 1], decrease_rate, time_varying_factor[t] * growth_rate, process_noise[t - 1]);
    
    if (debug) {
      print("log x[", t, "] = (", expected_x[t], ", ", x[t], "), x[", t, "] = (", exp(expected_x[t]), ", ", exp(x[t]), ")");
    }
  }
  
  return (expected_x, x);
}

/**
 * Parallel computation of tumor growth states for multiple patients using map_rect
 * 
 * This function computes the evolution of tumor states over time for multiple patients,
 * optionally using parallel processing via map_rect with configurable sharding.
 * 
 * @param visit_pos Position array for patient visits (size n_patients+1)
 * @param t_visits Array of time points for all visits across all patients
 * @param rho Vector of GP length scale parameters for each patient
 * @param delta Small value for numerical stability in GP calculations
 * @param process_sd Process noise standard deviations [decrease, growth]
 * @param L_process_corr Cholesky factor of 2x2 process correlation matrix
 * @param raw_process_noise Raw process noise values (n_visits-1 x 2 per patient)
 * @param lognormal_noise Flag for lognormal vs normal process noise
 * @param independ_long_process_noise Flag for independent longitudinal process noise
 * @param independ_cross_process_noise Flag for independent cross-component process noise
 * @param initial_states Matrix of initial states for all patients (n_patients x 2)
 * @param decrease_rate Vector of tumor decrease rates for each patient
 * @param growth_rate Vector of tumor growth rates for each patient
 * @param growth_lag Vector of growth lag parameters for each patient
 * @param growth_transition_rate Growth transition rate parameter (shared)
 * @param n_shards Number of shards for parallel processing (1 = serial)
 * @param debug Flag to enable debug printing
 * @return Matrix of computed states for all visits (n_total_visits x 2)
 * 
 * Packing Structure for Parallel Execution:
 * =========================================
 * 
 * For n_shards > 1, data is packed into shards where each shard processes
 * approximately n_patients/n_shards patients. The packing structure is:
 * 
 * x_is[shard] (integer data):
 * ┌─────────────────────┬─────────────────────┬─────┬─────────────────────┐
 * │  Common Section     │   Patient 1 Data    │ ... │   Patient k Data    │
 * └─────────────────────┴─────────────────────┴─────┴─────────────────────┘
 * 
 * Common Section:
 * [0]: Number of patients in shard
 * [1:n]: x_is patient positions array
 * [n+1:m]: theta patient positions array  
 * [m+1]: debug flag
 * [m+2]: lognormal_noise flag
 * 
 * Per-Patient Data (repeated for each patient in shard):
 * - x_is positions metadata
 * - Process noise flags
 * - Visit time points
 * - Theta positions metadata
 * 
 * thetas[shard] (real parameters):
 * ┌─────────────────┬─────┬─────────────────┐
 * │  Patient 1      │ ... │  Patient k      │
 * └─────────────────┴─────┴─────────────────┘
 * 
 * Per-Patient Parameters:
 * - GP length scale (rho)
 * - Process SDs (2 values)
 * - Process correlation
 * - Raw process noise (flattened)
 * - Initial states (2 values)
 * - Rate parameters (decrease, growth, lag, transition)
 */
matrix calc_states(
  data array[] int visit_pos, data array[] int t_visits, vector rho, data real delta, vector process_sd, matrix L_process_corr, 
  matrix raw_process_noise, data int lognormal_noise, data int independ_long_process_noise, data int independ_cross_process_noise, 
  matrix initial_states, vector decrease_rate, vector growth_rate, 
  vector growth_lag, real growth_transition_rate, data int n_shards, data int debug
) {
  int n_patients = size(visit_pos) - 1;
  
  // Create position array for time points with one fewer elements per patient
  array[n_patients + 1] int visit_m1_pos = create_pos(get_pos_size(visit_pos), -1);  
  
  assert_equal(visit_pos[n_patients + 1], size(t_visits) + 1);
  
  matrix[size(t_visits), 2] states;
 
  if (n_shards == 1) { 
    for (i in 1:n_patients) {
      int visit_start, visit_end, n_visits;
      (visit_start, visit_end) = get_pos(visit_pos, i);
      n_visits = get_pos_size(visit_pos, i);
      
      states[visit_start:visit_end] = calc_patient_states(
        initial_states[i], get_int_sub_array(t_visits, visit_pos, i), decrease_rate[i], growth_rate[i], growth_lag[i], growth_transition_rate, 
        get_sub_vert_matrix(raw_process_noise, visit_m1_pos, i), lognormal_noise, rho[i], delta, process_sd, L_process_corr,
        independ_long_process_noise, independ_cross_process_noise
      ).2;
    }
  } else {
    array[n_shards] int x_is_shard_sizes = zeros_int_array(n_shards);
    array[n_shards] int theta_shard_sizes = zeros_int_array(n_shards);
    
    int theta_pos_size = 6;
    array[n_patients, theta_pos_size + 1] int theta_pos;
    
    int x_pos_size = 3;
    array[n_patients, x_pos_size + 1] int x_is_pos;
    
    int base_shard_size = n_patients %/% n_shards;
    int remainder_patients = n_patients % n_shards; 
    array[n_shards] int patients_in_shard = rep_array(base_shard_size, n_shards);
    array[n_patients] int patient_shard, patient_shard_idx;
    
    if (remainder_patients > 0) {
      patients_in_shard[:remainder_patients] = rep_array(base_shard_size + 1, remainder_patients);
    }
    
    array[14] int shard_common_pos = create_pos({
      1, // Shard ID
      1, // num of patients per shard
      1, // max number of patients in a shard  
      1, // Size of this pos
      base_shard_size + (remainder_patients > 0) + 3, // common + x_is patient pos
      base_shard_size + (remainder_patients > 0) + 1, // theta patient pos
      1,                    // total state size for this shard
      1,                    // debug flag
      1,                     // lognormal noise flag
      2,                    // independence flags
      1,                    // x_i_pos size
      1,                    // theta_pos size
      1                     // shard_common_pos_size
    });
    int shard_common_size = get_pos_total_size(shard_common_pos);
    
    array[n_patients] int shard_patient_idx = zeros_int_array(n_patients);
    
    array[n_shards, base_shard_size + (remainder_patients > 0)] int x_is_shard_patient_size = rep_array(0, n_shards, base_shard_size + (remainder_patients > 0)), 
                                                                    theta_shard_patient_size = rep_array(0, n_shards, base_shard_size + (remainder_patients > 0));
   
    array[n_shards] int total_state_size = zeros_int_array(n_shards);
    int curr_patient_shard = 1;
    
    // Calculate indices for each patient based on their number of visits
    for (i in 1:n_patients) {
      if (((i - 1) %/% sum(patients_in_shard[:curr_patient_shard])) > 0) {
        curr_patient_shard += 1;
      }
      
      patient_shard[i] = curr_patient_shard;
      shard_patient_idx[i] = i - (curr_patient_shard > 1 ? sum(patients_in_shard[:(curr_patient_shard - 1)]) : 0); 
      
      int n_visits = get_pos_size(visit_pos, i);
      
      total_state_size[curr_patient_shard] += n_visits;
     
      x_is_pos[i] = create_pos({
        x_pos_size + 1, // the x_is_pos 
        n_visits, // time points
        theta_pos_size + 1   // theta positions
      });
      
      theta_pos[i] = create_pos({ 
        1, // rho (length scale)
        2, // process_sd
        1, // L_process_corr[2, 1] 
        (n_visits - 1) * 2, // raw_process_noise for this patient
        2, // initial states 
        4 // rates (2), lag, transition
      }); 
      
      int curr_x_i_size = get_pos_total_size(x_is_pos[i]), curr_theta_size = get_pos_total_size(theta_pos[i]);
      x_is_shard_patient_size[patient_shard[i], shard_patient_idx[i]] = curr_x_i_size; 
      theta_shard_patient_size[patient_shard[i], shard_patient_idx[i]] = curr_theta_size; 
      
      x_is_shard_sizes[patient_shard[i]] += curr_x_i_size;
      theta_shard_sizes[patient_shard[i]] += curr_theta_size;
    }
    
    array[n_shards, 2 + base_shard_size + (remainder_patients > 0) + 1] int x_is_shard_patient_pos;
    array[n_shards, base_shard_size + (remainder_patients > 0) + 1] int theta_shard_patient_pos;
    
    int max_x_is_shard_size = max(x_is_shard_sizes), max_theta_shard_size = max(theta_shard_sizes);
    
    vector[2] phi = process_sd; // Shared parameters vector
    
    // Initialize parameter arrays for map_rect
    int x_is_size = shard_common_size + size(shard_common_pos) + max_x_is_shard_size;
    array[n_shards, x_is_size] int x_is = rep_array(-1111, n_shards, x_is_size);
    array[n_shards] vector[max_theta_shard_size] thetas = rep_array(rep_vector(negative_infinity(), max_theta_shard_size), n_shards);
    
    for (d in 1:n_shards) { // Shard common
      x_is_shard_patient_pos[d] = create_pos(append_array({ shard_common_size, size(shard_common_pos) }, x_is_shard_patient_size[d]));
      theta_shard_patient_pos[d] = create_pos(theta_shard_patient_size[d]);
      
      int x_is_pos_start, x_is_pos_end, theta_pos_start, theta_pos_end;
      (x_is_pos_start, x_is_pos_end) = get_pos(shard_common_pos, 5);
      (theta_pos_start, theta_pos_end) = get_pos(shard_common_pos, 6);
     
      x_is[d, 1] = d; 
      x_is[d, 2] = 0;
      x_is[d, 3] = base_shard_size + (remainder_patients > 0);
      x_is[d, 4] = shard_common_size;
      
      x_is[d, x_is_pos_start:x_is_pos_end] = x_is_shard_patient_pos[d]; 
      x_is[d, theta_pos_start:theta_pos_end] = theta_shard_patient_pos[d]; 
      x_is[d, theta_pos_end + 1] = total_state_size[d];
      x_is[d, theta_pos_end + 2] = debug;
      x_is[d, theta_pos_end + 3] = lognormal_noise;
      x_is[d, theta_pos_end + 4] = independ_long_process_noise;
      x_is[d, theta_pos_end + 5] = independ_cross_process_noise;
      x_is[d, theta_pos_end + 6] = x_pos_size;
      x_is[d, theta_pos_end + 7] = theta_pos_size;
      x_is[d, theta_pos_end + 8] = size(shard_common_pos) - 1;
     
      int shard_common_pos_start, shard_common_pos_end;
      (shard_common_pos_start, shard_common_pos_end) = get_pos(x_is_shard_patient_pos[d], 2);
      
      x_is[d, shard_common_pos_start:shard_common_pos_end] = shard_common_pos;
    }
    
    // Fill the arrays for each patient
    for (i in 1:n_patients) {
      x_is[patient_shard[i], 2] += 1; // increment patients in shard count
      
      int visit_start, visit_end;
      (visit_start, visit_end) = get_pos(visit_pos, i);
      
      int visit_m1_start, visit_m1_end;
      (visit_m1_start, visit_m1_end) = get_pos(visit_m1_pos, i);
      
      int n_visits = get_pos_size(visit_pos, i);
      int n_visits_m1 = n_visits - 1;
      
      int x_is_patient_offset = get_pos(x_is_shard_patient_pos[patient_shard[i]], shard_patient_idx[i] + 2).1 - 1;
      
      int x_pos_start, x_pos_end;
      (x_pos_start, x_pos_end) = get_offset_pos(x_is_pos[i], 1, x_is_patient_offset);
      
      // Fill x_is with integer data
      x_is[patient_shard[i], x_pos_start:x_pos_end] = x_is_pos[i];
     
      int packed_visit_start, packed_visit_end;
      (packed_visit_start, packed_visit_end) = get_offset_pos(x_is_pos[i], 2, x_is_patient_offset);
      
      if (packed_visit_start > x_is_size || packed_visit_end > x_is_size) {
        fatal_error("packed_visit_start = ", packed_visit_start, ", x_is_pos[i] = ", x_is_pos[i]);
      }
      
      x_is[patient_shard[i], packed_visit_start:packed_visit_end] = get_int_sub_array(t_visits, visit_pos, i); // Visit times
      
      int packed_theta_pos_start, packed_theta_pos_end;
      (packed_theta_pos_start, packed_theta_pos_end) = get_offset_pos(x_is_pos[i], 3, x_is_patient_offset);
      
      x_is[patient_shard[i], packed_theta_pos_start:packed_theta_pos_end] = theta_pos[i];
      
      int theta_patient_offset = get_pos(theta_shard_patient_pos[patient_shard[i]], shard_patient_idx[i]).1 - 1;
      
      // Fill thetas with patient-specific parameters
      thetas[patient_shard[i], theta_patient_offset + 1] = rho[i]; // GP length scale parameter
      thetas[patient_shard[i], (2 + theta_patient_offset):(3 + theta_patient_offset)] = process_sd; // Process noise std
      thetas[patient_shard[i], theta_patient_offset + 4] = L_process_corr[2, 1]; // Correlation cholesky factor
      
      int proc_noise_start, proc_noise_end;
      (proc_noise_start, proc_noise_end) = get_offset_pos(theta_pos[i], 4, theta_patient_offset);
      
      // Extract raw process noise for this patient (reshape from matrix)
      // Vectorized approach using pre-calculated position indices
      thetas[patient_shard[i], proc_noise_start:proc_noise_end] = to_vector(get_sub_vert_matrix(raw_process_noise, visit_m1_pos, i));
      
      int init_start, init_end;
      (init_start, init_end) = get_offset_pos(theta_pos[i], 5, theta_patient_offset);
      
      // Initial state - use position utility functions for consistent access
      thetas[patient_shard[i], init_start:init_end] = initial_states[i]';
      
      int rates_start, rates_end;
      (rates_start, rates_end) = get_offset_pos(theta_pos[i], 6, theta_patient_offset);
      
      thetas[patient_shard[i], rates_start:rates_end] = [ decrease_rate[i], growth_rate[i], growth_lag[i], growth_transition_rate ]';
    }
    
    // Call map_rect to process patients in parallel
    states = to_matrix(map_rect(calc_patient_states_rect, phi, thetas, rep_array({ delta }, n_shards), x_is), size(t_visits), 2, 0);
  }

  return states;  
}

matrix calc_states(
  data array[] int visit_pos, data array[] int t_visits, vector rho, data real delta, vector process_sd, matrix L_process_corr, 
  matrix raw_process_noise, data int lognormal_noise, data int independ_long_process_noise, data int independ_cross_process_noise, 
  matrix initial_states, vector decrease_rate, vector growth_rate, 
  vector growth_lag, real growth_transition_rate
) {
  return calc_states(
    visit_pos, t_visits, rho, delta, process_sd, L_process_corr, raw_process_noise, lognormal_noise, independ_long_process_noise, independ_cross_process_noise, initial_states,
    decrease_rate, growth_rate, growth_lag, growth_transition_rate, 1, 0
  );
}


matrix calc_states(
  data array[] int visit_pos, data array[] int t_visits, vector rho, data real delta, vector process_sd, matrix L_process_corr, 
  data int independ_long_process_noise, data int independ_cross_process_noise, 
  matrix initial_states, vector decrease_rate, vector growth_rate, 
  vector growth_lag, real growth_transition_rate, data int n_shards, data int debug
) {
  int n_patients = size(visit_pos) - 1;
  
  return calc_states(
    visit_pos, t_visits, rho, delta, process_sd, L_process_corr, rep_matrix(0, size(t_visits) - n_patients, 2), 0, independ_long_process_noise, independ_cross_process_noise, initial_states,
    decrease_rate, growth_rate, growth_lag, growth_transition_rate, n_shards, debug
  );
}

/**
 * Stan function to compute states for multiple patients in a shard (used with map_rect)
 * 
 * @param phi Shared parameters (process_sd)
 * @param theta Shard-specific parameters for all patients in the shard
 * @param x_r Real data (delta)
 * @param x_i Integer data packed as: [n_patients, x_is_pos, theta_pos, debug, lognormal_noise, patient_data...]
 * @return Vector of computed states for all patients in this shard
 */
vector calc_patient_states_rect(vector phi, vector theta, data array[] real x_r, data array[] int x_i) {
  // Extract common metadata from the beginning of x_i
  int shard_id = x_i[1];
  int n_patients_in_shard = x_i[2];
  int max_patients_in_shard = x_i[3];
  int shard_common_size = x_i[4];
  
  int shard_common_pos_size = x_i[shard_common_size]; 
  
  array[shard_common_pos_size + 1] int shard_common_pos = validate_pos(segment(x_i, shard_common_size + 1, shard_common_pos_size + 1));
  
  // Extract position arrays from common section based on shard_common_pos structure
  array[max_patients_in_shard + 3] int x_is_patient_pos = validate_pos(get_int_sub_array(x_i, shard_common_pos, 5)); 
  
  // theta_patient_pos has n_patients + 1 elements
  array[max_patients_in_shard + 1] int theta_patient_pos = validate_pos(get_int_sub_array(x_i, shard_common_pos, 6));
  
  // Get total state size and flags from end of common section
  int total_state_size = x_i[shard_common_pos[7]];
  int debug = x_i[shard_common_pos[8]]; 
  int lognormal_noise = x_i[shard_common_pos[9]];
  int independ_long_process_noise = x_i[shard_common_pos[10]];
  int independ_cross_process_noise = x_i[shard_common_pos[10] + 1];
  int x_pos_size = x_i[shard_common_pos[11]];
  int theta_pos_size = x_i[shard_common_pos[12]];
  
  // Parse real data
  real delta = x_r[1]; // Numerical stability factor
  
  // Process noise std from shared parameters
  vector[2] process_sd = phi;
  
  // Pre-allocate results vector using the pre-calculated size
  vector[total_state_size * 2] all_states;
  int state_offset = 0;
  
  // Process each patient in the shard
  for (p in 1:n_patients_in_shard) {
    // Get patient-specific data positions
    int x_i_patient_start, x_i_patient_end;
    (x_i_patient_start, x_i_patient_end) = get_pos(x_is_patient_pos, p + 2);  // +2 to skip common section and its pos
    
    int theta_patient_start, theta_patient_end;
    (theta_patient_start, theta_patient_end) = get_pos(theta_patient_pos, p);
    
    array[x_pos_size + 1] int x_i_pos = segment(x_i, x_i_patient_start, x_pos_size + 1);
    
    // Get visit information
    int n_visits = get_pos_size(x_i_pos, 2);
    int n_visits_m1 = n_visits - 1;
    
    // Extract time points
    int visits_start, visits_end;
    (visits_start, visits_end) = get_offset_pos(x_i_pos, 2, x_i_patient_start - 1);
    array[n_visits] int time_points = x_i[visits_start:visits_end];
    
    // Extract theta positions
    int theta_pos_start, theta_pos_end;
    (theta_pos_start, theta_pos_end) = get_offset_pos(x_i_pos, 3, x_i_patient_start - 1);
    array[theta_pos_size + 1] int theta_pos = x_i[theta_pos_start:theta_pos_end];
    
    // Extract patient parameters from theta
    real rho = theta[theta_patient_start];
    
    // Process noise std already in phi
    // Skip the process_sd values in theta (positions 2-3)
    
    // Correlation matrix
    matrix[2, 2] L_process_corr = rep_matrix(0, 2, 2);
    L_process_corr[1, 1] = 1.0;
    L_process_corr[2, 1] = theta[theta_patient_start + 3];
    L_process_corr[2, 2] = sqrt(1 - square(L_process_corr[2, 1]));
    
    // Extract process noise
    int noise_start, noise_end;
    (noise_start, noise_end) = get_offset_pos(theta_pos, 4, theta_patient_start - 1);
    matrix[n_visits_m1, 2] raw_process_noise = to_matrix(theta[noise_start:noise_end], n_visits_m1, 2);
    
    // Extract initial states
    int init_start, init_end;
    (init_start, init_end) = get_offset_pos(theta_pos, 5, theta_patient_start - 1);
    row_vector[2] initial_states = theta[init_start:init_end]';
    
    // Extract rates
    int rates_start, rates_end;
    (rates_start, rates_end) = get_offset_pos(theta_pos, 6, theta_patient_start - 1);
    vector[4] rates = theta[rates_start:rates_end];
    
    real decrease_rate = rates[1];
    real growth_rate = rates[2];
    real growth_lag = rates[3];
    real growth_transition_rate = rates[4];
    
    // Calculate states for this patient
    matrix[n_visits, 2] expected_states, states;
    
    (expected_states, states) = calc_patient_states(
      initial_states, time_points, decrease_rate, growth_rate, growth_lag, growth_transition_rate,
      raw_process_noise, lognormal_noise, rho, delta, process_sd, L_process_corr, 
      independ_long_process_noise, independ_cross_process_noise
    );
    
    // Append this patient's states to results
    vector[n_visits * 2] patient_states_vec = to_vector(states');
    all_states[(state_offset + 1):(state_offset + n_visits * 2)] = patient_states_vec;
    state_offset += n_visits * 2;
  }
  
  return all_states;
}

tuple(matrix, matrix) calc_patient_states(
  row_vector initial_state, array[] int time_points, real decrease_rate, real growth_rate, real growth_lag, real growth_transition_rate,
  matrix raw_process_noise, int lognormal_noise, real rho, real delta, vector process_sd, matrix L_process_corr, int independ_long_process_noise, int independ_cross_process_noise
) {
  return calc_patient_states(
    initial_state, time_points, decrease_rate, growth_rate, growth_lag, growth_transition_rate, raw_process_noise, lognormal_noise, rho, delta, process_sd, L_process_corr,
    independ_long_process_noise, independ_cross_process_noise, 0
  );
}

tuple(matrix, matrix) calc_patient_states(
  row_vector initial_state, array[] int time_points, real decrease_rate, real growth_rate, real growth_lag, real growth_transition_rate,
  real rho, real delta, vector process_sd, matrix L_process_corr, int independ_long_process_noise, int independ_cross_process_noise, int debug
) {
  int n_visits = size(time_points);
  int n_visits_m1 = n_visits - 1;
  
  return sf_log_space_trajectory_ncp(initial_state, time_points, decrease_rate, growth_rate, growth_lag, growth_transition_rate, rep_matrix(0, n_visits_m1, 2), debug);
}
  
tuple(matrix, matrix) calc_patient_states(
  row_vector initial_state, array[] int time_points, real decrease_rate, real growth_rate, real growth_lag, real growth_transition_rate,
  matrix raw_process_noise, int lognormal_noise, real rho, real delta, vector process_sd, matrix L_process_corr, 
  int independ_long_process_noise, int independ_cross_process_noise, int debug
) {
  int n_visits = size(time_points);
  int n_visits_m1 = n_visits - 1;
  
  matrix[n_visits_m1, 2] process_noise = calc_patient_process_noise(
    raw_process_noise, lognormal_noise, time_points, rho, delta, process_sd, L_process_corr, independ_long_process_noise, independ_cross_process_noise
  );
  
  return sf_log_space_trajectory_ncp(initial_state, time_points, decrease_rate, growth_rate, growth_lag, growth_transition_rate, process_noise, debug);
}

matrix calc_patient_process_noise(
  matrix raw_process_noise, int lognormal_noise, array[] int time_points, real rho, real delta, vector process_sd, matrix L_process_corr, 
  int independ_long_process_noise, int independ_cross_process_noise
) {
  int n_visits_m1 = size(time_points) - 1;
  matrix[n_visits_m1, 2] noise; 
  
  if (!independ_long_process_noise) {
    // Use Gaussian process to generate correlated noise with absolute time points
    noise = calc_gp_pred(time_points[2:], rho, delta, process_sd, L_process_corr, raw_process_noise, 1);
  } else {
    if (independ_cross_process_noise) {
      noise = raw_process_noise;
    } else {
      // Apply correlation between components
      noise = raw_process_noise * L_process_corr'; 
    }
    
    // Scale process noise with sqrt(delta_t)
    // Vectorized calculation of time differences
    vector[n_visits_m1] delta_t = to_vector(time_points[2:]) - to_vector(time_points[:n_visits_m1]);
    
    // Create scaling matrix with vectorized operations - single line using outer product
    matrix[n_visits_m1, 2] scaling_matrix = sqrt(delta_t) * process_sd';
    
    // Apply scaling with element-wise multiplication
    noise = noise .* scaling_matrix;
  }
 
  if (lognormal_noise) { 
    noise[, 1] = -exp(noise[, 1]);
    noise[, 2] = exp(noise[, 2]);
  }
  
  return noise;
}

/**
 * Generate patient states including process noise, SLD trajectories, and forecast states
 */
tuple(
  matrix,  // forecast_patient_states for this patient
  vector,  // rep_patient_log_sld for this patient
  vector   // forecast_patient_log_sld for this patient
) generate_patient_states_rng(
  matrix patient_states,              // states for this patient's visits
  array[] real forecast_time,         // forecast time points (should be real, not int)
  real patient_log_decrease_rate,     // patient_log_decrease_rate
  real patient_log_growth_rate,       // patient_log_growth_rate
  real sum_tumor_size_baseline,       // baseline tumor size
  // Forecast configuration parameters
  real forecast_growth_lag,           // growth lag for forecast (was hardcoded to negative_infinity())
  real forecast_growth_transition,    // growth transition for forecast (was hardcoded to 1)
  matrix forecast_process_noise,      // process noise for forecast (was generated internally)
  // Global parameters
  real measure_sd
) {
  int n_patient_visits = rows(patient_states);
  // Calculate replicated SLD for observed visits
  vector[n_patient_visits] rep_log_sld = zeros_vector(n_patient_visits);
  
  // First visit uses baseline measurement
  rep_log_sld[1] = log(sum_tumor_size_baseline);
  
  // Subsequent visits are generated from states with measurement noise
  if (n_patient_visits > 1) {
    rep_log_sld[2:] = to_vector(normal_rng(
      to_vector(log_sum_exp(patient_states[2:, 1], patient_states[2:, 2])) + 
        log(sum_tumor_size_baseline),
      rep_vector(measure_sd, n_patient_visits - 1)
    ));
  }

  int forecast_size = size(forecast_time) - 1;
  
  // Calculate forecast states using provided parameters
  matrix[forecast_size, 2] forecast_expected_states, forecast_states;
  (forecast_expected_states, forecast_states) = sf_log_space_trajectory_ncp(
    patient_states[n_patient_visits],
    forecast_time,
    exp(patient_log_decrease_rate), exp(patient_log_growth_rate),
    forecast_growth_lag, forecast_growth_transition, // Use parameters instead of hardcoded values
    forecast_process_noise, // Use provided process noise
    0 // No debug
  );
  
  // Calculate forecast SLD with measurement noise
  vector[forecast_size] forecast_log_sld = zeros_vector(forecast_size);
  if (forecast_size > 0) {
    forecast_log_sld = to_vector(normal_rng(
      to_vector(log_sum_exp(forecast_states[, 1], forecast_states[, 2])) + log(sum_tumor_size_baseline),
      rep_vector(measure_sd, forecast_size)
    ));
  }
  
  return (forecast_states, rep_log_sld, forecast_log_sld);
}

matrix multi_normal_rng(
  matrix y_obs,                 // Observed values [n_obs, 2]
  array[] int time_obs,        // Observed time points
  array[] int time_pred,       // Prediction time points
  real time_rho,                // Temporal length scale
  vector process_sd,            // Process SDs [2]
  matrix L_process_corr,        // Cholesky of process correlation [2, 2]
  real delta                    // Small value for numerical stability
) {
  int n_obs = size(time_obs);
  int n_pred = size(time_pred);
  
  // Calculate temporal covariance matrices
  matrix[n_obs, n_obs] L_K_obs_obs = gp_exp_quad_cholesky_cov(time_obs, 1.0, time_rho, delta);
  matrix[n_pred, n_obs] K_pred_obs = gp_exp_quad_cov(time_pred, time_obs, 1.0, time_rho);
  matrix[n_pred, n_pred] K_pred_pred = gp_exp_quad_cov(time_pred, 1.0, time_rho, delta);
  
  // Create conditional mean matrix
  matrix[n_pred, 2] mu_cond;
  mu_cond[, 1] = gp_conditional_mean(y_obs[, 1], L_K_obs_obs, K_pred_obs); 
  mu_cond[, 2] = gp_conditional_mean(y_obs[, 2], L_K_obs_obs, K_pred_obs); 
  
  matrix[n_pred, n_pred] K_cond = gp_conditional_cov(L_K_obs_obs, K_pred_obs, K_pred_pred, delta);
  matrix[n_pred, n_pred] L_K_cond = cholesky_decompose(K_cond);
  
  // Generate standard normal random values
  matrix[n_pred, 2] eta_raw = to_matrix(to_vector(normal_rng(zeros_vector(n_pred * 2), rep_vector(1, n_pred * 2))), n_pred, 2);
  
  // Create the sample using the separable structure
  return mu_cond + L_K_cond * eta_raw * diag_pre_multiply(process_sd, L_process_corr)';
}

matrix multi_normal_rng(
  int n_pred,
  vector process_sd,           // Process SDs [2]
  matrix L_process_corr        // Cholesky of process correlation [2, 2]
) {
  // Generate standard normal random values
  matrix[n_pred, 2] eta_raw = to_matrix(to_vector(normal_rng(zeros_vector(n_pred * 2), rep_vector(1, n_pred * 2))), n_pred, 2);
  
  // Create the sample using the separable structure
  return eta_raw * diag_pre_multiply(process_sd, L_process_corr)';
}

void assert_matching_states(
  matrix states, row_vector initial_states, array[] int time_points, real decrease_rate, real growth_rate, real growth_lag, real growth_transit_rate,
  matrix process_noise, int debug 
) {
  int n_visits = size(time_points);
  
  matrix[n_visits, 2] curr_obs_states = sf_log_space_trajectory_ncp(
      initial_states,
      time_points,
      decrease_rate, growth_rate,
      growth_lag, growth_transit_rate, 
      process_noise,
      debug
    ).2;
    
  matrix[n_visits, 2] rate_diff = states - curr_obs_states;

  for (t in 1:n_visits) {
    int dec = abs(rate_diff[t, 1]) > 1e-6;
    int gro = abs(rate_diff[t, 2]) > 1e-6;

    if (dec || gro) {
      // print("initial_state = ", [ patient_log_decrease_prop[i], patient_log_growth_prop[i] ], ", time_points = ", get_int_sub_array(t_patient_visits, patient_visit_pos, i),
      //       ", dec rate = ", exp(patient_log_decrease_rate[i]), ", gro rate = ", exp(patient_log_growth_rate[i]), ", lag = ", exp(patient_log_growth_lag[i]), ", transit = ", exp(pop_log_growth_transition_rate));

      fatal_error("t = ", t, ", dec = ", dec, ", gro = ", gro, ", states[t, 1] = ", states[t, 1], ", curr_obs_states[t, 1] = ", curr_obs_states[t, 1],
      ", states[t, 2] = ", states[t, 2], ", curr_obs_states[t, 2] = ", curr_obs_states[t, 2]);
    }
  }
}