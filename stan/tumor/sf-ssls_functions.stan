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

/* Parallel computation of states for multiple patients using map_rect
 * 
 * @param visit_pos Position array for patient visits
 * @param t_visits Array of time points for all patients
 * @param rho Vector of GP length scale parameters for each patient
 * @param delta Small value for numerical stability
 * @param process_sd Process noise standard deviations [2]
 * @param L_process_corr Cholesky factor of process correlation matrix
 * @param raw_process_noise Raw process noise values
 * @param independ_long_process_noise Flag for independent longitudinal process noise
 * @param independ_cross_process_noise Flag for independent cross-component process noise
 * @param initial_states Matrix of initial states for all patients
 * @param decrease_rate Vector of tumor decrease rates
 * @param growth_rate Vector of tumor growth rates
 * @param growth_lag Vector of growth lag parameters
 * @param growth_transition_rate Growth transition rate parameter
 * @param debug Print debug information
 * @return Vector of computed states for all patients
 */
matrix calc_states(
  data array[] int visit_pos, data array[] int t_visits, vector rho, data real delta, vector process_sd, matrix L_process_corr, 
  matrix raw_process_noise, data int independ_long_process_noise, data int independ_cross_process_noise, 
  matrix initial_states, vector decrease_rate, vector growth_rate, 
  vector growth_lag, real growth_transition_rate, int parallel, int debug
) {
  int n_patients = size(visit_pos) - 1;
  // Create position array for time points with one fewer elements per patient
  array[n_patients + 1] int visit_m1_pos = create_pos(get_pos_size(visit_pos), -1);  
  
  assert_equal(visit_pos[n_patients + 1], size(t_visits) + 1);
  
  // Arrays to track indices for various components
  array[n_patients] int visit_end_idx;
  array[n_patients] int noise_end_idx;
  array[n_patients] int initial_states_idx;
  array[n_patients] int rates_idx;
  
  int theta_pos_size = 6;
  array[n_patients, theta_pos_size + 1] int theta_pos;
  
  int x_pos_size = 5;
  array[n_patients, x_pos_size + 1] int x_is_pos;
  
  int max_x_is_size = 0, max_theta_size = 0;
  
  // Calculate indices for each patient based on their number of visits
  for (i in 1:n_patients) {
    int n_visits = get_pos_size(visit_pos, i);
   
    x_is_pos[i] = create_pos({
      x_pos_size + 1, // the x_is_pos 
      2,              // independ flags  
      get_pos_size(visit_pos, i), // time points
      theta_pos_size + 1,   // theta positions
      1                     // debug flag
    });
    
    theta_pos[i] = create_pos({ 
      1, // rho (length scale)
      2, // process_sd
      1, // L_process_corr[2, 1] 
      (n_visits - 1) * 2, // raw_process_noise for this patient
      2, // initial states 
      4 
    }); // rates (2), lag, transition
    
    max_x_is_size = max(max_x_is_size, get_pos_total_size(x_is_pos[i]));
    max_theta_size = max(max_theta_size, get_pos_total_size(theta_pos[i]));
  }
  
  vector[2] phi = process_sd; // Shared parameters vector
  
  // Initialize parameter arrays for map_rect
  array[n_patients] vector[max_theta_size] thetas = rep_array(rep_vector(negative_infinity(), max_theta_size), n_patients);
  array[n_patients, max_x_is_size] int x_is = rep_array(-1111, n_patients, max_x_is_size);
  
  // Fill the arrays for each patient
  for (i in 1:n_patients) {
    int visit_start, visit_end;
    (visit_start, visit_end) = get_pos(visit_pos, i);
    
    int visit_m1_start, visit_m1_end;
    (visit_m1_start, visit_m1_end) = get_pos(visit_m1_pos, i);
    
    int n_visits = get_pos_size(visit_pos, i);
    int n_visits_m1 = n_visits - 1;
    
    int x_pos_start, x_pos_end;
    (x_pos_start, x_pos_end) = get_pos(x_is_pos[i], 1);
    
    // Fill x_is with integer data
    x_is[i, x_pos_start:x_pos_end] = x_is_pos[i];
   
    int flags_start, flags_end;
    (flags_start, flags_end) = get_pos(x_is_pos[i], 2);
    
    x_is[i, flags_start:flags_end] = { independ_long_process_noise, independ_cross_process_noise };
    
    int packed_visit_start, packed_visit_end;
    (packed_visit_start, packed_visit_end) = get_pos(x_is_pos[i], 3);
    
    x_is[i, packed_visit_start:packed_visit_end] = get_int_sub_array(t_visits, visit_pos, i); // Visit times
    
    int packed_theta_pos_start, packed_theta_pos_end;
    (packed_theta_pos_start, packed_theta_pos_end) = get_pos(x_is_pos[i], 4);
    
    x_is[i, packed_theta_pos_start:packed_theta_pos_end] = theta_pos[i];
    x_is[i, packed_theta_pos_end + 1] = i == 1 ? debug : 0; 
    
    // Fill thetas with patient-specific parameters
    thetas[i, 1] = rho[i]; // GP length scale parameter
    thetas[i, 2:3] = process_sd; // Process noise std
    thetas[i, 4] = L_process_corr[2, 1]; // Correlation cholesky factor
    
    int proc_noise_start, proc_noise_end;
    (proc_noise_start, proc_noise_end) = get_pos(theta_pos[i], 4);
    
    // Extract raw process noise for this patient (reshape from matrix)
    // Vectorized approach using pre-calculated position indices
    thetas[i, proc_noise_start:proc_noise_end] = to_vector(get_sub_vert_matrix(raw_process_noise, visit_m1_pos, i));
    
    int init_start, init_end;
    (init_start, init_end) = get_pos(theta_pos[i], 5);
    
    // Initial state - use position utility functions for consistent access
    thetas[i, init_start:init_end] = initial_states[i]';
    
    int rates_start, rates_end;
    (rates_start, rates_end) = get_pos(theta_pos[i], 6);
    
    thetas[i, rates_start:rates_end] = [ decrease_rate[i], growth_rate[i], growth_lag[i], growth_transition_rate ]';
    
    if (theta_pos[i, theta_pos_size + 1] - 1 > max_theta_size) {
      fatal_error(i, ": pos exceeds expected max size for theta: theta_pos[i, <end>] - 1 = ", theta_pos[i, theta_pos_size + 1] - 1, ", max_theta_size = ", max_theta_size);
    }
  }
  
  matrix[size(t_visits), 2] states;
  
  if (parallel) {
    // Call map_rect to process patients in parallel
    states = to_matrix(map_rect(calc_patient_states, phi, thetas, rep_array({ delta }, n_patients), x_is), size(t_visits), 2, 0);
  } else {
    for (i in 1:n_patients) {
      int visit_start, visit_end, n_visits;
      (visit_start, visit_end) = get_pos(visit_pos, i);
      n_visits = get_pos_size(visit_pos, i);
      
      states[visit_start:visit_end] = to_matrix(calc_patient_states(phi, thetas[i], { delta }, x_is[i]), n_visits, 2, 0);
    }
  }

  return states;  
}

matrix calc_states(
  data array[] int visit_pos, data array[] int t_visits, vector rho, data real delta, vector process_sd, matrix L_process_corr, 
  matrix raw_process_noise, data int independ_long_process_noise, data int independ_cross_process_noise, 
  matrix initial_states, vector decrease_rate, vector growth_rate, 
  vector growth_lag, real growth_transition_rate
) {
  return calc_states(
    visit_pos, t_visits, rho, delta, process_sd, L_process_corr, raw_process_noise, independ_long_process_noise, independ_cross_process_noise, initial_states,
    decrease_rate, growth_rate, growth_lag, growth_transition_rate, 1, 0
  );
}

/**
 * Stan function to compute states for a single patient (used with map_rect)
 * 
 * @param phi Shared parameters (process_sd)
 * @param theta Patient-specific parameters
 * @param x_r Real data (delta)
 * @param x_i Integer data (visit times, flags)
 * @return Vector of computed states for this patient
 */
vector calc_patient_states(vector phi, vector theta, data array[] real x_r, data array[] int x_i) {
  int x_pos_size = 5;
  array[x_pos_size + 1] int x_i_pos = x_i[:(x_pos_size + 1)];
  
  int flags_start, flags_end;
  (flags_start, flags_end) = get_pos(x_i_pos, 2);
  
  int independ_long_process_noise = x_i[flags_start];
  int independ_cross_process_noise = x_i[flags_end];
  
  // Parse integer data
  int n_visits = get_pos_size(x_i_pos, 3); // Number of visits
  int n_visits_m1 = n_visits - 1;
  
  array[n_visits] int time_points = get_int_sub_array(x_i, x_i_pos, 3); 
 
  int theta_pos_size = get_pos_size(x_i_pos, 4) - 1;
  array[theta_pos_size + 1] int theta_pos = get_int_sub_array(x_i, x_i_pos, 4);
 
  int debug_start = get_pos(x_i_pos, 5).1; 
  int debug = x_i[debug_start];
 
  // Parse real data
  real delta = x_r[1]; // Numerical stability factor
 
  // Extract parameters from theta
  real rho = theta[1]; // GP length scale
  vector[2] process_sd = phi; // Process noise std (from shared parameters)
  
  // Reconstruct Cholesky factor of correlation matrix
  matrix[2, 2] L_process_corr = rep_matrix(0, 2, 2);
  L_process_corr[1, 1] = 1.0;
  L_process_corr[2, 1] = theta[4];
  L_process_corr[2, 2] = sqrt(1 - square(L_process_corr[2, 1])); 
  
  // Extract raw process noise and construct matrix
  matrix[n_visits_m1, 2] raw_process_noise = to_matrix(get_sub_vector(theta, theta_pos, 4), n_visits_m1, 2);
  matrix[n_visits_m1, 2] process_noise = calc_patient_process_noise(
    raw_process_noise, time_points, rho, delta, process_sd, L_process_corr, independ_long_process_noise, independ_cross_process_noise
  );
  
  // Extract initial state
  row_vector[2] initial_state = get_sub_row_vector(theta, theta_pos, 5);
 
  vector[4] rates = get_sub_vector(theta, theta_pos, 6); 
  // Extract rates and other parameters
  real decrease_rate = rates[1];
  real growth_rate = rates[2];
  real growth_lag = rates[3];
  real growth_transition_rate = rates[4];
    
  // Calculate states using Stein-Fojo log-space state space model
  matrix[n_visits, 2] expected_states;
  matrix[n_visits, 2] states;
  
  (expected_states, states) = sf_log_space_trajectory_ncp(
    initial_state, time_points,
    decrease_rate, growth_rate, growth_lag, growth_transition_rate,
    process_noise
  );
  
  
  if (debug) {
    print("initial_state = ", initial_state, ", time_points = ", time_points, ", decrease_rate = ", decrease_rate, ", growth_rate = ", growth_rate, 
          ", growth_lag = ", growth_lag, ", growth_transition_rate = ", growth_transition_rate);
    print("expected_states = ", expected_states, ", states = ", states);
  }
  
  // Convert to vector for map_rect output - vectorized approach
  // to_vector converts matrix to column-major vector
  return to_vector(states');
}

matrix calc_patient_process_noise(
  matrix raw_process_noise, array[] int time_points, real rho, real delta, vector process_sd, matrix L_process_corr, int independ_long_process_noise, int independ_cross_process_noise
) {
  if (!independ_long_process_noise) {
    // Use Gaussian process to generate correlated noise with absolute time points
    return calc_gp_pred(time_points[2:], rho, delta, process_sd, L_process_corr, raw_process_noise, 1);
  } else {
    int n_visits_m1 = size(time_points) - 1;
    matrix[n_visits_m1, 2] process_noise;
    
    if (independ_cross_process_noise) {
      process_noise = raw_process_noise;
    } else {
      // Apply correlation between components
      process_noise = raw_process_noise * L_process_corr'; 
    }
    
    // Scale process noise with sqrt(delta_t)
    // Vectorized calculation of time differences
    vector[n_visits_m1] delta_t = to_vector(time_points[2:]) - to_vector(time_points[:n_visits_m1]);
    
    // Create scaling matrix with vectorized operations - single line using outer product
    matrix[n_visits_m1, 2] scaling_matrix = sqrt(delta_t) * process_sd';
    
    // Apply scaling with element-wise multiplication
    return process_noise .* scaling_matrix;
  }
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
  
  // // Process the first dimension
  // (mu_cond[,1], K_cond) = gp_conditional(y_obs[,1], K_obs_obs, K_pred_obs, K_pred_pred, delta);
  // 
  // // Process the second dimension (reusing the same K matrices)
  // mu_cond[,2] = gp_conditional(y_obs[,2], K_obs_obs, K_pred_obs, K_pred_pred, delta).1;
  // 
  // Get Cholesky of temporal covariance
  matrix[n_pred, n_pred] L_K_cond = cholesky_decompose(K_cond);
  
  // Generate standard normal random values
  matrix[n_pred, 2] eta_raw = to_matrix(to_vector(normal_rng(zeros_vector(n_pred * 2), rep_vector(1, n_pred * 2))), n_pred, 2);
  
  // Create the sample using the separable structure
  matrix[n_pred, 2] sample = mu_cond + L_K_cond * eta_raw * diag_pre_multiply(process_sd, L_process_corr)';
  
  return sample;
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