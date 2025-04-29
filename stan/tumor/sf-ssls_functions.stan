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

tuple(row_vector, row_vector) sf_log_space_transition_ncp(row_vector raw_next_x, row_vector current_x, 
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
  
  vector[T] time_varying_factor = get_growth_lag_factor(times, growth_lag, transition_rate); //  inv_logit((times[t] - growth_lag) / transition_rate);
  
  // State transitions
  for (t in 2:T) {
    (expected_x[t], x[t]) = sf_log_space_transition_ncp(raw_x[t - 1], x[t - 1], times[t], times[t - 1],
                                                        decrease_rate, time_varying_factor[t] * growth_rate,
                                                        process_sd, L_process_corr);
  }
  
  return (expected_x, x);
}

tuple(matrix, matrix) sf_log_space_trajectory_ncp(
  matrix raw_x, row_vector x0, array[] real times,
  real decrease_rate, real growth_rate, real growth_lag, real transition_rate,
  matrix process_noise
) {
  int T = rows(raw_x) + 1;
  matrix[T, 2] x;
  matrix[T, 2] expected_x;
  x[1] = x0;
  expected_x[1] = x0;
  
  vector[T] time_varying_factor = get_growth_lag_factor(times, growth_lag, transition_rate); //  inv_logit((times[t] - growth_lag) / transition_rate);
  
  // State transitions
  for (t in 2:T) {
    (expected_x[t], x[t]) = sf_log_space_transition_ncp(raw_x[t - 1], x[t - 1], times[t], times[t - 1],
                                                        decrease_rate, time_varying_factor[t] * growth_rate,
                                                        process_noise[t - 1]);
  }
  
  return (expected_x, x);
}