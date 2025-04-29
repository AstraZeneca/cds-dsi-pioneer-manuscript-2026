// Helper functions for numerical stability and code reuse

/**
 * Ensure a matrix is symmetric
 * 
 * @param cov Matrix to symmetrize
 * @return Symmetrized matrix
 */
matrix ensure_symmetry(matrix cov) {
  return 0.5 * (cov + cov');
}

/**
 * Ensure a matrix is positive definite by setting minimum values on diagonal
 * 
 * @param cov Matrix to ensure positive definiteness
 * @param epsilon Small positive value to use as minimum diagonal value
 * @return Modified matrix with positive diagonal elements
 */
matrix ensure_positive_definite(matrix cov, real epsilon) {
  matrix[rows(cov), cols(cov)] result = cov;
  for (i in 1:rows(cov)) {
    if (result[i, i] <= 0) {
      result[i, i] = epsilon;
    }
  }
  return result;
}

/**
 * Compute conditional mean for backward sampling
 * 
 * @param filtered_mean Current filtered state mean
 * @param next_sample Next state sample in backward pass
 * @param F State transition matrix
 * @param smoothing_gain Smoothing gain matrix
 * @return Conditional mean vector
 */
vector compute_conditional_mean(vector filtered_mean, vector next_sample, matrix F, matrix smoothing_gain) {
  return filtered_mean + smoothing_gain * (next_sample - F * filtered_mean);
}

/**
 * Compute smoothing gain matrix
 * 
 * @param P_t Current state covariance
 * @param F State transition matrix
 * @param P_pred Predicted state covariance
 * @return Smoothing gain matrix
 */
matrix compute_smoothing_gain(matrix P_t, matrix F, matrix P_pred) {
  // Ensure numerical symmetry
  matrix[rows(P_pred), cols(P_pred)] sym_P_pred = ensure_symmetry(P_pred);
  
  // Calculate smoothing gain more stably using solver instead of direct inversion
  matrix[rows(P_t), cols(P_t)] C = mdivide_left_spd(sym_P_pred, F * P_t');
  return P_t * F' * C';
}

// Create state transition matrix using proper matrix exponentiation
matrix sf_create_log_F(real decrease_rate, real grow_rate, real delta_t) {
  // Create the coefficient matrix for the ODE system
  matrix[2, 2] A;
  A[1, 1] = -decrease_rate;
  A[1, 2] = 0;
  A[2, 1] = 0;
  A[2, 2] = grow_rate;
  
  return A * delta_t; 
}

matrix sf_create_log_F(real decrease_rate, real grow_rate) {
  return sf_create_log_F(decrease_rate, grow_rate, 1.0);
}

// Create state transition matrix using proper matrix exponentiation
matrix sf_create_F(real decrease_rate, real grow_rate, real delta_t) {
  // Calculate matrix exponential: F = exp(A*delta_t)
  return matrix_exp(sf_create_log_F(decrease_rate, grow_rate, delta_t)); 
}

matrix sf_create_F(real decrease_rate, real grow_rate) {
  return sf_create_F(decrease_rate, grow_rate, 1.0);
}

/**
* Forward Kalman filter for the Stein-Fojo (SF) tumor growth model
*
* This function implements a numerically stable Kalman filter for the
* SF model of tumor growth. It returns both filtered means and covariances 
* which can be used for state estimation or as inputs to backward sampling.
*
* @param normalized_y Array of normalized tumor measurements
* @param time_points Array of observation time points
* @param decrease_rate Tumor regression rate constant (d)
* @param growth_rate Tumor growth rate constant (g)
* @param decrease_process_sd Process noise standard deviation for decrease component
* @param growth_process_sd Process noise standard deviation for growth component
* @param measure_sd Measurement noise standard deviation
*
* @return A tuple containing:
*   - array[] vector: filtered state means
*   - array[] matrix: filtered state covariances
*/
tuple(matrix, array[] matrix) sf_forward_filter(
  vector x_0,
  vector normalized_y, 
  array[] int time_points, 
  real decrease_rate, 
  real growth_rate, 
  real decrease_process_sd, 
  real growth_process_sd, 
  real measure_sd,
  int log_space
) {
  real decrease_process_var = decrease_process_sd^2;
  real growth_process_var = growth_process_sd^2;
  real measure_var = measure_sd^2;
  
  int N = size(time_points);
  int W = max(time_points) - min(time_points) + 1;
  matrix[2, W] filtered_means;
  array[W] matrix[2, 2] filtered_covs;
  array[N] int obs_idx;
  
  // Process noise covariance
  matrix[2, 2] Q = diag_matrix([decrease_process_var, growth_process_var]');
  
  // Initialize
  obs_idx[1] = 1;
  filtered_means[, 1] = x_0;
  filtered_covs[1] = Q;
  
  // Measurement matrix
  row_vector[2] H = [1, 1];
  
  int obs_pos = 1;
  
  // Forward pass with storage
  for (t in 2:W) {
    if (obs_pos + 1 > N) {
    }
    
    // print(t, ": time_points = ", time_points, ", W = ", W, ", obs_pos = ", obs_pos);
    
    int is_observed = time_points[obs_pos + 1] - time_points[1] + 1 == t;
    
    // real delta_t =  time_points[t] - time_points[t-1];
    matrix[2, 2] F = log_space ? sf_create_log_F(decrease_rate, growth_rate) : sf_create_F(decrease_rate, growth_rate);
    
    // Predict
    vector[2] pred_mean = F * filtered_means[, t-1];
    matrix[2, 2] pred_cov = F * filtered_covs[t-1] * F' + Q;
    pred_cov = ensure_symmetry(pred_cov);
    
    // Update
    // S should be a scalar or 1×1 matrix
    real S = (H * pred_cov * H') + measure_var + 1e-10;
    // K is calculated with scalar division
    vector[2] K = (pred_cov * H') / S;
    
    filtered_means[, t] = pred_mean;
   
    if (is_observed) {
      obs_pos += 1;
      obs_idx[obs_pos] = t;
      
      real pred_measurement = log_space ? log_sum_exp(pred_mean) : sum(pred_mean);
      real innovation = normalized_y[obs_pos] - pred_measurement;
      filtered_means[, t] += K * innovation;
    } 
    
    // Joseph form for covariance update (more stable)
    matrix[2, 2] I_KH = identity_matrix(2) - K * H;
    filtered_covs[t] = I_KH * pred_cov * I_KH' + K * measure_var * K';
    filtered_covs[t] = ensure_symmetry(filtered_covs[t]);
  }
  
  assert_equal(obs_pos, N);
  
  return (filtered_means[, obs_idx], filtered_covs[obs_idx]);
} 


tuple(matrix, array[] matrix) sf_forward_filter(
  vector x_0,
  vector normalized_y, 
  array[] int time_points, 
  real decrease_rate, 
  real growth_rate, 
  real decrease_process_sd, 
  real growth_process_sd, 
  real measure_sd
) {
  return sf_forward_filter(x_0, normalized_y, time_points, decrease_rate, growth_rate, decrease_process_sd, growth_process_sd, measure_sd, 0);
}

/**
 * Backward sampling for the Stein-Fojo (SF) tumor growth model
 *
 * This function implements the backward sampling step of the Forward-Filtering 
 * Backward-Sampling (FFBS) algorithm for the SF tumor growth model. It draws samples 
 * from the joint smoothing distribution by recursively sampling from conditional 
 * distributions, moving backward in time.
 *
 * @param filtered_means Array of state vector means from the forward filter
 * @param filtered_covs Array of state covariance matrices from the forward filter
 * @param time_points Array of observation time points
 * @param decrease_rate Tumor regression rate constant (d)
 * @param growth_rate Tumor growth rate constant (g)
 * @param decrease_process_sd Process noise standard deviation for decrease component
 * @param growth_process_sd Process noise standard deviation for growth component
 *
 * @return Array of sampled state vectors from the joint smoothing distribution
 */
matrix sf_backward_sample_rng(
  matrix filtered_means, 
  array[] matrix filtered_covs, 
  array[] real time_points, 
  real decrease_rate, 
  real growth_rate, 
  real decrease_process_sd, 
  real growth_process_sd
) {
  real decrease_process_var = decrease_process_sd^2;
  real growth_process_var = growth_process_sd^2;
  
  int N = size(time_points);
  matrix[2, N] smoothed_samples;
  
  // Process noise covariance
  matrix[2, 2] Q = diag_matrix([decrease_process_var, growth_process_var]');
  
  // Sample final state from filtered distribution
  smoothed_samples[, N] = multi_normal_rng(filtered_means[, N], filtered_covs[N]);
  
  // Backward pass
  for (t_idx in 1:(N - 1)) {
    int t = N - t_idx;
    real delta_t = time_points[t+1] - time_points[t];
    matrix[2, 2] F = sf_create_F(decrease_rate, growth_rate, delta_t);
    
    // Compute predicted covariance
    matrix[2, 2] P_pred = F * filtered_covs[t] * F' + Q;
    
    // Calculate smoothing gain
    matrix[2, 2] smoothing_gain = compute_smoothing_gain(filtered_covs[t], F, P_pred);
    
    // Compute conditional mean
    vector[2] cond_mean = compute_conditional_mean(filtered_means[, t], smoothed_samples[, t+1], F, smoothing_gain);
    
    // Compute conditional covariance
    matrix[2, 2] cond_cov = filtered_covs[t] - smoothing_gain * P_pred * smoothing_gain';
    
    // Ensure symmetry and positive definiteness
    cond_cov = ensure_symmetry(cond_cov);
    cond_cov = ensure_positive_definite(cond_cov, 1e-6);
    
    // Sample
    smoothed_samples[, t] = multi_normal_rng(cond_mean, cond_cov);
  }
  
  return smoothed_samples;
}

 /**
 * Forward Filtering Backward Sampling for the Stein-Fojo (SF) tumor growth model
 *
 * This function implements the complete FFBS algorithm for the SF model by:
 * 1. Calling sf_forward_filter() to perform Kalman filtering
 * 2. Using the filtered results to perform backward sampling
 *
 * @param normalized_y Array of normalized tumor measurements
 * @param time_points Array of observation time points
 * @param decrease_rate Tumor regression rate constant (d)
 * @param growth_rate Tumor growth rate constant (g)
 * @param decrease_process_sd Process noise standard deviation for decrease component
 * @param growth_process_sd Process noise standard deviation for growth component
 * @param measure_sd Measurement noise standard deviation
 *
 * @return Array of sampled state vectors from the joint smoothing distribution
 */
matrix sf_forward_filter_backward_sample_rng(
  vector x_0,
  vector normalized_y, 
  array[] int time_points, 
  real decrease_rate, 
  real growth_rate, 
  real decrease_process_sd, 
  real growth_process_sd, 
  real measure_sd
) {
  matrix[2, size(time_points)] filtered_means;
  array[size(time_points)] matrix[2, 2] filtered_covs;
  
  // Call forward filter to get filtered means and covariances
  (filtered_means, filtered_covs) = sf_forward_filter(x_0, normalized_y, time_points, decrease_rate, growth_rate, decrease_process_sd, growth_process_sd, measure_sd);
  
  // Perform backward sampling using filtered results
  return sf_backward_sample_rng(filtered_means, filtered_covs, time_points, decrease_rate, growth_rate, 
                              decrease_process_sd, growth_process_sd);
}  

/**
 * Forward Kalman filter for the Stein-Fojo (SF) tumor growth model using square root form
 *
 * This implementation uses the square root formulation of the Kalman filter
 * for improved numerical stability. It directly propagates Cholesky factors
 * of covariance matrices rather than the covariances themselves.
 *
 * @param normalized_y Vector of normalized tumor measurements
 * @param time_points Vector of observation time points
 * @param decrease_rate Tumor regression rate constant (d)
 * @param growth_rate Tumor growth rate constant (g)
 * @param process_sd Process noise standard deviation (not variance)
 * @param measure_sd Measurement noise standard deviation (not variance)
 *
 * @return A tuple containing:
 *   - matrix[2, N]: filtered state means
 *   - array[N] matrix[2, 2]: filtered state covariance square roots (Cholesky factors)
 */
tuple(matrix, array[] matrix) sf_forward_filter_sqrt(
    vector normalized_y, 
    array[] real time_points, 
    real decrease_rate, 
    real growth_rate, 
    real decrease_process_sd, 
    real growth_process_sd, 
    real measure_sd
) {
  int N = num_elements(time_points);
  matrix[2, N] filtered_means;
  array[N] matrix[2, 2] filtered_sqrt_covs;  // Cholesky factors of covariances
  
  // Initialize
  filtered_means[, 1] = [1.0, 1.0]';
  filtered_sqrt_covs[1] = diag_matrix([decrease_process_sd, growth_process_sd]');
  
  // Measurement matrix
  row_vector[2] H = [1, 1];
  
  // Process noise - store square root directly
  matrix[2, 2] sqrt_Q = diag_matrix([decrease_process_sd, growth_process_sd]');
  
  // Forward pass with square root covariance propagation
  for (t in 2:N) {
    real delta_t = time_points[t] - time_points[t-1];
    matrix[2, 2] F = sf_create_F(decrease_rate, growth_rate, delta_t);
    
    // --- Prediction step using QR decomposition ---
    
    // Form augmented matrix [F*S_{t-1}, sqrt_Q]
    matrix[2, 4] aug_mat = append_col(F * filtered_sqrt_covs[t-1], sqrt_Q);
    
    // QR decomposition for predicted square root covariance
    matrix[2, 2] pred_sqrt_cov = qr_thin_R(aug_mat)[1:2, 1:2];
    
    // Ensure upper triangular (numerical stability)
    for (i in 1:2) {
      for (j in 1:i-1) {
        pred_sqrt_cov[i, j] = 0;
      }
      // Ensure positive diagonal (for proper Cholesky factor)
      if (pred_sqrt_cov[i, i] < 0) {
        pred_sqrt_cov[i, ] = -pred_sqrt_cov[i, ];
      }
    }
    
    // Predict mean normally
    vector[2] pred_mean = F * filtered_means[, t-1];
    
    // --- Update step using square root formulation ---
    
    // Predicted measurement and innovation
    real pred_measurement = pred_mean[1] + pred_mean[2];
    real innovation = normalized_y[t] - pred_measurement;
    
    // Compute S = sqrt(H * P * H' + R)
    // First calculate H * pred_sqrt_cov efficiently
    row_vector[2] HP_sqrt = H * pred_sqrt_cov;
    real S = sqrt(dot_product(HP_sqrt, HP_sqrt) + measure_sd * measure_sd);
    
    // Compute Kalman gain
    vector[2] K = (pred_sqrt_cov * pred_sqrt_cov' * H') / S;
    
    // Update state estimate
    filtered_means[, t] = pred_mean + K * innovation;
    
    // --- Update covariance using square root form ---
    // Using Cholesky update/downdate approach
    
    // Scale innovation vector for rank-1 update
    vector[2] v = K * (measure_sd / S);
    
    // Compute square root of (I - K*H)*P*(I - K*H)' + K*R*K'
    // Using efficient rank-1 downdate (Joseph form square root equivalent)
    matrix[2, 2] U = pred_sqrt_cov;
    
    // Rank-1 downdate via Cholesky update on U'*U matrix
    for (j in 1:2) {
      real s = v[j] / U[j, j];
      real c, r;
      
      // Handle extreme or invalid cases
      if (abs(s) > 1.0) {
        // Case where |v[j]| > |U[j,j]| - not valid for a proper Givens rotation
        // Use regularization to ensure numerical stability
        s = s < 0 ? -0.9999 : 0.9999;  // Ternary operator to preserve sign
        c = sqrt(1.0 - s^2);   // Ensure c^2 + s^2 = 1
        r = U[j, j] * c;       // Derive r from the constraint
      } else {
        // Normal case, use stable computation
        c = sqrt(1.0 - s^2);
        r = U[j, j] * c;
      } 
      
      // Apply Givens rotation to relevant parts of U
      U[j, j] = r;
      
      if (j < 2) {
        row_vector[2-j] z = U[j, (j+1):2]; 
        row_vector[2-j] y = v[(j+1):2]';
        U[j, (j+1):2] = c * z - s * y;
        v[(j+1):2] = v[(j+1):2] - v[j] * U[j, (j+1):2]' / U[j, j];
      }
    }
    
    filtered_sqrt_covs[t] = U;
    
    // Ensure upper triangular form
    for (i in 1:2) {
      for (j in 1:i-1) {
        filtered_sqrt_covs[t][i, j] = 0;
      }
    }
  }
  
  return (filtered_means, filtered_sqrt_covs);
}

/**
 * Backward sampling for the Stein-Fojo (SF) tumor growth model using square root form
 * 
 * This function implements backward sampling using square root forms
 * of covariance matrices. It takes filtered means and covariance 
 * square roots from a forward pass and produces state trajectories 
 * sampled from the full posterior distribution.
 *
 * @param filtered_means Matrix of state means from the forward filter
 * @param filtered_sqrt_covs Array of state covariance square roots from the forward filter
 * @param time_points Vector of observation time points
 * @param decrease_rate Tumor regression rate constant (d)
 * @param growth_rate Tumor growth rate constant (g)
 * @param decrease_process_sd Process noise standard deviation for decrease component
 * @param growth_process_sd Process noise standard deviation for growth component
 *
 * @return Matrix of sampled state vectors from the joint smoothing distribution
 */
matrix sf_backward_sample_sqrt_rng(
    matrix filtered_means,
    array[] matrix filtered_sqrt_covs,
    array[] real time_points,
    real decrease_rate,
    real growth_rate,
    real decrease_process_sd,
    real growth_process_sd
) {
  int N = num_elements(time_points);
  matrix[2, N] smoothed_samples;
  
  // Process noise square root
  matrix[2, 2] sqrt_Q = diag_matrix([decrease_process_sd, growth_process_sd]');
  
  // Sample final state from filtered distribution
  matrix[2, 2] cov_N = filtered_sqrt_covs[N] * filtered_sqrt_covs[N]';
  smoothed_samples[, N] = multi_normal_rng(filtered_means[, N], cov_N);
  
  // Backward sampling pass
  for (t_idx in 1:(N-1)) {
    int t = N - t_idx; 
    real delta_t = time_points[t+1] - time_points[t];
    matrix[2, 2] F = sf_create_F(decrease_rate, growth_rate, delta_t);
    
    // Convert square root covariance to full covariance for this step
    matrix[2, 2] P_t = filtered_sqrt_covs[t] * filtered_sqrt_covs[t]';
    
    // Compute predicted covariance
    matrix[2, 2] P_pred = F * P_t * F' + sqrt_Q * sqrt_Q';
    
    // Calculate smoothing gain
    matrix[2, 2] smoothing_gain = compute_smoothing_gain(P_t, F, P_pred);
    
    // Compute conditional mean
    vector[2] cond_mean = compute_conditional_mean(filtered_means[, t], smoothed_samples[, t+1], F, smoothing_gain);
    
    // Compute conditional covariance
    matrix[2, 2] cond_cov = P_t - smoothing_gain * P_pred * smoothing_gain';
    
    // Ensure symmetry and positive definiteness
    cond_cov = ensure_symmetry(cond_cov);
    cond_cov = ensure_positive_definite(cond_cov, 1e-6);
    
    // Sample
    smoothed_samples[, t] = multi_normal_rng(cond_mean, cond_cov);
  }
  
  return smoothed_samples;
}

/**
 * Forward Filtering Backward Sampling using square root form for the SF model
 *
 * This function implements the complete FFBS algorithm with square root
 * formulation by first filtering then sampling.
 *
 * @param normalized_y Vector of normalized tumor measurements
 * @param time_points Vector of observation time points
 * @param decrease_rate Tumor regression rate constant (d)
 * @param growth_rate Tumor growth rate constant (g)
 * @param decrease_process_sd Process noise standard deviation for decrease component
 * @param growth_process_sd Process noise standard deviation for growth component
 * @param measure_sd Measurement noise standard deviation
 *
 * @return Matrix of sampled state vectors from the joint smoothing distribution
 */
matrix sf_forward_filter_backward_sample_sqrt_rng(
    vector normalized_y, 
    array[] real time_points, 
    real decrease_rate, 
    real growth_rate, 
    real decrease_process_sd, 
    real growth_process_sd, 
    real measure_sd
) {
  int N = num_elements(time_points);
  matrix[2, N] filtered_means;
  array[N] matrix[2, 2] filtered_sqrt_covs;
  
  // Get filtered estimates using square root filter
  (filtered_means, filtered_sqrt_covs) = sf_forward_filter_sqrt(
    normalized_y, time_points, decrease_rate, growth_rate, decrease_process_sd, growth_process_sd, measure_sd);
  
  // Perform backward sampling using filtered results
  return sf_backward_sample_sqrt_rng(
    filtered_means, filtered_sqrt_covs, time_points, decrease_rate, growth_rate, decrease_process_sd, growth_process_sd);
}

/**
 * Generate forecasted tumor measurements and states for the SF model
 *
 * This function first performs FFBS on observed data, then forecasts
 * future states and measurements beyond the observation period. It uses
 * the standard Kalman filter implementation.
 *
 * @param normalized_y Vector of normalized tumor measurements
 * @param time_points Vector of observation time points
 * @param forecast_time_points Vector of future time points for forecasting
 * @param decrease_rate Tumor regression rate constant (d)
 * @param growth_rate Tumor growth rate constant (g)
 * @param decrease_process_sd Process noise standard deviation for decrease component
 * @param growth_process_sd Process noise standard deviation for growth component
 * @param measure_sd Measurement noise standard deviation
 *
 * @return A tuple containing:
 *   - matrix[2, N+K]: states for observed and forecast periods
 *   - vector[N+K]: normalized y for observed and forecast periods
 */
tuple(matrix, vector) sf_forecast_rng(
  vector x_0,
  vector normalized_y,
  array[] int time_points,
  array[] real forecast_time_points,
  real decrease_rate,
  real growth_rate,
  real decrease_process_sd,
  real growth_process_sd,
  real measure_sd
) {
  if (min(forecast_time_points) <= max(time_points)) {
    fatal_error("All future time points must be after the last observed time point. Last observed = ", 
                max(time_points), " and first forecast = ", min(forecast_time_points));
  }

  int N = num_elements(time_points);
  int K = num_elements(forecast_time_points);
  int total_len = N + K;
  
  // Combined timepoints for reference
  array[total_len] real all_time_points;
  all_time_points[1:N] = time_points;
  all_time_points[(N+1):total_len] = forecast_time_points;
  
  // Output matrices
  matrix[2, total_len] states;
  vector[total_len] predictions;
  
  // First perform FFBS on the observed data using standard implementation
  matrix[2, N] filtered_means;
  array[N] matrix[2, 2] filtered_covs;
  
  // Forward filter observed data
  (filtered_means, filtered_covs) = sf_forward_filter(
      x_0, normalized_y, time_points, decrease_rate, growth_rate, decrease_process_sd, growth_process_sd, measure_sd);
  
  // Backward sample observed data
  matrix[2, N] smoothed_states = sf_backward_sample_rng(
      filtered_means, filtered_covs, time_points, decrease_rate, growth_rate, decrease_process_sd, growth_process_sd);
  
  // Store smoothed states for observed period
  states[, 1:N] = smoothed_states;
  predictions[1:N] = normalized_y;
  
  // Process noise covariance for forecasting
  matrix[2, 2] Q = diag_matrix([decrease_process_sd^2, growth_process_sd^2]');
  
  // Start forecasting from the last smoothed state
  vector[2] current_state = smoothed_states[, N];
  
  // Generate forecasts
  for (k in 1:K) {
      // Calculate time step size between forecast points
      real delta_t;
      if (k == 1) {
          delta_t = forecast_time_points[k] - time_points[N];
      } else {
          delta_t = forecast_time_points[k] - forecast_time_points[k-1];
      }
      
      // State transition matrix for this time step
      matrix[2, 2] F = sf_create_F(decrease_rate, growth_rate, delta_t);
      
      // Generate process noise
      vector[2] noise = multi_normal_rng([0, 0]', Q);
      
      // Forward predict state
      current_state = F * current_state + noise;
      
      // Store forecasted state
      states[, N+k] = current_state;
      
      // Generate measurement (normalized y value)
      // In the SF model, measurement = regression + growth
      real measurement_mean = current_state[1] + current_state[2];
      
      if (is_inf(measurement_mean)) {
        print("[k = ", k, "] current_state = ", current_state, ", F = ", F);
      }
      
      predictions[N+k] = normal_rng(measurement_mean, measure_sd);
  }
  
  return (states, predictions);
}

/**
 * Generate forecasted tumor measurements and states for the SF model
 *
 * This function first performs FFBS on observed data, then forecasts
 * future states and measurements beyond the observation period. It uses
 * the square root formulation for numerical stability.
 *
 * @param normalized_y Vector of normalized tumor measurements
 * @param time_points Vector of observation time points
 * @param forecast_time_points Vector of future time points for forecasting
 * @param decrease_rate Tumor regression rate constant (d)
 * @param growth_rate Tumor growth rate constant (g)
 * @param decrease_process_sd Process noise standard deviation for decrease component
 * @param growth_process_sd Process noise standard deviation for growth component
 * @param measure_sd Measurement noise standard deviation
 *
 * @return A tuple containing:
 *   - matrix[2, N+K]: states for observed and forecast periods
 *   - vector[N+K]: normalized y for observed and forecast periods
 */
tuple(matrix, vector) sf_forecast_sqrt_rng(
    vector normalized_y,
    array[] real time_points,
    array[] real forecast_time_points,
    real decrease_rate,
    real growth_rate,
    real decrease_process_sd,
    real growth_process_sd,
    real measure_sd
) {
  if (min(forecast_time_points) <= max(time_points)) {
    fatal_error("All future time points must be after the last observed time point. Last observed = ", 
                max(time_points), " and first forecast = ", min(forecast_time_points));
  }

  int N = num_elements(time_points);
  int K = num_elements(forecast_time_points);
  int total_len = N + K;
  
  // Combined timepoints for reference
  array[total_len] real all_time_points;
  all_time_points[1:N] = time_points;
  all_time_points[(N+1):total_len] = forecast_time_points;
  
  // Output matrices
  matrix[2, total_len] states;
  vector[total_len] predictions;
  
  // First perform FFBS on the observed data using square root implementation
  matrix[2, N] filtered_means;
  array[N] matrix[2, 2] filtered_sqrt_covs;
  
  // Forward filter observed data
  (filtered_means, filtered_sqrt_covs) = sf_forward_filter_sqrt(
      normalized_y, time_points, decrease_rate, growth_rate, decrease_process_sd, growth_process_sd, measure_sd);
  
  // Backward sample observed data
  matrix[2, N] smoothed_states = sf_backward_sample_sqrt_rng(
      filtered_means, filtered_sqrt_covs, time_points, decrease_rate, growth_rate, decrease_process_sd, growth_process_sd);
  
  // Store smoothed states for observed period
  states[, 1:N] = smoothed_states;
  predictions[1:N] = normalized_y;
  
  // Process noise for forecasting
  matrix[2, 2] sqrt_Q = diag_matrix([decrease_process_sd, growth_process_sd]');
  
  // Start forecasting from the last smoothed state
  vector[2] current_state = smoothed_states[, N];
  
  // Generate forecasts
  for (k in 1:K) {
      // Calculate time step size between forecast points
      real delta_t;
      if (k == 1) {
          delta_t = forecast_time_points[k] - time_points[N];
      } else {
          delta_t = forecast_time_points[k] - forecast_time_points[k-1];
      }
      
      // State transition matrix for this time step
      matrix[2, 2] F = sf_create_F(decrease_rate, growth_rate, delta_t);
      
      // Generate process noise
      vector[2] noise = multi_normal_rng([0, 0]', sqrt_Q * sqrt_Q');
      
      // Forward predict state
      current_state = F * current_state + noise;
      
      // Store forecasted state
      states[, N+k] = current_state;
      
      // Generate measurement (normalized y value)
      // In the SF model, measurement = regression + growth
      real measurement_mean = current_state[1] + current_state[2];
      predictions[N+k] = normal_rng(measurement_mean, measure_sd);
  }
  
  return (states, predictions);
}