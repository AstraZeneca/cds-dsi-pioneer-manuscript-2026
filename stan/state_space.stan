/**
 * Computes the discretized process noise covariance matrix for SDE representation
 *
 * This function uses the Van Loan method to compute the process noise covariance 
 * matrix for the discrete-time state-space model derived from a continuous-time SDE.
 * The formula is based on the solution to the matrix differential equation:
 * dQ(t)/dt = F*Q(t) + Q(t)*F' + L*q*L'
 *
 * @param F     The system dynamics matrix in the SDE dx/dt = Fx + Ldw
 * @param L     The noise input vector in the SDE dx/dt = Fx + Ldw
 * @param q     The spectral density of the continuous-time white noise
 * @param dt    The time step for discretization
 * @return      The discretized process noise covariance matrix with regularization
 */
matrix compute_Q(matrix F, vector L, real q, real dt) {
  int d = rows(F);
  
  // Form the expanded matrix for Van Loan method
  // [F, L*q*L'; 0, -F']
  matrix[2*d, 2*d] FQ;
  
  // Fill the blocks properly
  FQ[1:d, 1:d] = F;
  FQ[1:d, (d+1):(2*d)] = L * q * L';
  FQ[(d+1):(2*d), 1:d] = rep_matrix(0, d, d);
  FQ[(d+1):(2*d), (d+1):(2*d)] = -F';
  
  // Compute the matrix exponential
  matrix[2*d, 2*d] expFQ = matrix_exp(FQ * dt);
  
  // Extract the relevant blocks to form Q
  // Q = Φ₁₂ * Φ₂₂⁻¹, where Φ = exp(FQ*dt)
  matrix[d, d] Q = block(expFQ, 1, d+1, d, d) * 
                   block(expFQ, d+1, d+1, d, d)';
                   
  // Ensure numerical stability by enforcing symmetry
  matrix[d, d] Q_sym = ((Q + Q') / 2);
  
  // Add small regularization to guarantee positive-definiteness
  return Q_sym + diag_matrix(rep_vector(1e-8, d));
}

/**
 * Generates a trajectory from a Matern 1/2 GP using state-space form
 *
 * The Matern 1/2 kernel corresponds to an Ornstein-Uhlenbeck process
 * and produces sample paths that are continuous but not differentiable.
 * The SDE representation is:
 *   dx/dt = (-1/l)*x + dw
 * where l is the length scale and dw is a white noise process.
 *
 * @param ts            Array of time points at which to generate the process
 * @param n             Number of time points
 * @param alpha         Standard deviation parameter (determines process amplitude)
 * @param length_scale  Length scale parameter (determines process smoothness)
 * @return              Array of state vectors at each time point
 */
array[] vector matern12_ss_rng(array[] real ts, int n, real alpha, real length_scale) {
  int d = 1;  // State dimension for Matern 1/2
  array[n] vector[d] x;  // State at each time point
  
  // Initialize state at first time point with stationary covariance
  // P₀ = α² is the variance in the stationary distribution
  matrix[d, d] P0;
  P0[1,1] = square(alpha);
  
  // Initial state is drawn from a Gaussian with zero mean
  vector[d] m0 = zeros_vector(d);
  x[1] = multi_normal_rng(m0, P0);
  
  // State-space matrices from SDE formulation
  // F = [-1/l] is the system matrix
  matrix[d, d] F;
  F[1,1] = -1/length_scale;
  
  // L = [1] is the noise input vector
  vector[d] L;
  L[1] = 1;
  
  // q = 2α²/l is the spectral density of the white noise
  real q = (2 / length_scale) * square(alpha);
  
  // Forward simulation through time using discretized state-space model
  for (i in 2:n) {
    real dt = ts[i] - ts[i-1];
    
    // State transition matrix: A = exp(F*dt)
    matrix[d, d] A = matrix_exp(F * dt);
    
    // Process noise covariance (analytical for Matern 1/2)
    // Q = α²(1-exp(-2dt/l))
    matrix[d, d] Q;
    Q[1,1] = (1 - exp((-2 * dt) / length_scale)) * (square(alpha) / 2);
    
    // Add small regularization for numerical stability
    matrix[d, d] Q_reg = Q + diag_matrix(rep_vector(1e-8, d));
    
    // Generate next state according to x(t+dt) = Ax(t) + q(t)
    // where q(t) ~ N(0, Q)
    x[i] = multi_normal_rng(A * x[i-1], Q_reg);
  }
  
  return x;
}

/**
 * Generates a trajectory from a Matern 3/2 GP using state-space form
 *
 * The Matern 3/2 kernel produces sample paths that are once differentiable.
 * Its SDE representation is a 2-dimensional system:
 *   dx₁/dt = x₂
 *   dx₂/dt = -(√3/l)²x₁ - 2(√3/l)x₂ + dw
 * where l is the length scale and dw is a white noise process.
 *
 * @param ts            Array of time points at which to generate the process
 * @param n             Number of time points
 * @param alpha         Standard deviation parameter (determines process amplitude)
 * @param length_scale  Length scale parameter (determines process smoothness)
 * @return              Array of state vectors at each time point
 */
array[] vector matern32_ss_rng(array[] real ts, int n, real alpha, real length_scale) {
  int d = 2;  // State dimension for Matern 3/2
  array[n] vector[d] x;  // State at each time point
  
  // Initialize state at first time point with stationary covariance
  // The stationary covariance has specific structure for Matern 3/2
  matrix[d, d] P0;
  P0[1,1] = square(alpha);
  P0[1,2] = 0;
  P0[2,1] = 0;
  P0[2,2] = (3/(length_scale^2)) * square(alpha);
  
  // Initial state is drawn from a Gaussian with zero mean
  vector[d] m0 = zeros_vector(d);
  x[1] = multi_normal_rng(m0, P0);
  
  // State-space matrices from SDE formulation
  // F = [0 1; -3/l² -2√3/l] is the system matrix
  matrix[d, d] F;
  F[1,1] = 0;
  F[1,2] = 1;
  F[2,1] = -3/(length_scale^2);
  F[2,2] = -2*sqrt(3)/length_scale;
  
  // L = [0; 1] is the noise input vector
  vector[d] L;
  L[1] = 0;
  L[2] = 1;
  
  // q = 12√3α²/l³ is the spectral density of the white noise
  real q = ((12*sqrt(3))/(length_scale^3))*square(alpha);
  
  // Forward simulation through time using discretized state-space model
  for (i in 2:n) {
    real dt = ts[i] - ts[i-1];
    
    // State transition matrix: A = exp(F*dt)
    matrix[d, d] A = matrix_exp(F * dt);
    
    // Process noise covariance matrix using Van Loan method
    matrix[d, d] Q = compute_Q(F, L, q, dt);
    
    // Generate next state according to x(t+dt) = Ax(t) + q(t)
    // where q(t) ~ N(0, Q)
    x[i] = multi_normal_rng(A * x[i-1], Q);
  }
  
  return x;
}

/**
 * Generates a trajectory from a Matern 5/2 GP using state-space form
 *
 * The Matern 5/2 kernel produces sample paths that are twice differentiable.
 * Its SDE representation is a 3-dimensional system:
 *   dx₁/dt = x₂
 *   dx₂/dt = x₃
 *   dx₃/dt = -λ³x₁ - 3λ²x₂ - 3λx₃ + dw
 * where λ = √5/l is a parameter based on the length scale and dw is white noise.
 *
 * @param ts            Array of time points at which to generate the process
 * @param n             Number of time points
 * @param alpha         Standard deviation parameter (determines process amplitude)
 * @param length_scale  Length scale parameter (determines process smoothness)
 * @return              Array of state vectors at each time point
 */
array[] vector matern52_ss_rng(array[] real ts, int n, real alpha, real length_scale) {
  int d = 3;  // State dimension for Matern 5/2
  array[n] vector[d] x;  // State at each time point
  
  // Parameters derived from the kernel parameters
  real lambda = sqrt(5.0) / length_scale;  // λ = √5/l
  real kappa = (5.0 / 3.0) * (square(alpha) / square(length_scale));  // κ = (5/3)α²/l²
  
  // Initialize state covariance matrix for stationary distribution
  // This specific form ensures the process has Matern 5/2 covariance
  matrix[d, d] P0;
  P0[1,1] = square(alpha);
  P0[1,2] = 0;
  P0[1,3] = -kappa;
  P0[2,1] = 0;
  P0[2,2] = kappa;
  P0[2,3] = 0;
  P0[3,1] = -kappa;
  P0[3,2] = 0;
  P0[3,3] = (25.0 * square(alpha)) / pow(length_scale, 4);
  
  // Initial state is drawn from a Gaussian with zero mean
  vector[d] m0 = zeros_vector(d);
  x[1] = multi_normal_rng(m0, P0);
  
  // State-space matrices derived from the SDE formulation
  // F matrix encodes the differential equation system
  matrix[d, d] F;
  F[1,1] = 0;
  F[1,2] = 1;
  F[1,3] = 0;
  F[2,1] = 0;
  F[2,2] = 0;
  F[2,3] = 1;
  F[3,1] = -pow(lambda, 3);
  F[3,2] = -(3.0 * square(lambda));
  F[3,3] = -(3 * lambda);
  
  // L = [0; 0; 1] is the noise input vector
  vector[d] L;
  L[1] = 0;
  L[2] = 0;
  L[3] = 1;
  
  // q is the spectral density of the white noise process
  // q = (400√5α²)/(3l⁵)
  real q = ((400.0 * sqrt(5.0)) / 3.0) * (square(alpha) / pow(length_scale, 5));
  
  // Forward simulation through time using discretized state-space model
  for (i in 2:n) {
    real dt = ts[i] - ts[i-1];
    
    // State transition matrix: A = exp(F*dt)
    matrix[d, d] A = matrix_exp(F * dt);
    
    // Compute process noise covariance using Van Loan method
    matrix[d, d] Q = compute_Q(F, L, q, dt);
    
    // Generate next state according to x(t+dt) = Ax(t) + q(t)
    // where q(t) ~ N(0, Q)
    x[i] = multi_normal_rng(A * x[i-1], Q);
  }
  
  return x;
}

