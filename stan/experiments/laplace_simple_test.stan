// Phase 2: Simple 1D Laplace experiment (built-in laplace_marginal_tol)
// ============================================================================
// Model: y_i ~ Normal(mu + sigma * z_i, obs_sd), z_i ~ Normal(0, 1)
//
// Three fitting modes controlled by `laplace_mode`:
//   0 = Full HMC (z_i are sampled parameters)
//   1 = Built-in Laplace marginalization (z_i integrated out)
//   2 = Exact marginalization (uses known y_i ~ N(mu, sqrt(obs_sd^2+sigma^2)))
//
// For this Gaussian-Gaussian model, modes 1 and 2 should give identical
// results, and mode 0 should match within MCSE.
// ============================================================================

functions {
  // Log-likelihood functor: y_i ~ Normal(mu + sigma * theta_i, obs_sd)
  real ll_fn(vector theta,
             real mu, real sigma, real obs_sd,
             vector y) {
    return normal_lpdf(y | mu + sigma * theta, obs_sd);
  }

  // Prior covariance functor: identity matrix (standard normal NCP prior)
  matrix K_fn(int N, int dummy) {
    return identity_matrix(N);
  }
}

data {
  int<lower=1> N;
  vector[N] y;
  real<lower=0> obs_sd;
  int<lower=0, upper=2> laplace_mode;
}

transformed data {
  vector[N] theta_0 = rep_vector(0.0, N);
  real tolerance = 1e-6;
  int max_num_steps = 100;
  int hessian_block_size = 1;
  int solver = 1;
  int max_steps_line_search = 0;
  int allow_fallthrough = 1;
}

parameters {
  real mu;
  real<lower=0> sigma;
  vector[laplace_mode == 0 ? N : 0] z;
}

model {
  mu ~ normal(0, 10);
  sigma ~ normal(0, 5);

  if (laplace_mode == 0) {
    // MODE 0: Full HMC
    z ~ std_normal();
    for (i in 1:N) {
      y[i] ~ normal(mu + sigma * z[i], obs_sd);
    }
  } else if (laplace_mode == 1) {
    // MODE 1: Built-in Laplace marginalization
    target += laplace_marginal_tol(
      ll_fn,
      (mu, sigma, obs_sd, y),
      hessian_block_size,
      K_fn,
      (N, N),
      (theta_0, tolerance, max_num_steps, solver, max_steps_line_search, allow_fallthrough)
    );
  } else {
    // MODE 2: Exact marginalization
    real marginal_sd = sqrt(square(obs_sd) + square(sigma));
    y ~ normal(mu, marginal_sd);
  }
}
