// Phase 1 (revised): Verify nightly stanc3 supports laplace_marginal_tol
// Using the correct API from stanc3 test/integration/good/code-gen/laplace_functionals.stan

functions {
  // Log-likelihood functor: N patients, 1 latent NCP param each
  // y_i ~ Normal(mu + sigma * theta_i, obs_sd)
  real ll_fn(vector theta,    // latent NCP params (length N)
             real mu, real sigma, real obs_sd,
             vector y) {
    return normal_lpdf(y | mu + sigma * theta, obs_sd);
  }

  // Prior covariance functor: identity matrix (standard normal NCP prior)
  // Takes (N, dummy) since K_fn args must be a tuple (≥2 elements)
  matrix K_fn(int N, int dummy) {
    return identity_matrix(N);
  }
}

data {
  int<lower=1> N;
  vector[N] y;
  real<lower=0> obs_sd;
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
}

model {
  mu ~ normal(0, 10);
  sigma ~ normal(0, 5);

  // Laplace marginal — integrates out theta analytically
  target += laplace_marginal_tol(
    ll_fn,
    (mu, sigma, obs_sd, y),
    hessian_block_size,
    K_fn,
    (N, N),
    (theta_0, tolerance, max_num_steps, solver, max_steps_line_search, allow_fallthrough)
  );
}
