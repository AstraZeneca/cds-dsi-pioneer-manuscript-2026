// Toy model to verify reduce_sum threading works.
//
// y_i ~ Normal(mu, sigma), priors on mu and sigma.
// Likelihood accumulated via reduce_sum over observations.

functions {
  real partial_sum(
    array[] real y_slice,
    int start, int end,
    real mu,
    real sigma
  ) {
    return normal_lpdf(y_slice | mu, sigma);
  }
}

data {
  int<lower=1> N;
  array[N] real y;
  int<lower=1> grainsize;
}

parameters {
  real mu;
  real<lower=0> sigma;
}

model {
  mu ~ normal(0, 10);
  sigma ~ exponential(1);
  target += reduce_sum(partial_sum, y, grainsize, mu, sigma);
}
