int<lower = 0> n_bootstrap_param;
array[n_bootstrap_param] int<lower = 1> prediction_week; // At what week are starting our prediction
array[n_bootstrap_param] real<lower = 0> recruit_lambda; // neg binom rate
real<lower = 0> recruit_phi; // neg binom dispersion
