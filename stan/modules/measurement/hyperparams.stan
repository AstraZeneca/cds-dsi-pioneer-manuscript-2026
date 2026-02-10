// Measurement error model hyperparameters
// These control the observation model noise (independent measurement error)

real<lower = 0> measure_sd_alpha;  // inv_gamma shape for measurement noise SD
real<lower = 0> measure_sd_beta;   // inv_gamma scale for measurement noise SD
