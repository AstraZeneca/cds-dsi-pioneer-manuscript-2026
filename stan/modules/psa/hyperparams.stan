// ============================================================================
// PSA MEASUREMENT HYPERPARAMETERS
// ============================================================================
// Prior hyperparameters for PSA observation model

// PSA measurement error prior hyperparameters
real<lower=0> measure_sd_psa_alpha;  // inv_gamma shape for PSA measurement noise SD
real<lower=0> measure_sd_psa_beta;   // inv_gamma scale for PSA measurement noise SD
