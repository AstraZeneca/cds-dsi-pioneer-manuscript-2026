// ============================================================================
// PSA MEASUREMENT PRIORS
// ============================================================================
// Prior distributions for PSA observation model parameters

// PSA measurement error prior
// inv_gamma keeps mass away from zero, avoiding geometry issues
measure_sd_psa ~ inv_gamma(measure_sd_psa_alpha, measure_sd_psa_beta);
