// ============================================================================
// TUMOR/SLD MEASUREMENT HYPERPARAMETERS
// ============================================================================
// Prior hyperparameters for SLD observation model

// SLD measurement error prior hyperparameters
real<lower=0> measure_sd_sld_alpha;  // inv_gamma shape for measurement noise SD
real<lower=0> measure_sd_sld_beta;   // inv_gamma scale for measurement noise SD
real<lower=2> measure_nu_sld;        // Student-t degrees of freedom (fixed); use ~5 for robust, higher → Gaussian
