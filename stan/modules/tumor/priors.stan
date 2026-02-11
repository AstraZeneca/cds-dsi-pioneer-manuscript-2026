// ============================================================================
// TUMOR/SLD MEASUREMENT PRIORS
// ============================================================================
// Prior distributions for SLD observation model parameters

// SLD measurement error prior
// inv_gamma keeps mass away from zero, avoiding geometry issues
measure_sd_sld ~ inv_gamma(measure_sd_sld_alpha, measure_sd_sld_beta);
