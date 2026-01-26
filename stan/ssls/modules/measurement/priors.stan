// Measurement error model priors
// inv_gamma keeps mass away from zero, avoiding geometry issues

measure_sd ~ inv_gamma(measure_sd_alpha, measure_sd_beta);
