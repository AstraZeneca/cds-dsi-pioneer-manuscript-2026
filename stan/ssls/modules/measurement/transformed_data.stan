// Measurement error model transformed data
// Limit of detection (LOD) for tumor size measurements

real<lower = 0> lod = 0.1; // cm - limit of detection for tumor size
real log_lod = log(0.1);   // Log of LOD for use in observation model
