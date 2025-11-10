// ============================================================================
// Other Events Model Flags
// ============================================================================

// Baseline hazard hierarchy
int<lower=0,upper=1> oe_enable_trial_baseline_hazard;

// Population-level covariate effects (slopes only, no intercepts)
int<lower=0,upper=1> oe_enable_pop_cov;
int<lower=0,upper=1> oe_enable_pop_tumor_cov;  // Time-varying log(SLD) from states

// Trial-level random slopes (deviations from population)
int<lower=0,upper=1> oe_enable_trial_cov;
int<lower=0,upper=1> oe_enable_trial_tumor_cov;  // Scaffolded for future use
