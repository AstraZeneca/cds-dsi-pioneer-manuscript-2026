// ============================================================================
// Other Events Model Flags
// ============================================================================

// Baseline hazard hierarchy (kept as trial-level for Phase 6A; GP doesn't generalize easily to N-level)
int<lower=0,upper=1> oe_enable_trial_baseline_hazard;

// Population-level covariate effects (slopes only, no intercepts)
int<lower=0,upper=1> oe_enable_pop_cov;
int<lower=0,upper=1> oe_enable_pop_tumor_cov;  // Time-varying log(SLD) from states

// Multi-level random slopes for non-tumor covariates (replaces oe_enable_trial_cov)
// Index 1 = first grouping level (e.g., trial), Index n_levels = patient level
array[n_levels] int<lower=0,upper=1> oe_enable_level_cov;

// Trial-level random slopes for tumor covariates (scaffolded for future use)
int<lower=0,upper=1> oe_enable_trial_tumor_cov;
