// ============================================================================
// Other Events Model Data
// ============================================================================
// Data dimensions and arrays specific to the other events model
//
// NOTE: This model treats "other events" as INDEPENDENT from target progression.
// Other events include: non-target progression, new lesions, death, dropout, etc.
// Both target progression AND other-events can be observed for the same patient.

// Other events progression-free survival (weeks from baseline)
// This is independent from target_pfs - both can be observed
array[n_patients] int<lower=0> other_events_pfs;

// Whether other events were right-censored (1=censored, 0=event observed)
array[n_patients] int<lower=0, upper=1> other_events_right_censored;

// Number of tumor-derived covariates for proportional hazards
// Tumor covariates are extracted from the state-space model:
//   1. log(SLD) - current tumor burden (time-varying)
//   2. log(decrease rate) - patient-specific regression rate (time-invariant, from tr + frac modules)
//   3. log(growth rate) - patient-specific growth rate (time-invariant, from tr + frac modules)
// Note: These are NOT passed as data but extracted from states_full_grid in transformed_parameters
// n_tumor_covar is kept for backward compatibility but currently unused (scaffolding)
int<lower=0> n_tumor_covar;

