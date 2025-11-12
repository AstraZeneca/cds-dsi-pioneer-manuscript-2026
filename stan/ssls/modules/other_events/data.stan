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
// Typically: baseline SLD, sum of tumor sizes, etc.
int<lower=0> n_tumor_covar;

// Tumor-derived covariate matrix [n_patients x n_tumor_covar]
// Note: This will be QR-decomposed in transformed_data if n_tumor_covar > 0
matrix[n_patients, n_tumor_covar] tumor_sum_covar;
