// ============================================================================
// Other Events Model Data
// ============================================================================
// Data dimensions and arrays specific to the other events (competing risks) model

// Number of tumor-derived covariates for proportional hazards
// Typically: baseline SLD, sum of tumor sizes, etc.
int<lower=0> n_tumor_covar;

// Tumor-derived covariate matrix [n_patients x n_tumor_covar]
// Note: This will be QR-decomposed in transformed_data if n_tumor_covar > 0
matrix[n_patients, n_tumor_covar] tumor_sum_covar;
