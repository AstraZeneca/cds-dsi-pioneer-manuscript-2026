// ============================================================================
// Other Events Model Transformed Data
// ============================================================================

// --- Independent Events Framework ---
// INDEPENDENT EVENTS MODEL:
// Target progression (from tumor dynamics) and other-events (from this model)
// are modeled as INDEPENDENT processes. Both can be observed for the same patient.
//
// Other events include: non-target progression, new lesions, death, dropout, etc.
// Some other events (e.g., non-target PD) may be non-terminal - patient continues.

array[n_patients] int<lower=0> ic_other_events_pfs; 
array[n_patients] int<lower=0> other_events_interval_censored;

for (i in 1:n_patients) {
  // Interval censoring only applies when other events observed (not censored)
  // Note: We assume interval_censored applies to whichever event occurred
  // If both events observed, we use the same interval uncertainty for other-events
  other_events_interval_censored[i] = other_events_right_censored[i] ? 0 : interval_censored[i];
  ic_other_events_pfs[i] = other_events_pfs[i] + other_events_interval_censored[i]; 
}

// --- QR Decomposition for Covariates ---
matrix[n_patients, n_tumor_covar] Q_tumor_sum_covar;
matrix[n_tumor_covar, n_tumor_covar] R_tumor_sum_covar;

if (n_tumor_covar > 0) {
  Q_tumor_sum_covar = qr_thin_Q(tumor_sum_covar) * sqrt(n_patients - 1);
  R_tumor_sum_covar = qr_thin_R(tumor_sum_covar) / sqrt(n_patients - 1);
} else {
  Q_tumor_sum_covar = rep_matrix(0, n_patients, 0);
  R_tumor_sum_covar = rep_matrix(0, 0, 0);
}

// Q_covar_design_matrix and R_covar_design_matrix already in _sf_transformed_data.stan
