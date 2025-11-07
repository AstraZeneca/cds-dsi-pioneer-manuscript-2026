// ============================================================================
// Other Events Model Transformed Data
// ============================================================================

// --- Competing Risks Censoring Logic ---
// COMPETING RISKS FRAMEWORK:
// Overall PFS is the minimum of:
//   1. Target lesion PFS (from tumor dynamics model)
//   2. Other events PFS (from this model)
//
// Other events include: non-target progression, new lesions, death, etc.

array[n_patients] int<lower=0> other_events_pfs, ic_other_events_pfs; 
array[n_patients] int<lower=0> other_events_interval_censored;
array[n_patients] int<lower=0,upper=1> other_events_right_censored = zeros_int_array(n_patients);

for (i in 1:n_patients) {
  // Other events PFS equals overall PFS time (minimum of all causes)
  other_events_pfs[i] = pfs[i];
  
  // Other events are censored if:
  // 1. Patient was overall censored (no event observed), OR
  // 2. Target lesions progressed first (competing event)
  other_events_right_censored[i] = right_censored[i] || (!target_right_censored[i] && pfs[i] >= target_pfs[i]);
  
  // Interval censoring only applies when other events observed (not censored)
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
