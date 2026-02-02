// ============================================================================
// Other Events Model Transformed Data
// ============================================================================

// --- Module-Specific Enabled Group Counts ---
// Compute enabled group count for oe slopes (no intercepts in this module)
int n_enabled_groups_oe_slope = compute_n_enabled_groups(
  n_groups_per_level, oe_enable_level_cov
);

// Position array for enabled slope levels only
array[n_levels + 1] int enabled_level_pos_oe_slope = create_enabled_pos(
  n_groups_per_level, oe_enable_level_cov
);

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

// ============================================================================
// SLD normalization constants for tumor covariates in proportional hazards
// ============================================================================

// Compute robust normalization statistics (median and IQR) from observed log(SLD)
// Using median/IQR instead of mean/SD for robustness to outliers and extreme values
// This avoids issues with zero or very small SLD values
real median_log_sld_obs;
real iqr_log_sld_obs;

// Patient-level baseline SLD (in cm) for converting normalized states to absolute SLD
vector[n_patients] log_baseline_sld;

{
  // Get all observed log(SLD) values across all patients and visits
  // Filter out zeros to avoid -Inf
  vector[sum(n_patient_visits)] log_sld_all_obs;
  int n_positive = 0;
  
  for (i in 1:sum(n_patient_visits)) {
    if (sum_tumor_size[i] > 0) {
      n_positive += 1;
      log_sld_all_obs[n_positive] = log(sum_tumor_size[i]);
    }
  }
  
  // Compute robust normalization using median and IQR from positive observations
  array[3] real quantiles_obs = quantile(log_sld_all_obs[1:n_positive], {0.25, 0.5, 0.75});
  median_log_sld_obs = quantiles_obs[2];
  iqr_log_sld_obs = quantiles_obs[3] - quantiles_obs[1];
  
  // Also store baseline SLD for each patient (for converting states)
  for (i in 1:n_patients) {
    int visit_start, visit_end;
    (visit_start, visit_end) = get_pos(patient_visit_pos, i);
    log_baseline_sld[i] = log(sum_tumor_size[visit_start]);
  }
}

