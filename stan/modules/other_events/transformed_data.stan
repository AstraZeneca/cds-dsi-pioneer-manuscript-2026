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

// Pre-computed flat indices for slopes
array[n_patients, n_levels] int patient_oe_slope_flat_idx;
{
  for (i in 1:n_patients) {
    for (lv in 1:n_levels) {
      if (oe_enable_level_cov[lv]) {
        patient_oe_slope_flat_idx[i, lv] =
          get_global_group_idx(enabled_level_pos_oe_slope, lv, patient_level_groups[i, lv]);
      } else {
        patient_oe_slope_flat_idx[i, lv] = 1;
      }
    }
  }
}

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
// BIOMARKER NORMALIZATION CONSTANTS FOR PROPORTIONAL HAZARDS
// ============================================================================
// These variables use generic "biomarker" naming to support both SLD and PSA.
// The actual values come from whichever observation module is active.
//
// IMPORTANT: The main model file must define `log_baseline_biomarker` before
// including this module. For SLD: log_baseline_biomarker = log_baseline_sld
// For PSA: log_baseline_biomarker = log_baseline_psa

// Compute robust normalization statistics (median and IQR) from observed log(biomarker)
// Using median/IQR instead of mean/SD for robustness to outliers and extreme values
// This avoids issues with zero or very small biomarker values
real median_log_biomarker_obs;
real iqr_log_biomarker_obs;

{
  // Get all observed log(biomarker) values across all patients and visits
  // Filter out zeros to avoid -Inf
  // Currently uses sum_tumor_size; future: could switch based on observation_type
  vector[sum(n_patient_visits)] log_biomarker_all_obs;
  int n_positive = 0;

  for (i in 1:sum(n_patient_visits)) {
    if (sum_tumor_size[i] > 0) {
      n_positive += 1;
      log_biomarker_all_obs[n_positive] = log(sum_tumor_size[i]);
    }
  }

  // Compute robust normalization using median and IQR from positive observations
  {
    array[3] real quantiles_obs = quantile(log_biomarker_all_obs[1:n_positive], {0.25, 0.5, 0.75});
    median_log_biomarker_obs = quantiles_obs[2];
    iqr_log_biomarker_obs = quantiles_obs[3] - quantiles_obs[1];
  }
}

