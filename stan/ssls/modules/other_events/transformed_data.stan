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
