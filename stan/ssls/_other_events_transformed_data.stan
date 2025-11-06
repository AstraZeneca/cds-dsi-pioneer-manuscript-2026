// Recreated from tumor/other_events_transformed_data.stan for modular ssls build.
// 
// COMPETING RISKS FRAMEWORK:
// This code implements a competing risks model where overall PFS is the minimum of:
//   1. Target lesion PFS: progression based on target lesion SLD (from tumor dynamics model)
//   2. Other events PFS: progression from ALL other causes
// 
// "Other events" include:
//   - Non-target lesion progression
//   - New lesion appearance  
//   - Death
//   - Any other progression cause not captured by target lesion SLD
//
// The other events hazard model (GP baseline hazard) captures everything except
// target lesion progression, creating a two-component competing risks model.

array[n_patients] int<lower = 0> other_events_pfs, ic_other_events_pfs; 
array[n_patients] int<lower = 0> other_events_interval_censored;
array[n_patients] int<lower = 0, upper = 1> other_events_right_censored = zeros_int_array(n_patients);

for (i in 1:n_patients) {
  // Other events PFS time equals overall PFS time (last visit without progression)
  // This is always equal to pfs[i] because pfs = min(target_pfs, true_other_events_pfs)
  other_events_pfs[i] = pfs[i];
  
  // Other events are censored if:
  // 1. Patient was censored overall (no event observed), OR
  // 2. Target lesions progressed first (competing event)
  other_events_right_censored[i] = right_censored[i] || (!target_right_censored[i] && pfs[i] >= target_pfs[i]);
  
  // Interval censoring only applies when other events are observed (not censored)
  // When censored, we observe survival up to pfs[i] with no interval uncertainty
  other_events_interval_censored[i] = other_events_right_censored[i] ? 0 : interval_censored[i];

  ic_other_events_pfs[i] = other_events_pfs[i] + other_events_interval_censored[i]; 
}
