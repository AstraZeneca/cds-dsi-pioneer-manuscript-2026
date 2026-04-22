// ============================================================================
// PSA MEASUREMENT DATA
// ============================================================================
// PSA-specific longitudinal measurements and outcomes.
// This module contains data aligned with the visit schedule defined in _base_data.stan.
//
// Array lengths: All measurement arrays have length sum(n_patient_visits) and are
// aligned with t_patient_visits from _base_data.stan.

// PSA values at each visit (ng/mL)
vector<lower=0>[sum(n_patient_visits)] psa_values;

// Indicator for whether PSA was measured at each visit (handles missing data)
array[sum(n_patient_visits)] int<lower=0, upper=1> psa_measured;

// ============================================================================
// PSA-BASED PROGRESSION OUTCOMES (PCWG3)
// ============================================================================

// PSA progression-free survival
array[n_patients] int<lower=0> psa_pfs;  // Weeks after baseline
array[n_patients] int<lower=0, upper=1> psa_right_censored;
array[n_patients] int<lower=0> psa_interval_censored;  // Weeks after psa_pfs that actual progression could have occurred

// PCWG3 response category at each visit
// 1 = Undetectable (CR equivalent)
// 2 = PSA50 (PR equivalent)
// 3 = Stable (SD equivalent)
// 4 = PSA-PD (confirmed progression)
// 5 = Not Evaluable
array[sum(n_patient_visits)] int<lower=1, upper=5> pcwg3_category;

// ============================================================================
// PSA-SPECIFIC THRESHOLDS
// ============================================================================

// PSA undetectable threshold (protocol: < 0.2 ng/mL = CR)
real<lower=0> psa_undetectable_threshold;
