// ============================================================================
// Multistate Hazard Model Data
// ============================================================================
// Patient-level data for multistate transitions

// --- Patient Final State ---
// Upper bound depends on reachable states:
//   - State 2 reachable if enable_ms_02=1 (0→2) or enable_ms_12=1 (1→2)
//   - SLD model (1,0,0): only states 0,1 → upper=1
//   - PSA model with death: states 0,1,2 → upper=2
// Note: ms_max_state computed in transformed_data
array[n_patients] int<lower=0> ms_final_state;
// 0 = no event (right-censored)
// 1 = progressed / PFS event (state depends on model)
// 2 = dead (only if enable_ms_02 or enable_ms_12)

// --- Transition Times ---
// Time to each transition (0 if not observed/applicable)
array[n_patients] int<lower=0> ms_time_01;  // Time to progression (0→1)
array[n_patients] int<lower=0> ms_time_02;  // Time to death w/o progression (0→2)
array[n_patients] int<lower=0> ms_time_12;  // Post-progression survival / sojourn time (1→2)

// --- Censoring Indicators ---
// 1 = censored for this transition, 0 = event observed
array[n_patients] int<lower=0, upper=1> ms_censored_01;
array[n_patients] int<lower=0, upper=1> ms_censored_02;
array[n_patients] int<lower=0, upper=1> ms_censored_12;

// --- Deterministic Progression Flag ---
// Was progression determined by the mechanistic model (PSA-PD/RECIST-PD)?
// If 1, no hazard contribution at T₀₁ from the stochastic 0→1 component
array[n_patients] int<lower=0, upper=1> ms_prog_deterministic;

// --- Time Grid Dimensions ---
// max_all_t is already defined in base data
// Max sojourn time grid (only needed if enable_ms_12=1)
int<lower=1> ms_max_sojourn_t;

// --- Covariate Dimensions ---
int<lower=0> n_time_varying_covar;    // Number of time-varying covariates
int<lower=0> n_time_invariant_covar;  // Number of time-invariant covariates
