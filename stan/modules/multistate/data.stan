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

// --- Observed 1→2 death calendar week ---
// Exact calendar week of death for progressed-then-died patients (0 otherwise).
// Used in GQ to avoid round-trip through pfs + pmax(1, death_week - pfs)
// which can overshoot when detection-adjustment pushes pfs past death_week.
array[n_patients] int<lower=0> ms_os_event_12;


// --- 0→3 Transition: Dropout ---
array[n_patients] int<lower=0> ms_time_03;  // Calendar week of dropout (= patient_max_t for all patients)

// --- 3→2 Transition: Off-trial death ---
array[n_patients] int<lower=0> ms_time_32;              // Sojourn time in state 3 until off-trial death, 0 if N/A
array[n_patients] int<lower=0, upper=1> ms_censored_32; // 1 = censored in state 3, 0 = off-trial death observed

// --- Deterministic Progression Flag ---
// Was progression determined by the mechanistic model (PSA-PD/RECIST-PD)?
// If 1, no hazard contribution at T₀₁ from the stochastic 0→1 component
array[n_patients] int<lower=0, upper=1> ms_prog_deterministic;

// --- Interval Censoring Gap (0→1 transition) ---
// Weeks between last clean assessment and progression detection.
// = 0 for: censored patients, RECIST-determined PD, death-without-PD.
// Used in transformed data to derive ms_ic_gap_01.
array[n_patients] int<lower=0> interval_censored;

// --- Time Grid Dimensions ---
// max_all_t is already defined in base data
// Max sojourn time grid (only needed if enable_ms_12=1)
int<lower=1> ms_max_sojourn_t;
// Max sojourn time grid for 3→2 (only needed if enable_ms_32=1)
int<lower=1> ms_max_sojourn_t_32;

// --- GP Knot Grid Resolution ---
// Number of weeks per GP knot (1 = weekly, 4 = 4-weekly, etc.)
// Coarser grids dramatically reduce Cholesky cost (n_knots^3).
// Recommended: 4 (48 knots from 191 weeks) for good performance.
int<lower=1> ms_gp_grid_step;

// --- Covariate Dimensions ---
int<lower=0> n_time_varying_covar;    // Number of time-varying covariates
int<lower=0> n_time_invariant_covar;  // Number of time-invariant covariates
