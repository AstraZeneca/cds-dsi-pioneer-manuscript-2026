// ============================================================================
// Multistate Hazard Model Flags
// ============================================================================
// Configurable transitions for flexible multistate modeling:
//   - SLD model: (1, 0, 0) - single 0→1 transition (PFS event)
//   - PSA model: (1, 1, 1) - full illness-death

// --- Transition Enable Flags ---
int<lower=0, upper=1> enable_ms_01;  // 0→1: Progression / PFS event
int<lower=0, upper=1> enable_ms_02;  // 0→2: Death without progression
int<lower=0, upper=1> enable_ms_12;  // 1→2: Post-progression death

// --- Time Scale for 1→2 Transition ---
// 0 = Markov (clock-forward t)
// 1 = Semi-Markov (sojourn time s = t - T₀₁)
// 2 = Extended (additive: GP_t(t) + GP_s(s))
int<lower=0, upper=2> ms_time_scale_12;

// Shared GP shape for 0→2 and 1→2 clock-forward death transitions.
// When enabled, alpha/rho/eta are shared; intercepts remain separate.
// Only valid when ms_time_scale_12 == 0 (Markov). Sojourn GP unaffected.
int<lower=0, upper=1> share_dead_gp_shape;

// --- Baseline Hazard Hierarchy ---
// Uses generic n_levels structure (consistent with tr, frac, init modules)
// 5-state flag: 0 = no level effect, 1 = FE intercept, 2 = RE intercept (NCP),
// 3 = RE + GP residual, 4 = RE intercept (CP)
array[n_levels] int<lower=0, upper=4> enable_ms_level_baseline_hazard;

// --- Dropout / Off-trial Transitions ---
// 0→3: Dropout (off-trial). Constant hazard (no GP).
int<lower=0, upper=1> enable_ms_03;

// 3→2: Off-trial death. Markovian GP baseline hazard (clock-forward).
int<lower=0, upper=1> enable_ms_32;

// --- Covariate Flags ---
int<lower=0, upper=1> enable_ms_pop_time_invariant_cov;  // Population-level time-invariant covariates
int<lower=0, upper=1> enable_ms_pop_time_varying_cov;    // Population-level time-varying covariates
array[n_levels] int<lower=0, upper=1> enable_ms_level_cov;  // Per-level random slopes for time-invariant

// --- Visit-Gated 0->1 Mode ---
// When enabled, 0->1 hazard accumulates only at assessment visit weeks
// (not every calendar week). Interval censoring for 0->1 is disabled.
// By default, time-varying covariates for 0->1 are built from observed PSA
// (data). Set enable_ms_visit_gated_latent_01=1 to use the modeled latent
// PSA trajectory instead (still evaluated only at visit weeks).
int<lower=0, upper=1> enable_ms_visit_gated_01;
int<lower=0, upper=1> enable_ms_visit_gated_latent_01;

// --- 0->2 Time-Varying Covariate ---
// When enabled, 0->2 uses modeled PSA trajectory as time-varying covariate.
// When disabled, 0->2 uses GP baseline + time-invariant covariates only.
int<lower=0, upper=1> enable_ms_02_time_varying_cov;

// --- 0->3 / 3->2 Population-Level Covariate Flags ---
// Per-transition switches that gate population covariate effects on top of the
// existing module-wide enable_ms_pop_time_invariant_cov / enable_ms_pop_time_varying_cov.
// 0->3 supports both time-invariant and time-varying (tumor bridge) covariates.
// 3->2 supports time-invariant only (sojourn-clock indexing for time-varying not implemented).
int<lower=0, upper=1> enable_ms_03_time_invariant_cov;
int<lower=0, upper=1> enable_ms_32_time_invariant_cov;
int<lower=0, upper=1> enable_ms_03_time_varying_cov;

// --- PSA-at-State-Entry Covariates ---
// When enabled, the last observed log-PSA before entering state 1 (1→2) or
// state 3 (3→2) is used as a time-invariant patient-level covariate, shifting
// the entire sojourn hazard up/down based on PSA burden at transition entry.
int<lower=0, upper=1> enable_ms_12_entry_covar;
int<lower=0, upper=1> enable_ms_32_entry_covar;
