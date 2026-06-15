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

// --- Baseline Hazard Hierarchy (decomposed, per-transition x per-level) ---
// The legacy single shared flag `enable_ms_level_baseline_hazard[n_levels]`
// (0-4) conflated three orthogonal decisions and forced every transition to
// share one mode. It is replaced by three orthogonal arrays indexed by
// (intercept slot, level). N_TRANS = 6 intercept slots, ordered:
//   1 = 01, 2 = 02, 3 = 03, 4 = 12_s (sojourn), 5 = 12_t (clock-fwd), 6 = 32.
// (Slot-index constants MS_SLOT_* are defined in transformed_data.stan.)
//
// 1. Hierarchical baseline GP residual per transition x level (was bundled in
//    legacy mode == 3 / RE_GP).
array[6, n_levels] int<lower=0, upper=1> enable_ms_level_gp;
//
// 2. Level intercept mode per transition x level (was bundled in legacy
//    mode 1/2/4). Codes: 0 = none, 1 = FE, 2 = RE (NCP), 3 = RE (CP).
array[6, n_levels] int<lower=0, upper=3> ms_level_intercept_mode;
//
// 3. Correlation grouping (cross-transition frailty). Same positive code at a
//    level => shared MVN/LKJ block; 0 = independent singleton (legacy
//    behavior). Phase 1 (decomposition) keeps this all-zero.
array[6, n_levels] int<lower=0> ms_level_intercept_corr_group;

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

// --- Time-Varying Covariate Basis Selector ---
// Selects the burden-coupling feature set for the time-varying covariates.
//   0 (default): legacy 3-feature basis [level, log_decrease_rate, log_growth_rate]
//                (n_time_varying_covar = 3). Features 2-3 are bi-exponential SSM
//                component parameters.
//   1          : (level, velocity) 2-feature basis [level, velocity]
//                (n_time_varying_covar = 2). Feature 2 = central-difference of the
//                latent log-burden trajectory. Required for the Laplace surrogate
//                (the bi-exponential rates have no analog under the quadratic
//                surrogate; velocity = d/dw of log-burden does). See
//                docs/superpowers/specs/2026-06-09-level-velocity-coupling-design.md
// The R config MUST set this together with a matching n_time_varying_covar:
//   (enable_ms_velocity_basis, n_time_varying_covar) in {(1,2), (0,3)}.
int<lower=0, upper=1> enable_ms_velocity_basis;

// --- 0->1 Baseline Log-Time Trend ---
// When enabled, the 0->1 population baseline log-hazard gets an additive
// monotone log-time trend term, a sibling to the GP residual in the baseline
// temporal unit: log_pop_lambda_01(t) = intercept + slope*g(t) + GP_resid(t),
// where g(t) = log(t) - ms_log_t_centering. OFF by default → size-0 parameter
// space, bit-identical to all existing fits. 0->1 transition only.
int<lower=0, upper=1> enable_ms_baseline_trend_01;  // 0→1 baseline log-time trend

// --- PSA-at-State-Entry Covariates ---
// When enabled, the last observed log-PSA before entering state 1 (1→2) or
// state 3 (3→2) is used as a time-invariant patient-level covariate, shifting
// the entire sojourn hazard up/down based on PSA burden at transition entry.
int<lower=0, upper=1> enable_ms_12_entry_covar;
int<lower=0, upper=1> enable_ms_32_entry_covar;
