# Multistate Hazard Model Implementation Plan

## Context

The PIONEER model currently uses an `other_events` module for non-target progression and death as a single competing risk. We need a **flexible multistate module** that can support:

- **Current SLD model**: Single transition (0→1) - progression/death as one PFS event, equivalent to current `other_events`
- **Future PSA model**: Full illness-death (0→1, 0→2, 1→2) - separates progression from death

The design specification is in `docs/plan/PSA_260209.qmd` Section 6 (@sec-multistate).

**Key Design Decisions**:

1. **Configurable transitions** via flags - not hardcoded to 3 transitions
2. PSA-PD/RECIST-PD is **deterministic** from SF model; λ₀₁ only for "other" progression
3. Unified covariate structure with transition-specific coefficients
4. Support multiple time scales for 1→2 (Markov, semi-Markov, extended with additive GPs)

---

## Module Structure

**Decision: Create new `stan/modules/multistate/` module that replaces `other_events`**

Rationale: The multistate model is a generalization. Current `other_events` is a special case (single 0→2 transition). After implementation, `other_events/` will be removed.

```
stan/modules/multistate/
├── flags.stan
├── data.stan
├── hyperparams.stan
├── transformed_data.stan
├── parameters.stan
├── transformed_parameters.stan
└── priors.stan

stan/multistate.stanfunctions  # Likelihood function
```

---

## Data Structure

### Flags (`stan/modules/multistate/flags.stan`)

```stan
// Which transitions are enabled
int<lower=0, upper=1> enable_ms_01;  // 0→1: Progression
int<lower=0, upper=1> enable_ms_02;  // 0→2: Death without progression
int<lower=0, upper=1> enable_ms_12;  // 1→2: Post-progression death

// Time scale for 1→2: 0=Markov (clock-forward t), 1=semi-Markov (sojourn s), 2=extended (additive)
int<lower=0, upper=2> ms_time_scale_12;

// Hierarchy flags - all use generic n_levels structure (consistent with tr, frac, init)
array[n_levels] int<lower=0, upper=1> enable_ms_level_baseline_hazard;  // Per-level GP for baseline hazard
int<lower=0, upper=1> enable_ms_pop_time_invariant_cov;                 // Population-level time-invariant covariates
int<lower=0, upper=1> enable_ms_pop_time_varying_cov;                   // Population-level time-varying covariates
array[n_levels] int<lower=0, upper=1> enable_ms_level_cov;              // Per-level random slopes for time-invariant
```

**Configuration Examples:**

| Model | enable_ms_01 | enable_ms_02 | enable_ms_12 | Description |
|-------|--------------|--------------|--------------|-------------|
| SLD (current) | 1 | 0 | 0 | Single 0→1 transition (PFS event), equivalent to `other_events` |
| PSA full | 1 | 1 | 1 | Full illness-death model |
| PSA no post-prog | 1 | 1 | 0 | Progression + death w/o progression, no 1→2 |

### Stan Data (`stan/modules/multistate/data.stan`)

```stan
// Patient state at end of observation
// Upper bound depends on reachable states:
//   - State 2 reachable if enable_ms_02=1 (0→2) or enable_ms_12=1 (1→2)
//   - SLD model (1,0,0): only states 0,1 → upper=1
//   - PSA model with death: states 0,1,2 → upper=2
int ms_max_state = (enable_ms_02 || enable_ms_12) ? 2 : 1;
array[n_patients] int<lower=0, upper=ms_max_state> ms_final_state;
// 0=no event, 1=progressed/PFS event, 2=dead (if reachable)

// Transition times (0 if not observed/applicable)
array[n_patients] int<lower=0> ms_time_01;  // Time to progression
array[n_patients] int<lower=0> ms_time_02;  // Time to death w/o progression
array[n_patients] int<lower=0> ms_time_12;  // Post-progression survival (sojourn)

// Censoring per transition
array[n_patients] int<lower=0, upper=1> ms_censored_01;
array[n_patients] int<lower=0, upper=1> ms_censored_02;
array[n_patients] int<lower=0, upper=1> ms_censored_12;

// Was progression deterministic (from PSA/RECIST)?
array[n_patients] int<lower=0, upper=1> ms_prog_deterministic;

// Max sojourn time grid (only needed if enable_ms_12=1)
int<lower=1> max_sojourn_t;
```

### Patient Scenario Encoding

| Scenario | final_state | time_01 | time_02 | time_12 | Requires |
|----------|-------------|---------|---------|---------|----------|
| Right-censored | 0 | 0 | 0 | 0 | — |
| PFS event (SLD) / Progressed (PSA) | 1 | T₀₁ | 0 | 0 | enable_ms_01 |
| Died w/o progression | 2 | 0 | T₀₂ | 0 | enable_ms_02 |
| Progressed, then died | 2 | T₀₁ | 0 | s | enable_ms_01 ∧ enable_ms_12 |

### Transformed Data (`stan/modules/multistate/transformed_data.stan`)

```stan
// Derived flags for conditional parameter sizing

// Time scale flags for 1→2
int need_12_s_gp = enable_ms_12 && (ms_time_scale_12 == 1 || ms_time_scale_12 == 2);
int need_12_t_gp = enable_ms_12 && (ms_time_scale_12 == 0 || ms_time_scale_12 == 2);

// Baseline hazard N-level hierarchy (same indexing pattern as tr, frac, init)
int n_enabled_groups_ms_baseline_01 = ...;  // Computed from enable_ms_level_baseline_hazard & enable_ms_01
int n_enabled_groups_ms_baseline_02 = ...;  // Computed from enable_ms_level_baseline_hazard & enable_ms_02
int n_enabled_groups_ms_baseline_12_s = ...; // For sojourn GP
int n_enabled_groups_ms_baseline_12_t = ...; // For clock-forward GP

// Covariate N-level hierarchy (same indexing pattern as tr, frac, init)
int n_enabled_groups_ms_slope = ...;  // Computed from enable_ms_level_cov
```

---

## Parameters

Parameters are declared **conditionally** based on enabled transitions. Stan requires fixed array sizes, so we use size 0 for disabled transitions.

Hierarchy uses generic N-level structure throughout (consistent with tr, frac, init modules):
- **Baseline hazard GP**: Population + N-level (via `enable_ms_level_baseline_hazard`)
- **Time-invariant covariate slopes**: Population + N-level random effects (via `enable_ms_level_cov`)
- **Time-varying covariate coefficients**: Population only (no hierarchy)

### Baseline Hazards (GP per transition)

Each enabled transition gets its own GP baseline hazard with optional N-level hierarchy:

```stan
// parameters.stan

// --- 0→1: Progression (only if enable_ms_01=1) ---
// Population level
real log_lambda_gp_01_pop_intercept[enable_ms_01 ? 1 : 0];
real<lower=0> log_lambda_gp_01_pop_alpha[enable_ms_01 ? 1 : 0];
real<lower=0> log_lambda_gp_01_pop_rho[enable_ms_01 ? 1 : 0];
row_vector[enable_ms_01 ? max_all_t : 0] log_lambda_gp_01_pop_eta;

// N-level hierarchy (optional per level, via enable_ms_level_baseline_hazard)
// GP hyperparameters per enabled level
array[n_levels] real<lower=0> log_lambda_gp_01_level_alpha;   // sized 0 if level disabled
array[n_levels] real<lower=0> log_lambda_gp_01_level_rho;
array[n_levels] real<lower=0> log_lambda_gp_01_level_intercept_sd;
// Raw effects: [n_enabled_groups_baseline, max_all_t] using compacted indexing
matrix[n_enabled_groups_ms_baseline_01, max_all_t] log_lambda_gp_01_level_eta;
vector[n_enabled_groups_ms_baseline_01] raw_log_lambda_gp_01_level_intercept;

// --- 0→2, 1→2: Same N-level pattern ---
// ... analogous structure for each transition ...

// --- 1→2: Post-progression death ---
// For semi-Markov or extended: GP on sojourn time s (with N-level hierarchy)
real log_lambda_gp_12_s_pop_intercept[need_12_s_gp ? 1 : 0];
real<lower=0> log_lambda_gp_12_s_pop_alpha[need_12_s_gp ? 1 : 0];
real<lower=0> log_lambda_gp_12_s_pop_rho[need_12_s_gp ? 1 : 0];
row_vector[need_12_s_gp ? max_sojourn_t : 0] log_lambda_gp_12_s_pop_eta;
// N-level hierarchy follows same pattern...

// For Markov or extended: GP on clock-forward time t
// Extended uses BOTH GPs additively: log λ₁₂(t,s) = GP_t(t) + GP_s(s)
real log_lambda_gp_12_t_pop_intercept[need_12_t_gp ? 1 : 0];
real<lower=0> log_lambda_gp_12_t_pop_alpha[need_12_t_gp ? 1 : 0];
real<lower=0> log_lambda_gp_12_t_pop_rho[need_12_t_gp ? 1 : 0];
row_vector[need_12_t_gp ? max_all_t : 0] log_lambda_gp_12_t_pop_eta;
// N-level hierarchy follows same pattern...
```

### Covariate Coefficients (transition-specific)

```stan
// --- Population-level time-varying coefficients (no hierarchy) ---
// e.g., biomarker burden, velocity - computed from state-space model
vector[enable_ms_01 && enable_ms_pop_time_varying_cov ? n_time_varying_covar : 0] ms_time_varying_coef_01;
vector[enable_ms_02 && enable_ms_pop_time_varying_cov ? n_time_varying_covar : 0] ms_time_varying_coef_02;
vector[enable_ms_12 && enable_ms_pop_time_varying_cov ? n_time_varying_covar : 0] ms_time_varying_coef_12;

// --- Population-level time-invariant coefficients (QR space) ---
// Note: Post-progression specific covariates (e.g., burden at progression, time-to-progression)
// can be computed in R and passed as regular time-invariant covariates for the 1→2 transition
// e.g., baseline characteristics
vector[enable_ms_01 && enable_ms_pop_time_invariant_cov ? n_time_invariant_covar : 0] ms_time_invariant_coef_qr_01;
vector[enable_ms_02 && enable_ms_pop_time_invariant_cov ? n_time_invariant_covar : 0] ms_time_invariant_coef_qr_02;
vector[enable_ms_12 && enable_ms_pop_time_invariant_cov ? n_time_invariant_covar : 0] ms_time_invariant_coef_qr_12;

// --- Multi-level random slopes for time-invariant covariates ---
// SD hyperparameters per transition per level
array[n_levels] vector<lower=0>[enable_ms_01 ? n_covar : 0] ms_sd_level_slope_01;
array[n_levels] vector<lower=0>[enable_ms_02 ? n_covar : 0] ms_sd_level_slope_02;
array[n_levels] vector<lower=0>[enable_ms_12 ? n_covar : 0] ms_sd_level_slope_12;
// Raw effects sized by enabled groups (same indexing pattern as other_events)
```

---

## Transformed Parameters

Conditional survival matrices are computed only for enabled transitions:

```stan
// transformed_parameters.stan

// 0→1 survival (PFS event - progression or death; this is the SLD baseline)
matrix[enable_ms_01 ? n_patients : 0, enable_ms_01 ? max_all_t : 0] log_cond_surv_01;
if (enable_ms_01) {
  // ... compute GP baseline + covariates for 0→1
}

// 0→2 survival (death without progression - PSA model only)
matrix[enable_ms_02 ? n_patients : 0, enable_ms_02 ? max_all_t : 0] log_cond_surv_02;
if (enable_ms_02) {
  // ... compute GP baseline + covariates for 0→2
}

// 1→2 survival (post-progression, on sojourn time grid)
matrix[enable_ms_12 ? n_patients : 0, enable_ms_12 ? max_sojourn_t : 0] log_cond_surv_12;
if (enable_ms_12) {
  // ... compute GP baseline + covariates for 1→2
}
```

---

## Likelihood

New function in `stan/multistate.stanfunctions`:

```stan
real calc_multistate_loglik(
  // Flags
  int enable_01, int enable_02, int enable_12,
  // Patient data
  array[] int final_state,
  array[] int time_01, time_02, time_12,
  array[] int censored_01, censored_02, censored_12,
  array[] int prog_deterministic,
  // Conditional survival matrices (can be empty if transition disabled)
  matrix log_cond_surv_01,  // [n_patients, max_all_t] or [0,0]
  matrix log_cond_surv_02,  // [n_patients, max_all_t] or [0,0]
  matrix log_cond_surv_12   // [n_patients, max_sojourn_t] or [0,0]
)
```

**Key logic:**

- Check flags before computing each transition's contribution
- If `enable_01=1`: compute 0→1 survival/hazard (PFS event)
- If `enable_02=1`: compute 0→2 survival/hazard (death without progression)
- If `enable_12=1`: compute post-progression survival on sojourn time
- If `prog_deterministic=1`: no 0→1 hazard contribution at T₀₁ (SF model determines it)
- If `prog_deterministic=0`: add log-hazard at T₀₁

**SLD mode (enable_01=1, enable_02=0, enable_12=0):** Single 0→1 transition - likelihood reduces to exactly the current `other_events` computation (PFS event = progression or death).

**PSA mode (enable_01=1, enable_02=1, enable_12=1):** Full illness-death with competing risks for 0→1 and 0→2, plus post-progression death.

---

## R Data Preparation

### Flag Configuration (`r/priors.R` or targets)

```r
# SLD model (current behavior - single 0→1 PFS transition)
sld_multistate_flags <- list(
  # Transition flags
  enable_ms_01 = 1L,
  enable_ms_02 = 0L,
  enable_ms_12 = 0L,
  ms_time_scale_12 = 0L,  # unused when enable_ms_12=0

  # Hierarchy flags - all use generic n_levels structure
  enable_ms_level_baseline_hazard = c(1L, 0L),  # e.g., level 1 enabled
  enable_ms_pop_time_invariant_cov = 1L,
  enable_ms_pop_time_varying_cov = 1L,
  enable_ms_level_cov = c(0L, 0L)
)

# PSA model (full illness-death)
psa_multistate_flags <- list(
  # Transition flags
  enable_ms_01 = 1L,
  enable_ms_02 = 1L,
  enable_ms_12 = 1L,
  ms_time_scale_12 = 1L,  # semi-Markov (clock-reset)

  # Hierarchy flags - all use generic n_levels structure
  enable_ms_level_baseline_hazard = c(1L, 0L),
  enable_ms_pop_time_invariant_cov = 1L,
  enable_ms_pop_time_varying_cov = 1L,
  enable_ms_level_cov = c(0L, 0L)
)
```

### Event Classification (`r/sclc/prepare_analysis_data.R`)

```r
prepare_multistate_data <- function(analysis_data, enable_01 = FALSE) {
  # For SLD model: no 0→1 data needed, just 0→2 (death/censoring)
  # For PSA model: full state encoding

  analysis_data |>
    mutate(
      # Final state encoding (same for both models)
      ms_final_state = case_when(
        death & !any_progression ~ 2L,
        death & any_progression ~ 2L,
        any_progression ~ 1L,
        TRUE ~ 0L
      ),

      # Transition times (zeros for disabled transitions)
      ms_time_01 = if_else(enable_01 & any_progression, prog_time, 0L),
      ms_time_02 = if_else(death & !any_progression, death_week, 0L),
      ms_time_12 = if_else(death & any_progression, death_week - prog_time, 0L),

      # Censoring flags
      ms_censored_01 = if_else(enable_01, as.integer(!any_progression), 0L),
      ms_censored_02 = as.integer(!death | any_progression),
      ms_censored_12 = if_else(any_progression, as.integer(!death), 0L),

      # Deterministic progression (PSA-PD/RECIST-PD from SF model)
      ms_prog_deterministic = if_else(enable_01, as.integer(target_pd), 0L)
    )
}
```

---

## Critical Files to Modify/Create

| File | Action | Description |
|------|--------|-------------|
| `stan/modules/multistate/*.stan` | CREATE | New module (7 files) |
| `stan/multistate.stanfunctions` | CREATE | `calc_multistate_loglik()` |
| `stan/sf-ssm-log-space.stan` | MODIFY | Replace other_events includes with multistate |
| `stan/sf-ssls-lfo.stan` | MODIFY | Same changes for LFO variant |
| `r/sclc/prepare_analysis_data.R` | MODIFY | Add `prepare_multistate_data()` |
| `r/priors.R` | MODIFY | Replace other_events priors with multistate |
| `r/initializers.R` | MODIFY | Replace other_events initializers |
| `stan/modules/other_events/` | DELETE | After multistate validated |

**Reference files (patterns to follow):**

- `stan/modules/other_events/transformed_parameters.stan` - covariate computation pattern
- `stan/pfs.stanfunctions` - `calc_pch_loglik()` as likelihood template

---

## Phased Implementation

### Phase 1: Minimal Multistate (0→1 only) - Backward Compatible

- Create module skeleton with flag infrastructure
- Implement 0→1 transition only (`enable_ms_01=1, enable_ms_02=0, enable_ms_12=0`)
- This mirrors current `other_events` exactly (PFS event = progression or death)
- Keep `other_events/` temporarily for comparison

**Configuration:** SLD model uses `(1, 0, 0)` flags - single 0→1 transition, equivalent to current behavior.

**Verification:** Multistate with SLD flags produces identical results to current `other_events`. Once validated, delete `other_events/` module.

### Phase 2: Death Without Progression (0→2)

- Add 0→2 transition for death without progression
- Now progression and death are separate events
- Configuration: `(1, 1, 0)` - progression + death w/o progression, no post-progression yet

**Verification:** Short sampling runs, check parameter recovery

### Phase 3: Post-Progression Death (1→2)

- Add 1→2 hazard on sojourn time grid
- Post-progression covariates (burden at prog, time to prog)
- Full multistate likelihood with `enable_ms_12=1`

**Configuration:** Full PSA model with `(1, 1, 1)` flags.

**Verification:** Short sampling runs, check parameter recovery

### Phase 4: Extended Time Scales & Cleanup

- Markov (clock-forward) option for 1→2 (`ms_time_scale_12=0`): hazard depends on t
- Extended (both clocks) option (`ms_time_scale_12=2`): additive GP model
  - log λ₁₂(t, s) = GP₁(t) + GP₂(s) - two independent 1D GPs, no interaction term
  - Future work: consider 2D GP for interactions
- Model comparison infrastructure
- Delete `stan/modules/other_events/` after full validation
- Update all references in main model files and R code

---

## Verification Plan

1. **Syntax check:** `stanc --include-paths=stan stan/sf-ssm-log-space.stan`

2. **Unit tests:** Create `tests/testthat/test-multistate-likelihood.R`

   - Test likelihood function with known scenarios
   - Test patient encoding edge cases
   - Test each flag configuration independently

3. **Short sampling:** 10 iterations with SLD flags `(1, 0, 0)`, check no runtime errors

4. **Backward compatibility:** Run SLD model with multistate module using `(1, 0, 0)` flags, compare to current `other_events` results

5. **Full run:** Compare PFS/OS estimates between old (`other_events`) and new (multistate with SLD flags)

6. **PSA mode verification:** After Phase 3, test full `(1, 1, 1)` configuration with synthetic data
