# PSA Adaptation Plan - Implementation (Sections 1-5)

## Context

The PIONEER model currently uses SLD (Sum of Longest Diameters) for solid tumor dynamics. We need to adapt it for metastatic prostate cancer using PSA as the primary biomarker. This plan covers Sections 1-5 of PSA_260209.md, **excluding the multi-state hazard model** (Section 6). The existing other_events module will be retained.

**Key insight**: With λ=1 (Phase 1), the model is mathematically identical to the current SLD model - we're essentially renaming variables and changing progression criteria from RECIST to PCWG3.

## Design Decisions

- **Backward compatibility**: Maintain both SLD and PSA codepaths via `observation_type` flag
- **PSA covariates**: Include nadir, nadir_ratio in Phase 1 (clinically important for PCWG3)
- **Test data**: PSA time course data (longitudinal PSA measurements per patient in weeks) - to be provided
- **Future-proof for λ**: Structure code so λ parameter can be introduced later (Phase 3)
- **Future-proof for dual observations**: Support observation_type=2 (both SLD + PSA) for estimating λ and φ (Phase 5)

---

## Stan Directory Structure (Updated 2026-02-10)

The Stan codebase has been reorganized into a modular structure. Key directories:

```
stan/
├── _base_data.stan              # Core patient/visit structure (shared by all biomarkers)
├── _base_transformed_data.stan  # Base preprocessing
├── sf-ssm-log-space.stan        # Main model entry point
├── sf-ssls-lfo.stan             # Leave-future-out variant
├── sf-ssls-lfo-endpoints.stan   # LFO endpoints variant
│
├── *.stanfunctions              # Function libraries (util, pos, gp, pfs, lfo)
│
└── modules/
    ├── state_space/             # State-space infrastructure (orchestrates parameter modules)
    │   ├── functions.stanfunctions  # Core SF dynamics (~1900 lines)
    │   ├── data.stan            # Outcome configuration (quantiles, timepoints)
    │   ├── transformed_data.stan
    │   ├── transformed_parameters.stan
    │   ├── generated_quantities.stan
    │   ├── checks.stan
    │   └── lfo_data.stan
    │
    ├── tumor/                   # Tumor/SLD-specific module
    │   ├── data.stan            # SLD measurements, RECIST, PFS outcomes
    │   ├── transformed_data.stan
    │   └── tumor.stanfunctions  # calc_log_sld_mean, calculate_target_recist
    │
    ├── tr/                      # Total rate parameters
    ├── frac/                    # Growth fraction parameters
    ├── init/                    # Initial state parameters
    ├── measurement/             # Measurement model (σ_meas)
    └── other_events/            # Non-tumor progression hazards
```

**Naming conventions**:
- `.stanfunctions` extension for files containing only function definitions
- `_` prefix for base files included at root level
- Module files follow 7-file pattern: `flags.stan`, `data.stan`, `hyperparams.stan`, `transformed_data.stan`, `parameters.stan`, `transformed_parameters.stan`, `priors.stan`

---

## Phase 1: Dual-Mode Model (SLD + PSA, λ=1)

### 1.1 New PSA Module Structure

Create a new `modules/psa/` module following the established pattern:

```
stan/modules/psa/
├── data.stan                # PSA measurements and outcomes
├── transformed_data.stan    # PSA normalization, visit indices
├── psa.stanfunctions        # calc_log_psa_mean, calculate_psa_category (PCWG3)
├── flags.stan               # fit_psa_data, observation_type flags
└── hyperparams.stan         # PSA measurement error priors, lambda_fixed
```

### 1.2 Stan Data Input Changes

**File**: `stan/modules/psa/flags.stan` (NEW)

```stan
// Observation type flag (controls which biomarker likelihoods are active)
// 0 = SLD only (tumor module)
// 1 = PSA only (psa module)
// 2 = Both SLD + PSA (joint modeling)
int<lower=0, upper=2> observation_type;

// Whether to fit PSA observation model
int<lower=0, upper=1> fit_psa_data;

// Lambda parameter control (for future Phase 3/5)
int<lower=0, upper=1> estimate_lambda;  // 0=fix at lambda_fixed, 1=estimate
```

**File**: `stan/modules/psa/hyperparams.stan` (NEW)

```stan
// Fixed lambda value when estimate_lambda=0 (default: 1.0)
real<lower=0> lambda_fixed;

// PSA measurement error prior hyperparameters
real<lower=0> measure_sd_psa_mean;
real<lower=0> measure_sd_psa_sd;
```

**File**: `stan/modules/psa/data.stan` (NEW)

```stan
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

// PCWG3 response category at each visit
// 1 = Undetectable (CR equivalent)
// 2 = PSA50 (PR equivalent)
// 3 = Stable (SD equivalent)
// 4 = PSA-PD (confirmed progression)
// 5 = Not Evaluable
array[sum(n_patient_visits)] int<lower=1, upper=5> pcwg3_category;

// PSA undetectable threshold (typically 0.1 ng/mL)
real<lower=0> psa_undetectable_threshold;
```

**Identifiability table** (from PSA_260209.md Section 2.4):

| observation_type | estimate_lambda | What's identifiable |
|------------------|-----------------|---------------------|
| 0 (SLD only) | N/A | B(t), φ, d, g |
| 1 (PSA only) | 0 (λ=1) | B(t) = B_PSA(t), φ, d, g |
| 1 (PSA only) | 0 (λ≠1 fixed) | B_PSA(t), α_PSA, d, g (NOT φ, λ separately) |
| 2 (both) | 1 | B(t), B_PSA(t), φ, λ, d, g |

### 1.3 Observation Model Changes

**File**: `stan/modules/psa/transformed_data.stan` (NEW)

```stan
// ============================================================================
// PSA NORMALIZATION AND PREPROCESSING
// ============================================================================

// Normalized PSA (relative to baseline)
vector[sum(n_patient_visits)] normalized_psa;

// Log-transformed PSA values
vector[sum(n_patient_visits)] log_psa_values;
array[n_patients] real log_baseline_psa;

// Compute normalized PSA per patient
for (i in 1:n_patients) {
  int visit_start, visit_end;
  (visit_start, visit_end) = get_pos(patient_visit_pos, i);

  // Baseline is first post-treatment visit (at screening boundary)
  int baseline_idx = visit_start + n_screening_visits[i] - 1;
  real baseline_psa = psa_values[baseline_idx];
  log_baseline_psa[i] = log(baseline_psa);

  for (v in visit_start:visit_end) {
    if (psa_measured[v]) {
      normalized_psa[v] = psa_values[v] / baseline_psa;
      log_psa_values[v] = log(psa_values[v]);
    }
  }
}
```

**File**: `stan/modules/psa/psa.stanfunctions` (NEW)

```stan
// ============================================================================
// PSA-SPECIFIC FUNCTIONS
// ============================================================================
// These functions are specific to PSA measurements. They convert
// the general 2-component log-space state (from state_space module) into
// PSA-specific quantities.
//
// Follows the same pattern as modules/tumor/tumor.stanfunctions

/**
 * Calculate mean log(PSA) from two-component log-space states
 *
 * When lambda=1, this is identical to calc_log_sld_mean.
 * When lambda!=1, applies the PSA-weighted burden formula:
 *   B_PSA(t) = [φ × e^(-d×t) + λ × (1-φ) × e^(g×t)] / [φ + λ × (1-φ)]
 *
 * @param patient_states Matrix of states [n_visits × 2]
 * @param psa_baseline Baseline PSA measurement
 * @param lambda PSA production ratio (resistant/sensitive), default 1.0
 * @param log_phi Log of sensitive fraction at baseline
 * @return Vector of log(PSA) means at each visit
 */
vector calc_log_psa_mean(matrix patient_states, real psa_baseline, real lambda, real log_phi) {
  assert_equal(cols(patient_states), 2);
  int n_visits = rows(patient_states);
  vector[n_visits] log_psa_mean;

  if (lambda == 1.0) {
    // Simplified: identical to SLD model
    log_psa_mean = to_vector(log_sum_exp(patient_states[, 1], patient_states[, 2]))
                   + log(psa_baseline);
  } else {
    // Full PSA-weighted burden formula
    real log_lambda = log(lambda);
    real log1m_phi = log1m_exp(log_phi);
    real normalizer = log_sum_exp(log_phi, log_lambda + log1m_phi);

    for (t in 1:n_visits) {
      real x_s = patient_states[t, 1];  // log(sensitive burden)
      real x_g = patient_states[t, 2];  // log(resistant burden)
      log_psa_mean[t] = log_sum_exp(x_s, log_lambda + x_g) - normalizer + log(psa_baseline);
    }
  }

  return log_psa_mean;
}

/**
 * Calculate PCWG3 PSA response category with confirmation logic
 *
 * Categories:
 *   1 = Undetectable (CR equivalent): PSA < threshold
 *   2 = PSA50 (PR equivalent): ≥50% decrease from baseline
 *   3 = Stable (SD equivalent): Neither PSA50 nor PSA-PD
 *   4 = PSA-PD: ≥25% AND ≥2 ng/mL from nadir, confirmed ≥3 weeks later
 *
 * @param psa_trajectory Absolute PSA values (ng/mL) at each visit
 * @param pre_nadir Pre-existing nadir from previous visits (0 if none)
 * @param n_screening Number of screening visits before treatment
 * @param psa_undetectable_threshold Threshold for undetectable (typically 0.1 ng/mL)
 * @return Array of PCWG3 categories for post-treatment visits
 */
array[] int calculate_psa_category(
  vector psa_trajectory,
  real pre_nadir,
  int n_screening,
  real psa_undetectable_threshold
) {
  int n = rows(psa_trajectory);
  assert_greater_than_or_equal(n, n_screening);
  int n_treat = n - n_screening;
  array[n_treat] int psa_status;

  real baseline = psa_trajectory[n_screening];
  real nadir = pre_nadir <= 0 ? baseline : pre_nadir;
  int pending_pd_week = -1;

  int CR = 1; int PSA50 = 2; int STABLE = 3; int PSA_PD = 4;

  for (t in (n_screening + 1):n) {
    real psa_t = psa_trajectory[t];
    nadir = fmin(nadir, psa_t);

    // Check PCWG3 progression: ≥25% AND ≥2 ng/mL from nadir
    int meets_pd_criteria = (nadir > 0) &&
                           ((psa_t - nadir) / nadir >= 0.25) &&
                           ((psa_t - nadir) >= 2.0);

    // Handle confirmation (simplified: check if criteria met for 3+ weeks)
    if (meets_pd_criteria) {
      if (pending_pd_week > 0 && (t - pending_pd_week) >= 3) {
        psa_status[t - n_screening] = PSA_PD;  // Confirmed
      } else {
        if (pending_pd_week < 0) pending_pd_week = t;
        // Not yet confirmed - check other categories
        if (psa_t < psa_undetectable_threshold) {
          psa_status[t - n_screening] = CR;
        } else if (baseline > 0 && (baseline - psa_t) / baseline >= 0.5) {
          psa_status[t - n_screening] = PSA50;
        } else {
          psa_status[t - n_screening] = STABLE;
        }
      }
    } else {
      pending_pd_week = -1;  // Reset confirmation
      if (psa_t < psa_undetectable_threshold) {
        psa_status[t - n_screening] = CR;
      } else if (baseline > 0 && (baseline - psa_t) / baseline >= 0.5) {
        psa_status[t - n_screening] = PSA50;
      } else {
        psa_status[t - n_screening] = STABLE;
      }
    }
  }
  return psa_status;
}
```

### 1.4 Measurement Model Changes

**File**: `stan/modules/measurement/parameters.stan` (ADD)

```stan
// PSA measurement error (when observation_type >= 1)
real<lower=0> measure_sd_psa;
```

**File**: `stan/modules/measurement/priors.stan` (ADD)

```stan
// PSA measurement error prior (when observation_type >= 1)
if (observation_type >= 1) {
  measure_sd_psa ~ normal(measure_sd_psa_mean, measure_sd_psa_sd);
}
```

### 1.5 Other Events Covariate Changes

**File**: `stan/modules/other_events/transformed_data.stan`

Variable renaming for biomarker-agnostic naming:
- `median_log_sld_obs` → `median_log_biomarker_obs`
- `iqr_log_sld_obs` → `iqr_log_biomarker_obs`
- `log_baseline_sld` → `log_baseline_biomarker`

**File**: `stan/modules/other_events/transformed_parameters.stan`

Variable renaming:
- `log_sld_normalized` → `log_biomarker_normalized`
- `log_sld_absolute` → `log_biomarker_absolute`
- `log_sld_standardized` → `log_biomarker_standardized`
- `sld_velocity` → `biomarker_velocity`

### 1.6 PSA-Derived Covariates (Nadir, Nadir Ratio)

**File**: `stan/modules/other_events/transformed_parameters.stan` (ADD)

Add after velocity computation:

```stan
// PSA-specific covariates (only when observation_type >= 1)
// These use model-predicted burden (not observed PSA) for smoothness
if (observation_type >= 1 && n_tumor_covar >= 5) {
  for (k in 1:n_causes) {
    for (i in 1:n_patients) {
      // Compute nadir and nadir ratio from predicted burden
      row_vector[max_all_t] nadir_ratio = rep_row_vector(1, max_all_t);

      real running_nadir = exp(log_biomarker_absolute[i, 1]);
      for (t in 1:max_all_t) {
        running_nadir = fmin(running_nadir, exp(log_biomarker_absolute[i, t]));
        nadir_ratio[t] = exp(log_biomarker_absolute[i, t]) / running_nadir;
      }

      // Add nadir ratio effect (covariate index 5)
      oe_time_varying_log_hazard_ratio[k, i] += oe_tumor_coef_pop[k][5] * log(nadir_ratio);
    }
  }
}
```

### 1.7 Main Model File Updates

**File**: `stan/sf-ssm-log-space.stan` (MODIFY)

Add PSA module includes:

```stan
functions {
  #include "util.stanfunctions"
  #include "pos.stanfunctions"
  #include "gp.stanfunctions"
  #include "pfs.stanfunctions"
  #include "lfo.stanfunctions"
  #include "modules/state_space/functions.stanfunctions"
  #include "modules/tumor/tumor.stanfunctions"
  #include "modules/psa/psa.stanfunctions"  // NEW
}

data {
  #include "_base_data.stan"
  #include "modules/tumor/data.stan"
  #include "modules/psa/data.stan"          // NEW
  #include "modules/psa/flags.stan"         // NEW
  #include "modules/psa/hyperparams.stan"   // NEW
  #include "modules/state_space/data.stan"
  // ... rest of includes
}

transformed data {
  #include "_base_transformed_data.stan"
  #include "modules/tumor/transformed_data.stan"
  #include "modules/psa/transformed_data.stan"  // NEW
  // ... rest of includes
}
```

**File**: `stan/sf-ssm-log-space.stan` model block (MODIFY)

Add conditional PSA likelihood:

```stan
model {
  // ... existing priors ...

  profile("loglik") {
    // SLD likelihood (observation_type 0 or 2)
    if (fit_tumor_data && observation_type != 1) {
      profile("tumor loglik") {
        for (i in 1:n_patients) {
          int visit_start, visit_end;
          (visit_start, visit_end) = get_pos(patient_visit_pos, i);
          normalized_sld[visit_start:visit_end] ~ sf_log_space_obs(
            states[visit_start:visit_end], measure_sd, log_lod - log_baseline_sld[i]
          );
        }
      }
    }

    // PSA likelihood (observation_type 1 or 2)
    if (fit_psa_data && observation_type >= 1) {
      profile("psa loglik") {
        for (i in 1:n_patients) {
          int visit_start, visit_end;
          (visit_start, visit_end) = get_pos(patient_visit_pos, i);
          for (v in visit_start:visit_end) {
            if (psa_measured[v]) {
              // When lambda=1, this is mathematically identical to SLD model
              real log_psa_mean = calc_log_psa_mean(
                states[v:v], psa_values[visit_start], lambda_fixed,
                patient_log_decrease_prop[i]
              )[1];
              log(normalized_psa[v]) ~ normal(log_psa_mean - log_baseline_psa[i], measure_sd_psa);
            }
          }
        }
      }
    }

    // ... other events likelihood ...
  }
}
```

### 1.8 R Data Preparation Changes

**File**: `r/sclc/prepare_analysis_data.R`

Add PSA data handling:

```r
prepare_tumor_stan_data <- function(
  analysis_data,
  ...,
  observation_type = 0,      # 0=SLD, 1=PSA, 2=both
  estimate_lambda = 0,       # 0=use lambda_fixed, 1=estimate (requires observation_type=2)
  lambda_fixed = 1.0         # Fixed value when estimate_lambda=0
) {
  # Validate: can only estimate lambda when both observations available
  if (estimate_lambda == 1 && observation_type != 2) {
    stop("estimate_lambda=1 requires observation_type=2 (both SLD and PSA)")
  }

  # SLD data (observation_type 0 or 2)
  if (observation_type != 1) {
    sld_values <- unnest(analysis_data, visit_data) |>
      pull(mmsumdiam) |> magrittr::divide_by(10)
    fit_tumor_data <- 1L
  } else {
    sld_values <- numeric(0)
    fit_tumor_data <- 0L
  }

  # PSA data (observation_type 1 or 2)
  if (observation_type >= 1) {
    psa_values <- unnest(analysis_data, visit_data) |>
      pull(psa)  # PSA in ng/mL
    psa_measured <- as.integer(!is.na(psa_values))
    psa_values[is.na(psa_values)] <- 0  # Fill NAs for Stan
    fit_psa_data <- 1L
  } else {
    psa_values <- numeric(0)
    psa_measured <- integer(0)
    fit_psa_data <- 0L
  }

  lst(
    observation_type = observation_type,
    estimate_lambda = estimate_lambda,
    lambda_fixed = lambda_fixed,
    fit_tumor_data = fit_tumor_data,
    fit_psa_data = fit_psa_data,
    sum_tumor_size = sld_values,
    psa_values = psa_values,
    psa_measured = psa_measured,
    psa_undetectable_threshold = 0.1,  # ng/mL
    ...
  )
}
```

**Phase 1 usage** (PSA only, λ=1):
```r
stan_data <- prepare_tumor_stan_data(
  analysis_data,
  observation_type = 1,      # PSA only
  estimate_lambda = 0,       # Don't estimate
  lambda_fixed = 1.0         # λ = 1 (mathematically identical to SLD model)
)
```

**Future Phase 5 usage** (both SLD + PSA, estimate λ):
```r
stan_data <- prepare_tumor_stan_data(
  analysis_data,
  observation_type = 2,      # Both SLD and PSA
  estimate_lambda = 1        # Estimate λ from data
)
```

---

## Future Phases (NOT in current scope, but code structured to support)

### Phase 3: Hierarchical λ Parameter (PSA-only, λ≠1)

When `observation_type=1` and `estimate_lambda=0` but `lambda_fixed≠1`:

**Use case**: Apply an informative prior or fixed value for λ based on external knowledge about PSA production rates in resistant vs sensitive cells.

**Stan parameter** (add to `stan/modules/psa/parameters.stan`):
```stan
// When estimate_lambda=1 (requires observation_type=2)
real<lower=0> lambda_psa;  // PSA production ratio (resistant/sensitive)
```

**Prior** (add to `r/priors.R`):
```stan
log(lambda_psa) ~ Normal(0, 0.5);  // Centered at λ=1, allows 0.6 to 1.6
```

**Important**: With PSA-only data, λ and φ are **not separately identifiable**. The model estimates α_PSA (PSA-weighted sensitive fraction) where:
```
α_PSA = φ / [φ + λ × (1 - φ)]
```

### Phase 5: Dual Observations (SLD + PSA, estimate λ and φ)

When `observation_type=2` and `estimate_lambda=1`:

**Use case**: Dataset has both imaging (SLD) and PSA measurements. This enables estimation of true tumor burden B(t), true sensitive fraction φ, AND the PSA production ratio λ.

**Observation model**:
```stan
// SLD directly observes true burden
log(SLD/SLD_base) ~ Normal(log(B(t)), σ_sld)

// PSA observes λ-weighted burden
log(PSA/PSA_base) ~ Normal(log(B_PSA(t)), σ_psa)

// Where B_PSA(t) = [φ × e^(-d×t) + λ × (1-φ) × e^(g×t)] / [φ + λ × (1-φ)]
```

**What becomes identifiable**:
- φ (true sensitive fraction) - from SLD
- λ (PSA production ratio) - from SLD+PSA comparison
- d, g (rates) - from both
- Separate σ_sld and σ_psa (measurement noise)

**Data requirements**:
- Same patients must have both SLD and PSA observations
- Observations don't need to be at same timepoints (model interpolates)

---

## Files to Modify

### Phase 1 (Current Scope)

| File | Changes |
|------|---------|
| `stan/modules/psa/data.stan` | **NEW**: PSA measurements, PCWG3 categories, outcomes |
| `stan/modules/psa/flags.stan` | **NEW**: observation_type, fit_psa_data, estimate_lambda |
| `stan/modules/psa/hyperparams.stan` | **NEW**: lambda_fixed, measure_sd_psa priors |
| `stan/modules/psa/transformed_data.stan` | **NEW**: PSA normalization, preprocessing |
| `stan/modules/psa/psa.stanfunctions` | **NEW**: calc_log_psa_mean, calculate_psa_category |
| `stan/modules/measurement/parameters.stan` | Add measure_sd_psa |
| `stan/modules/measurement/priors.stan` | Add measure_sd_psa prior |
| `stan/modules/other_events/transformed_data.stan` | Rename sld → biomarker variables |
| `stan/modules/other_events/transformed_parameters.stan` | Rename + add nadir covariates |
| `stan/sf-ssm-log-space.stan` | Add PSA module includes + conditional likelihood |
| `stan/sf-ssls-lfo.stan` | Add PSA module includes |
| `stan/sf-ssls-lfo-endpoints.stan` | Add PSA module includes |
| `stan/_endpoints_generated_quantities.stan` | Use PCWG3 for PSA mode |
| `r/sclc/prepare_analysis_data.R` | Add PSA data handling |
| `r/priors.R` | Add priors for PSA measurement error |

### Future Phases (Structure now, implement later)

| File | Phase | Changes |
|------|-------|---------|
| `stan/modules/psa/parameters.stan` | 3, 5 | Add `lambda_psa` parameter |
| `stan/modules/psa/priors.stan` | 3, 5 | Add `log(lambda_psa) ~ Normal(0, 0.5)` prior |
| `stan/modules/psa/psa.stanfunctions` | 3, 5 | Update calc_log_psa_mean for λ≠1 |
| `r/priors.R` | 3, 5 | Add `lambda_psa` prior hyperparameters |
| `r/initializers.R` | 3, 5 | Add `lambda_psa` initializer |

---

## Staged Testing Plan

### Stage 1: Syntax Validation (No data required)

```bash
# Check Stan syntax after changes (use latest cmdstan)
~/.cmdstan/cmdstan-2.38.0/bin/stanc --include-paths=stan stan/sf-ssm-log-space.stan
```

**Pass criteria**: No syntax errors

### Stage 2: Unit Tests for PCWG3 Functions

Create `tests/testthat/test-stan-pcwg3.R`:

```r
test_that("PCWG3 PSA categories computed correctly", {
  # Test PSA50 response
  psa_trajectory <- c(100, 100, 45, 40, 35)  # >50% decline
  result <- calculate_psa_category_r(psa_trajectory, n_screening = 2)
  expect_equal(result[3], 2)  # PSA50

  # Test PSA-PD (confirmed progression)
  psa_trajectory <- c(100, 100, 50, 40, 55, 60)  # nadir=40, then rise
  result <- calculate_psa_category_r(psa_trajectory, n_screening = 2)
  # Week 5: (55-40)/40 = 37.5% > 25%, 55-40=15 > 2 ng/mL
  # Week 6: confirmed (>3 weeks)
  expect_equal(result[5], 4)  # PSA-PD
})
```

**Pass criteria**: All PCWG3 logic tests pass

### Stage 3: Data Preparation Test

```r
# Test with mock PSA data
mock_psa_data <- tibble(
  psa = c(100, 80, 60, 50, 55, 70),
  week = c(-2, 0, 4, 8, 12, 16)
)

stan_data <- prepare_tumor_stan_data(
  analysis_data = mock_analysis_data,
  observation_type = 1,    # PSA mode
  estimate_lambda = 0,     # Don't estimate (Phase 1)
  lambda_fixed = 1.0       # λ = 1
)

expect_true("psa_values" %in% names(stan_data))
expect_equal(length(stan_data$psa_values), 6)
expect_equal(stan_data$observation_type, 1)
expect_equal(stan_data$estimate_lambda, 0)
expect_equal(stan_data$lambda_fixed, 1.0)

# Test validation: estimate_lambda=1 should fail without observation_type=2
expect_error(
  prepare_tumor_stan_data(mock_analysis_data, observation_type = 1, estimate_lambda = 1),
  "requires observation_type=2"
)
```

**Pass criteria**: Stan data structure correct for PSA, lambda parameters validated

### Stage 4: Model Compilation

```bash
# Compile modified model via targets pipeline
./sclc_targets.sh -k -m '"tumor_ssls_exe_hash"'
```

**Pass criteria**: Model compiles without errors

### Stage 5: Short Sampling Test (10 iterations)

```r
# Quick validation that model runs
model$sample(
  data = psa_stan_data,
  chains = 1,
  iter_warmup = 10,
  iter_sampling = 10
)
```

**Pass criteria**: No runtime errors, parameters initialized correctly

### Stage 6: Backward Compatibility Test

```r
# Run with observation_type=0 (SLD mode)
# Should produce identical results to current model
stan_data_sld <- prepare_tumor_stan_data(observation_type = 0)
fit_sld <- model$sample(data = stan_data_sld, ...)

# Compare with baseline run
expect_equal(
  fit_sld$summary("tr_intercept_pop")$mean,
  baseline_fit$summary("tr_intercept_pop")$mean,
  tolerance = 0.01
)
```

**Pass criteria**: SLD mode unchanged from baseline

### Stage 7: Full PSA Data Run

```r
# Load PSA time course data (to be provided)
# Expected format: longitudinal PSA measurements per patient
# Columns: patient_id, week, psa (ng/mL), [covariates]
psa_data <- read_csv("path/to/psa_timecourse.csv")

# Prepare analysis data
psa_analysis_data <- prepare_psa_analysis_data(psa_data)

# Full MCMC run
fit_psa <- model$sample(
  data = prepare_tumor_stan_data(psa_analysis_data, observation_type = 1),
  chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000
)

# Check convergence
expect_true(all(fit_psa$summary()$rhat < 1.1))
```

**Pass criteria**: Model converges on PSA time course data

---

## Verification Checklist

### Phase 1 (Current Scope)

- [ ] Stan syntax check passes
- [ ] PCWG3 unit tests pass
- [ ] R data preparation handles PSA (observation_type=1)
- [ ] R data preparation validates estimate_lambda requires observation_type=2
- [ ] Model compiles with new data fields
- [ ] Short sampling runs without errors (observation_type=1, lambda_fixed=1)
- [ ] SLD mode backward compatible (observation_type=0)
- [ ] Full PSA run converges

### Future Phases (Verify structure only)

- [ ] Data structure supports observation_type=2 (both SLD + PSA)
- [ ] estimate_lambda flag present in Stan data
- [ ] lambda_fixed value passed through correctly
- [ ] Code comments/TODO markers in place for Phase 3/5 implementation

---

## Document History

- **2026-02-10 (Implementation)**: Actual implementation deviated from plan in several key ways:

  **1. No flags.stan in PSA module**
  - Plan: Create `stan/modules/psa/flags.stan` with observation_type, fit_psa_data, estimate_lambda flags
  - Actual: Did NOT create flags.stan. Flags like observation_type belong at a higher level (main model files), not inside pluggable observation modules. This makes modules truly pluggable.

  **2. Removed measurement/ module entirely**
  - Plan: Add measure_sd_psa to existing `stan/modules/measurement/` module
  - Actual: Removed measurement/ module completely. Measurement error now lives inside each observation module:
    - `tumor/hyperparams.stan`, `tumor/parameters.stan`, `tumor/priors.stan` for SLD
    - `psa/hyperparams.stan`, `psa/parameters.stan`, `psa/priors.stan` for PSA
  - Rationale: Each observation module should be self-contained with its own measurement model

  **3. Simplified calc_log_psa_mean for λ=1**
  - Plan: Function signature included lambda and log_phi parameters
  - Actual: For Phase 1 (λ=1), simplified to `calc_log_psa_mean(patient_states, psa_baseline)` without lambda/phi args
  - Rationale: λ=1 makes model mathematically identical to SLD; no need for extra parameters

  **4. Added psa_interval_censored**
  - Plan: Only psa_pfs and psa_right_censored in data.stan
  - Actual: Also added `psa_interval_censored` for consistency with other censoring patterns

  **5. Main model files not yet modified**
  - Plan: Update sf-ssm-log-space.stan, sf-ssls-lfo.stan, sf-ssls-lfo-endpoints.stan with PSA includes
  - Actual: Only updated to use tumor/ module (moved from measurement/). PSA includes deferred to Phase 2.

  **Files created during implementation:**
  - `stan/modules/psa/data.stan` ✓
  - `stan/modules/psa/transformed_data.stan` ✓
  - `stan/modules/psa/psa.stanfunctions` ✓
  - `stan/modules/psa/hyperparams.stan` ✓
  - `stan/modules/psa/parameters.stan` ✓
  - `stan/modules/psa/priors.stan` ✓
  - `stan/modules/tumor/hyperparams.stan` ✓ (moved from measurement/)
  - `stan/modules/tumor/parameters.stan` ✓ (moved from measurement/)
  - `stan/modules/tumor/priors.stan` ✓ (moved from measurement/)

  **Files removed:**
  - `stan/modules/measurement/` (entire directory)

- **2026-02-10**: Updated file paths and module structure to reflect Stan reorganization
  - Changed from `stan/tumor/base_data.stan` to `stan/_base_data.stan`
  - Changed from `stan/ssls/_sf_*.stan` to `stan/modules/state_space/*.stan`
  - Changed from `stan/pcwg3.stanfunctions` to `stan/modules/psa/psa.stanfunctions`
  - Added complete PSA module structure following established patterns
  - Updated cmdstan version in test commands (2.37.0 → 2.38.0)
- **2026-02-09**: Initial implementation plan created
