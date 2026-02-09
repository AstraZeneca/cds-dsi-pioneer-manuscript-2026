# PSA Adaptation Plan - Implementation (Sections 1-5)

## Context

The PIONEER model currently uses SLD (Sum of Longest Diameters) for solid tumor dynamics. We need to adapt it for metastatic prostate cancer using PSA as the primary biomarker. This plan covers Sections 1-5 of PSA_260209.md, **excluding the multi-state hazard model** (Section 6). The existing other_events module will be retained.

**Key insight**: With λ=1 (Phase 1), the model is mathematically identical to the current SLD model - we're essentially renaming variables and changing progression criteria from RECIST to PCWG3.

## Design Decisions

- **Backward compatibility**: Maintain both SLD and PSA codepaths via `observation_type` flag
- **PSA covariates**: Include nadir, nadir_ratio in Phase 1 (clinically important for PCWG3)
- **Test data**: PSA time course data (longitudinal PSA measurements per patient in weeks) - to be provided

---

## Phase 1: Dual-Mode Model (SLD + PSA, λ=1)

### 1.1 Stan Data Input Changes

**File**: `stan/tumor/base_data.stan` (or equivalent data declaration)

```stan
// Add observation type flag
int<lower=0, upper=2> observation_type;  // 0=SLD, 1=PSA, 2=both

// Rename for clarity (backward compatible)
// sum_tumor_size → sum_biomarker (or keep sum_tumor_size for SLD, add psa_values for PSA)
vector[sum(n_patient_visits)] biomarker_values;  // PSA in ng/mL or SLD in cm
```

**Changes needed**:
- Add `observation_type` flag to Stan data
- Keep `sum_tumor_size` for SLD compatibility
- Add `psa_values` array for PSA data (when observation_type >= 1)

### 1.2 Observation Model Changes

**File**: `stan/ssls/_sf_transformed_data.stan` (lines 12-18)

Current:
```stan
normalized_sld[visit_pos:visit_end] = sum_tumor_size[visit_pos:visit_end] / sum_tumor_size[visit_pos];
```

Change to:
```stan
// Generalized normalization (works for both SLD and PSA)
if (observation_type == 0) {
  // SLD mode (current behavior)
  normalized_biomarker[visit_pos:visit_end] = sum_tumor_size[...] / sum_tumor_size[visit_pos];
} else {
  // PSA mode
  normalized_biomarker[visit_pos:visit_end] = psa_values[...] / psa_values[visit_pos];
}
```

**File**: `stan/ssls/_sf_functions.stan` (line 114, `sf_log_space_obs_lpdf`)

No changes needed - function already works with normalized observations.

### 1.3 PCWG3 Progression Criteria

**New file**: `stan/pcwg3.stanfunctions`

Create PCWG3 equivalent of `stan/recist.stanfunctions`:

```stan
/**
 * Calculate PCWG3 PSA response category
 * Returns: 1=CR (undetectable), 2=PSA50 (PR), 3=Stable (SD), 4=PSA-PD
 */
array[] int calculate_psa_category(
  vector psa_trajectory,      // Absolute PSA values (ng/mL)
  int n_screening,
  real psa_undetectable_threshold  // 0.1 ng/mL
) {
  int n = rows(psa_trajectory);
  int n_treat = n - n_screening;
  array[n_treat] int psa_status;

  real baseline = psa_trajectory[n_screening];
  real nadir = baseline;
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
        } else if ((baseline - psa_t) / baseline >= 0.5) {
          psa_status[t - n_screening] = PSA50;
        } else {
          psa_status[t - n_screening] = STABLE;
        }
      }
    } else {
      pending_pd_week = -1;  // Reset confirmation
      if (psa_t < psa_undetectable_threshold) {
        psa_status[t - n_screening] = CR;
      } else if ((baseline - psa_t) / baseline >= 0.5) {
        psa_status[t - n_screening] = PSA50;
      } else {
        psa_status[t - n_screening] = STABLE;
      }
    }
  }
  return psa_status;
}
```

### 1.4 Other Events Covariate Changes

**File**: `stan/ssls/modules/other_events/transformed_data.stan` (lines 51-87)

Variable renaming for clarity:
- `median_log_sld_obs` → `median_log_biomarker_obs`
- `iqr_log_sld_obs` → `iqr_log_biomarker_obs`
- `log_baseline_sld` → `log_baseline_biomarker`

**File**: `stan/ssls/modules/other_events/transformed_parameters.stan` (lines 60-137)

Variable renaming:
- `log_sld_normalized` → `log_biomarker_normalized`
- `log_sld_absolute` → `log_biomarker_absolute`
- `log_sld_standardized` → `log_biomarker_standardized`
- `sld_velocity` → `biomarker_velocity`

### 1.5 R Data Preparation Changes

**File**: `r/sclc/prepare_analysis_data.R`

Add PSA data handling:

```r
prepare_tumor_stan_data <- function(analysis_data, ..., observation_type = 0) {
  # observation_type: 0=SLD, 1=PSA, 2=both

  if (observation_type == 0) {
    # Current SLD logic
    biomarker_values <- unnest(analysis_data, visit_data) |>
      pull(mmsumdiam) |> divide_by(10)
  } else {
    # PSA logic
    biomarker_values <- unnest(analysis_data, visit_data) |>
      pull(psa)  # PSA in ng/mL
  }

  lst(
    observation_type = observation_type,
    sum_tumor_size = if (observation_type == 0) biomarker_values else numeric(0),
    psa_values = if (observation_type >= 1) biomarker_values else numeric(0),
    ...
  )
}
```

### 1.6 PSA-Derived Covariates (Nadir, Nadir Ratio)

**File**: `stan/ssls/modules/other_events/transformed_parameters.stan`

Add after velocity computation (around line 134):

```stan
// PSA-specific covariates (only when observation_type >= 1)
if (observation_type >= 1 && n_tumor_covar >= 5) {
  for (k in 1:n_causes) {
    for (i in 1:n_patients) {
      // Compute nadir and nadir ratio
      row_vector[max_all_t] nadir = rep_row_vector(1e10, max_all_t);
      row_vector[max_all_t] nadir_ratio = rep_row_vector(1, max_all_t);

      real running_nadir = exp(log_biomarker_absolute[1]);
      for (t in 1:max_all_t) {
        running_nadir = fmin(running_nadir, exp(log_biomarker_absolute[t]));
        nadir[t] = running_nadir;
        nadir_ratio[t] = exp(log_biomarker_absolute[t]) / running_nadir;
      }

      // Add nadir ratio effect (covariate index 5)
      oe_time_varying_log_hazard_ratio[k, i] += oe_tumor_coef_pop[k][5] * log(nadir_ratio);
    }
  }
}
```

**File**: `r/sclc/prepare_analysis_data.R` (line 134)

```r
n_tumor_covar <- if (observation_type == 0) 3 else 5
# PSA adds: nadir, nadir_ratio (time_to_nadir can be derived)
```

---

## Future Phase: Hierarchical λ Parameter (NOT in current scope)

This phase is **deferred**. When implemented:
- Add `lambda_psa` parameter with prior `log(λ) ~ Normal(0, 0.5)`
- Modify observation model to use B_PSA(t) formula
- Only identifiable when both PSA and SLD available

---

## Files to Modify

| File | Changes |
|------|---------|
| `stan/tumor/base_data.stan` | Add `observation_type`, `psa_values` |
| `stan/ssls/_sf_transformed_data.stan` | Conditional normalization (SLD vs PSA) |
| `stan/ssls/modules/other_events/transformed_data.stan` | Rename sld → biomarker |
| `stan/ssls/modules/other_events/transformed_parameters.stan` | Rename + add nadir covariates |
| `stan/pcwg3.stanfunctions` | **NEW**: PCWG3 category functions |
| `stan/ssls/_endpoints_generated_quantities.stan` | Use PCWG3 for PSA mode |
| `r/sclc/prepare_analysis_data.R` | Add PSA data handling |
| `r/priors.R` | Add priors for new PSA covariates |

---

## Staged Testing Plan

### Stage 1: Syntax Validation (No data required)

```bash
# Check Stan syntax after changes
~/.cmdstan/cmdstan-2.37.0/bin/stanc --include-paths=stan,stan/ssls stan/ssls/sf-ssm-log-space.stan
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
  observation_type = 1  # PSA mode
)

expect_true("psa_values" %in% names(stan_data))
expect_equal(length(stan_data$psa_values), 6)
```

**Pass criteria**: Stan data structure correct for PSA

### Stage 4: Model Compilation

```bash
# Compile modified model
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

- [ ] Stan syntax check passes
- [ ] PCWG3 unit tests pass
- [ ] R data preparation handles PSA
- [ ] Model compiles
- [ ] Short sampling runs without errors
- [ ] SLD mode backward compatible
- [ ] Full PSA run converges
