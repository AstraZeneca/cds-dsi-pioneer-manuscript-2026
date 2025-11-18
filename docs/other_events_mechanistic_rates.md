# Other Events Model: Mechanistic Tumor Rate Covariates

## Overview

This document describes the addition of patient-specific mechanistic tumor dynamics rates as covariates in the other events (non-target progression, death) hazard model.

**Date Implemented:** November 13, 2025  
**Branch:** karim/non-target

## Motivation

Previously, the other events model only used current tumor burden (log(SLD)) as a time-varying tumor covariate. However, the mechanistic state-space model for tumor dynamics provides rich patient-specific information about:

1. **Decrease rate** (`patient_log_decrease_rate`): How fast the tumor regresses (treatment effect)
2. **Growth rate** (`patient_log_growth_rate`): How fast the tumor grows (disease aggressiveness)
3. **SLD velocity**: How fast the tumor size is changing over time (first derivative)

These rates and derivatives capture underlying tumor biology and treatment response that may be predictive of other events independent of current tumor size.

## Implementation Approach

The patient-level mechanistic rates (`patient_log_decrease_rate` and `patient_log_growth_rate`) are already computed in `_sf_transformed_parameters.stan` and are **used directly** in the other events module. No additional data passing or grid expansion is needed.

### Variables Used

The other events model accesses these variables directly:
- `patient_log_decrease_rate[n_patients]` - patient-specific log decrease rates
- `patient_log_growth_rate[n_patients]` - patient-specific log growth rates
- `states_full_grid[2][n_patients, max_t_width]` - time-varying tumor states for log(SLD)

All are already available in the transformed parameters block scope.

### Files Modified

#### Stan Files

1. **`stan/ssls/modules/other_events/transformed_data.stan`**
   - Sets `n_tumor_covar = 4` (hardcoded in transformed data block)
   - Defines the 4 tumor covariates to be used

2. **`stan/ssls/modules/other_events/transformed_parameters.stan`**
   - Uses `patient_log_decrease_rate[i]` and `patient_log_growth_rate[i]` directly
   - Computes SLD velocity as weekly change in log(SLD)
   - Z-score normalizes all covariates using population mean and SD
   - Broadcasts scalar rate effects to all time points (time-invariant)
   - Applies coefficients to compute hazard ratio contributions

3. **`stan/ssls/modules/other_events/data.stan`**
   - Removed `n_tumor_covar` from data block (now in transformed_data)
   - `tumor_sum_covar` is scaffolded but unused

4. **`stan/ssls/modules/other_events/parameters.stan`**
   - `oe_tumor_coef_pop` sized based on `n_tumor_covar` from transformed_data
   - Supports 4 tumor coefficients per cause

#### R Files

5. **`r/sclc/prepare_analysis_data.R`**
   - Set `n_tumor_covar = 4` to specify 4 tumor coefficients for priors
   - Added documentation that covariates are extracted from model parameters in Stan
   - `tumor_sum_covar` remains empty (scaffolded for potential future use)

6. **`r/priors.R`**
   - `oe_tumor_coef_pop_mean` and `oe_tumor_coef_pop_sd` use `stan_data$n_tumor_covar`
   - With `n_tumor_covar = 4`, provides priors for all 4 coefficients
   - Weakly informative: mean=0, sd=0.5 for all coefficients

## Tumor Covariate Coefficients

The other events model now estimates **4 population-level coefficients per cause**:

| Index | Covariate | Type | Interpretation | Expected Sign |
|-------|-----------|------|----------------|---------------|
| 1 | log(SLD) | Time-varying | Current tumor burden | Positive (+) |
| 2 | log(decrease rate) | Time-invariant | Treatment response rate | Negative (-) |
| 3 | log(growth rate) | Time-invariant | Disease aggressiveness | Positive (+) |
| 4 | SLD velocity | Time-varying | Rate of tumor change | Positive (+) |

### Z-Score Normalization

All four covariates are z-scored normalized:
- **log(SLD)**: Normalized using mean/SD of ALL observed SLD values across all patients and visits
- **log(decrease rate)**: Normalized using mean/SD of patient-level rates
- **log(growth rate)**: Normalized using mean/SD of patient-level rates
- **SLD velocity**: Normalized using mean/SD of velocity values (excluding first time point)

This makes coefficients interpretable as **log hazard ratio per 1-SD change** in the covariate.

## Prior Specifications

Weakly informative priors on all four coefficients (per cause):

```r
# With n_tumor_covar = 4, this creates 4 coefficients per cause:
oe_tumor_coef_pop_mean = array(rep(0, 4), dim = c(n_causes, 4))
oe_tumor_coef_pop_sd   = array(rep(0.5, 4), dim = c(n_causes, 4))
```

Expected effects based on clinical reasoning:
- **log(SLD)**: Positive (+) - Higher tumor burden increases risk
- **log(decrease rate)**: Negative (-) - Faster regression suggests good response, lower risk
- **log(growth rate)**: Positive (+) - Faster growth suggests aggressive disease, higher risk
- **SLD velocity**: Positive (+) - Increasing tumor size suggests worsening disease, higher risk

These allow the data to dominate while providing gentle regularization.

## Time-Varying vs Time-Invariant

- **log(SLD)**: Time-varying (changes at each week based on model predictions)
- **log(decrease rate)**: Time-invariant (constant patient-specific value, broadcasted to all time points)
- **log(growth rate)**: Time-invariant (constant patient-specific value, broadcasted to all time points)
- **SLD velocity**: Time-varying (computed as weekly change in log(SLD))

The rates are intrinsic patient-level characteristics estimated from the tumor dynamics model, while SLD and velocity evolve over time.

## Backward Compatibility

The implementation maintains backward compatibility:

1. If `oe_tumor_coef_pop` has size 1 (only log(SLD) coefficient), the model runs as before
2. If `oe_tumor_coef_pop` has size >= 3, the rate effects are included
3. If `oe_tumor_coef_pop` has size >= 4, the velocity effect is also included
4. The `n_tumor_covar` is set in both R (for priors) and Stan transformed_data (hardcoded)

## Usage in Targets Workflow

No changes needed to `targets/sclc_targets.R` - the new covariates are automatically available from the mechanistic model.

To use only log(SLD), set in R:
```r
n_tumor_covar <- 1  # in prepare_tumor_stan_data
```

To use SLD + rates (no velocity), set:
```r
n_tumor_covar <- 3
```

To use all covariates including velocity (default):
```r
n_tumor_covar <- 4
```

## Testing Recommendations

1. **Fit with only log(SLD)** (size 1) to confirm backward compatibility
2. **Fit with all 3 covariates** to assess incremental value of mechanistic rates
3. **Compare via LOO-CV** to quantify predictive improvement
4. **Check coefficient signs** match prior expectations
5. **Examine patient-level predictions** for clinical plausibility

## Future Extensions

Potential additional tumor covariates that could be added to `states_full_grid`:

- SLD velocity (first derivative): `d(log_SLD)/dt`
- SLD acceleration (second derivative): `d²(log_SLD)/dt²`
- Growth fraction: proportion in growing state
- Resistance score: increase in growth rate after nadir
- Time since nadir SLD

## References

- **Naming conventions:** `docs/multi_level_hierarchy_design.md`
- **Other events design:** `docs/other_events_covariate_plan.md` (if it exists)
- **State space model:** `sld_state_space_model.md`

---

**Maintainer:** Update this document when modifying tumor covariate structure
**Last Updated:** November 13, 2025
