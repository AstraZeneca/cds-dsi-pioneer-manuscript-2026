# SF State-Space Refactoring - Complete

**Date**: 2026-02-10

## Summary

Separated general state-space model functions from tumor-specific implementations, recognizing that the Stein-Fojo (SF) framework is biomarker-agnostic and can be reused for PSA and other biomarkers.

## Key Insight

The SF model is a **general two-component log-space state-space model** with:
- Component 1: Decreasing dynamics (treatment effect)
- Component 2: Increasing dynamics (progression/resistance)

This applies to:
- Tumor burden (SLD) ✓ Current implementation
- PSA dynamics (planned)
- Any biomarker with competing decrease/growth dynamics

## What Changed

### File Structure

**Before:**
```
stan/
└── _sf_functions.stan (1926 lines)
    ├── General state-space dynamics
    └── Tumor-specific wrappers (mixed together)
```

**After:**
```
stan/
├── sf_state_space.stan (1921 lines) - GENERAL STATE-SPACE MODEL
│   ├── sf_log_space_transition
│   ├── sf_log_space_trajectory_ncp
│   ├── calc_states
│   ├── generate_patient_states_rng
│   └── ... (all core dynamics - biomarker agnostic)
│
└── modules/tumor/
    └── functions.stan (27 lines) - TUMOR-SPECIFIC WRAPPERS
        └── calc_log_sld_mean(states → SLD)
```

**Future (ready to add):**
```
modules/psa/
└── functions.stan - PSA-SPECIFIC WRAPPERS
    └── calc_log_psa_mean(states → PSA)
```

## Files Modified

1. **stan/_sf_functions.stan → stan/sf_state_space.stan**
   - Renamed to clarify it's a general state-space model
   - Added header documentation explaining biomarker-agnostic design
   - Removed `calc_log_sld_mean` (moved to tumor module)
   - 1926 → 1921 lines

2. **stan/modules/tumor/functions.stan** (NEW)
   - Created tumor-specific wrappers
   - Contains `calc_log_sld_mean` that converts states → SLD
   - 27 lines
   - Future PSA module will follow same pattern

3. **stan/sf-ssm-log-space.stan**
   - Changed `#include "_sf_functions.stan"` → `#include "sf_state_space.stan"`
   - Added `#include "modules/tumor/functions.stan"`

4. **stan/sf-ssls-lfo.stan**
   - Changed `#include "_sf_functions.stan"` → `#include "sf_state_space.stan"`
   - Added `#include "modules/tumor/functions.stan"`

5. **stan/sf-ssls-lfo-endpoints.stan**
   - Changed `#include "_sf_functions.stan"` → `#include "sf_state_space.stan"`
   - Added `#include "modules/tumor/functions.stan"`

## What Functions Are General?

### Core State-Space Dynamics (General - in sf_state_space.stan)
- `sf_log_space_transition` - State transitions between time points
- `sf_log_space_transition_ncp` - Non-centered transitions
- `sf_log_space_transition_lpdf` - Transition probability density
- `sf_log_space_trajectory_ncp` - Full trajectory generation
- `sf_log_space_trajectory_ncp_vectorized` - Vectorized computation
- `sf_log_space_obs_lpdf` - Observation probability density
- `calc_states` - Calculate state trajectories
- `calc_patient_states` - Patient-specific state calculations
- `calc_patient_states_rect` - For map_rect parallelization
- `calc_patient_process_noise` - Process noise calculation
- `generate_patient_states_rng` - Generate state trajectories
- `generate_patient_states_with_means_rng` - Generate with means
- `calculate_all_patients_endpoints_rng` - Simulate endpoints
- `aggregate_trial_metrics` - Aggregate by trial
- `aggregate_conditional_group_metrics` - Aggregate by groups
- `get_growth_lag_factor` - Growth lag calculations
- `multi_normal_rng` - Multivariate normal generation

### Biomarker-Specific (in modules/<biomarker>/functions.stan)
- **Tumor**: `calc_log_sld_mean(states, baseline) → log(SLD)`
- **PSA** (future): `calc_log_psa_mean(states, baseline) → log(PSA)`

The key is that states are just `[n_visits × 2]` matrices. What those states *mean* (SLD, PSA, etc.) is defined by the biomarker-specific wrapper functions.

## Design Pattern for Adding PSA

When adding PSA, follow this pattern:

```stan
// modules/psa/functions.stan
vector calc_log_psa_mean(matrix patient_states, real psa_baseline) {
  assert_equal(cols(patient_states), 2);

  // Same formula: total = exp(component1) + exp(component2)
  return to_vector(log_sum_exp(patient_states[, 1], patient_states[, 2])) + log(psa_baseline);
}
```

Then use the same general functions from `sf_state_space.stan`:
- `sf_log_space_trajectory_ncp()` to generate PSA states
- `calc_log_psa_mean()` to convert states → PSA values
- Same state-space infrastructure, different biomarker!

## Verification

All three models compile successfully:

```bash
✅ stan/sf-ssm-log-space.stan
✅ stan/sf-ssls-lfo.stan
✅ stan/sf-ssls-lfo-endpoints.stan
```

Verified with:
```bash
~/.cmdstan/cmdstan-2.38.0/bin/stanc --include-paths=stan stan/<model>.stan
```

## Benefits

1. **Code Reuse**: PSA can use the same 1900+ lines of state-space logic
2. **Clear Separation**: General dynamics vs. biomarker-specific conversions
3. **Modular Design**: Each biomarker module contains only its specific logic
4. **Maintainability**: Bug fixes to state-space model benefit all biomarkers
5. **Conceptual Clarity**: SF is a model structure, not a biomarker

## Why This Matters for PSA

When you add PSA:
- **No need to rewrite state-space dynamics** - use `sf_state_space.stan` as-is
- **Just write 20 lines** in `modules/psa/functions.stan` for PSA-specific conversions
- **Joint modeling ready**: Tumor and PSA can share the same underlying state-space framework

The SF model is like a "template" for modeling longitudinal biomarkers with treatment effect + resistance dynamics.

## Backward Compatibility Note

Some function parameters still use tumor-specific names (e.g., `sum_tumor_size_baseline`) for backward compatibility. When adding PSA, these could optionally be renamed to more general terms (e.g., `baseline_value`), but it's not necessary since the functions are biomarker-agnostic already.

## Git Commit

Files changed:
- Renamed: stan/_sf_functions.stan → stan/sf_state_space.stan
- Modified: stan/sf_state_space.stan (removed calc_log_sld_mean, added header)
- Created: stan/modules/tumor/functions.stan
- Modified: stan/sf-ssm-log-space.stan (updated includes)
- Modified: stan/sf-ssls-lfo.stan (updated includes)
- Modified: stan/sf-ssls-lfo-endpoints.stan (updated includes)
- Created: stan/SF_REFACTORING_COMPLETE.md (this file)
