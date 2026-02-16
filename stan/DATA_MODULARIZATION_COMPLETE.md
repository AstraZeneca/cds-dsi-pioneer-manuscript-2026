# Data Modularization - Complete

**Date**: 2026-02-10

## Summary

Successfully separated tumor-specific data from base patient/visit structure to enable modular addition of PSA and other biomarkers.

## What Changed

### File Structure

**Before:**
```
stan/
├── base_data.stan (99 lines)
│   ├── Base structure (patients, visits, hierarchy)
│   └── Tumor data mixed in (sum_tumor_size, recist, pfs, etc.)
```

**After:**
```
stan/
├── base_data.stan (62 lines) - BASE ONLY
│   ├── Trials, patients, hierarchy
│   ├── Visit schedules (shared by all biomarkers)
│   └── Patient covariates
└── modules/
    └── tumor/ - TUMOR MODULE (following standard module pattern)
        ├── data.stan (57 lines)
        │   ├── sum_tumor_size (SLD measurements)
        │   ├── recist (response categories)
        │   ├── PFS outcomes (pfs, right_censored, interval_censored)
        │   ├── Target-lesion PFS (target_pfs, target_right_censored)
        │   ├── death_week
        │   └── Model config (fit_tumor_data, sf_rep_T, debug, n_shards)
        └── transformed_data.stan (63 lines)
            ├── log_sum_tumor_size, post_treat_sld
            ├── Population/patient visit indices
            └── Forecast visit calculations
```

### Key Separation Principle

**Shared Visit Schedule**: The visit schedule (`n_patient_visits`, `t_patient_visits`) is defined in `base_data.stan` and SHARED by all biomarkers:

```
base_data.stan:
  array[n_patients] int<lower = 1> n_patient_visits;
  array[sum(n_patient_visits)] int t_patient_visits;

modules/tumor/data.stan:
  vector<lower = 0>[sum(n_patient_visits)] sum_tumor_size;  // ← Same length

modules/psa/data.stan (future):
  vector<lower = 0>[sum(n_patient_visits)] psa_values;     // ← Same length
```

This enables:
- Joint modeling of multiple biomarkers
- Easy indexing: visit `v` has `tumor_size[v]`, `psa[v]`, etc.
- Missing data handling via indicators

## Files Modified

### Stan Models

1. **stan/base_data.stan**
   - Removed tumor-specific section (lines 75-99)
   - Removed `sum_tumor_size` declaration (line 51)
   - Updated comments to reference modular structure
   - Now 62 lines (was 99)

2. **stan/modules/tumor/data.stan** (NEW)
   - Created modular tumor data file (follows standard module pattern)
   - Well-documented sections: measurements, outcomes, config
   - 57 lines

3. **stan/modules/tumor/transformed_data.stan** (MOVED)
   - Moved from `stan/tumor/tumor_transformed_data.stan`
   - Renamed to follow standard module naming convention
   - No content changes
   - 63 lines

4. **stan/sf-ssm-log-space.stan**
   - Added `#include "modules/tumor/data.stan"` after base_data
   - Changed `#include "tumor/tumor_transformed_data.stan"` to `#include "modules/tumor/transformed_data.stan"`

5. **stan/sf-ssls-lfo.stan**
   - Added `#include "modules/tumor/data.stan"` after base_data
   - Changed `#include "tumor/tumor_transformed_data.stan"` to `#include "modules/tumor/transformed_data.stan"`

6. **stan/sf-ssls-lfo-endpoints.stan**
   - Added `#include "modules/tumor/data.stan"` after base_data
   - Changed `#include "tumor/tumor_transformed_data.stan"` to `#include "modules/tumor/transformed_data.stan"`

### Documentation

6. **stan/DATA_REFACTORING_PLAN.md** (NEW)
   - Detailed refactoring plan with rationale
   - Template for future PSA addition

7. **stan/DATA_MODULARIZATION_COMPLETE.md** (this file)

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

1. **Cleaner Separation**: Base structure vs. biomarker measurements clearly separated
2. **Modular Addition**: Adding PSA just requires:
   ```stan
   data {
     #include "base_data.stan"
     #include "modules/tumor/data.stan"
     #include "modules/psa/data.stan"      // ← Just add this
     // ...
   }
   ```
3. **Shared Infrastructure**: All biomarkers use common visit schedule
4. **Joint Modeling Ready**: Multiple longitudinal biomarkers naturally aligned
5. **Backward Compatible**: No changes to data structure, just organization

## Data Flow

```
R data preparation
    ↓
Creates lists:
    - base_data: list(n_patients, n_patient_visits, t_patient_visits, ...)
    - tumor_data: list(sum_tumor_size, recist, pfs, ...)
    - psa_data: list(psa_values, psa_pfs, ...)  [future]
    ↓
Stan model reads:
    #include "base_data.stan"      → base_data list
    #include "tumor/data.stan"     → tumor_data list
    #include "psa/data.stan"       → psa_data list [future]
```

## R Code Impact

The R data preparation functions will need to separate data into:
- `base_data` - Core structure
- `tumor_data` - Tumor measurements

However, if passing as a combined list, Stan will still work as all variables are declared.

## Next Steps for PSA Addition

1. Create `stan/modules/psa/data.stan` using modules/tumor/data.stan as template
2. Add `#include "modules/psa/data.stan"` to model files
3. Create `stan/modules/psa/transformed_data.stan` for PSA-specific processing
4. Add PSA state-space dynamics in parameters/transformed_parameters
5. Add PSA measurement model

The modular structure is now ready for this!

## Module Structure Consistency

The tumor module now follows the same pattern as other modules:

```
stan/modules/
├── tr/ - Tumor regression parameters
├── frac/ - Growth fraction parameters
├── init/ - Initial state parameters
├── measurement/ - Measurement model parameters
├── other_events/ - Non-tumor progression events
│   └── data.stan (has module-specific data)
└── tumor/ - Tumor measurement data (NEW)
    ├── data.stan (tumor-specific data inputs)
    └── transformed_data.stan (tumor-specific computations)
```

Future PSA module will follow the same pattern:
```
stan/modules/psa/
├── data.stan
├── transformed_data.stan
└── [other module files as needed]
```

## Git Commit

Files changed:
- Modified: stan/base_data.stan
- Modified: stan/sf-ssm-log-space.stan
- Modified: stan/sf-ssls-lfo.stan
- Modified: stan/sf-ssls-lfo-endpoints.stan
- Created: stan/tumor/data.stan
- Created: stan/DATA_REFACTORING_PLAN.md
- Created: stan/DATA_MODULARIZATION_COMPLETE.md

Lines changed:
- base_data.stan: 99 → 62 lines (-37 lines)
- tumor/data.stan: 0 → 57 lines (+57 lines)
- Net change: +20 lines (due to documentation)
