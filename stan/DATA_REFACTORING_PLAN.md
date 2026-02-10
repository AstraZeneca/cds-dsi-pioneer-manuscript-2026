# Data Modularization Plan

## Goal
Separate tumor-specific longitudinal data from base patient/visit structure to enable easy addition of PSA and other biomarkers.

## Current Structure
```
base_data.stan (99 lines)
├── Lines 1-74: Base patient/visit structure
└── Lines 75-99: TUMOR-SPECIFIC DATA section (already marked!)
    ├── sum_tumor_size (SLD measurements)
    ├── recist (response categories)
    ├── pfs, right_censored, interval_censored
    ├── target_pfs, target_right_censored
    ├── death_week
    ├── fit_tumor_data, sf_rep_T, debug, n_shards
    └── n_covar, covar_design_matrix
```

## Proposed Structure

### base_data.stan (BASE ONLY)
**Core patient/visit structure - shared by all biomarkers**

```stan
// Trial and Patient Hierarchy
int<lower = 1> n_trials;
int<lower = 0> n_patients;
array[n_patients] int<lower = 1, upper = n_trials> patient_trial;

// Multi-level hierarchy (trial → region → site → patient, etc.)
int<lower = 1> n_levels;
array[n_levels] int<lower = 1> n_groups_per_level;
array[n_patients, n_levels] int<lower = 1> patient_level_groups;

// Visit schedule (SHARED across all biomarkers)
array[n_patients] int<lower = 1> n_patient_visits;
array[sum(n_patient_visits)] int t_patient_visits;       // Weeks
array[sum(n_patient_visits)] int t_patient_visits_day;   // Days

// Calendar information
array[n_patients] int<lower = 1> calendar_week;
array[n_patients] int<lower = 1> calendar_day;

// Time horizon extension
int<lower = 1> extend_max_all_t;

// Patient-level covariates (SHARED)
int<lower = 0> n_covar;
matrix[n_patients, n_covar] covar_design_matrix;
```

### modules/tumor/data.stan (NEW FILE)
**Tumor/SLD-specific measurements and outcomes**

```stan
// Tumor measurement data (aligned with visit schedule from base_data)
vector<lower = 0>[sum(n_patient_visits)] sum_tumor_size; // SLD in cm
array[sum(n_patient_visits)] int<lower = 1, upper = 5> recist;

// Tumor-based progression outcomes
array[n_patients] int<lower = 0> pfs;  // Progression-free survival (weeks)
array[n_patients] int<lower = 0, upper = 1> right_censored;
array[n_patients] int<lower = 0> interval_censored;

// Target-lesion-only progression
array[n_patients] int<lower = 0> target_pfs;
array[n_patients] int<lower = 0, upper = 1> target_right_censored;

// Death information
array[n_patients] int<lower = 0> death_week;

// Model configuration
int<lower = 0, upper = 1> fit_tumor_data;
int<lower = 0> sf_rep_T;
int<lower = 0, upper = 1> debug;
int<lower = 1, upper = n_patients> n_shards;
```

### Future: modules/psa/data.stan (TEMPLATE)
**PSA-specific measurements and outcomes**

```stan
// PSA measurement data (aligned with visit schedule from base_data)
vector<lower = 0>[sum(n_patient_visits)] psa_values;    // ng/mL
array[sum(n_patient_visits)] int<lower = 0, upper = 1> psa_measured;  // Indicator

// PSA-based progression outcomes
array[n_patients] int<lower = 0> psa_pfs;  // PSA progression (weeks)
array[n_patients] int<lower = 0, upper = 1> psa_right_censored;

// Model configuration
int<lower = 0, upper = 1> fit_psa_data;
```

## Key Design Principles

1. **Shared Visit Schedule**: `n_patient_visits`, `t_patient_visits` are in base_data and shared by ALL biomarkers
   - Tumor SLD measured at these visits
   - PSA measured at these visits (with indicator for missing)
   - Allows joint modeling of multiple biomarkers

2. **Ragged Array Alignment**: All longitudinal measurements use `sum(n_patient_visits)` length
   - Easy to index: visit `v` has tumor_size[v], psa[v], etc.
   - Missing data handled via indicators (e.g., `psa_measured[v]`)

3. **Separate Progression Outcomes**: Each biomarker can have its own PFS
   - `pfs` - tumor-based progression (RECIST)
   - `psa_pfs` - PSA-based progression (future)
   - Composite outcome can be computed in transformed data

4. **Modular Fitting**: Each biomarker has a `fit_*_data` flag
   - Enables tumor-only, PSA-only, or joint models
   - Easy to turn on/off components

## Files to Modify

1. **stan/base_data.stan** - Remove tumor section (lines 75-99)
2. **stan/modules/tumor/data.stan** - Create new file with tumor data
3. **stan/modules/tumor/transformed_data.stan** - Move from stan/tumor/tumor_transformed_data.stan
4. **stan/sf-ssm-log-space.stan** - Add `#include "modules/tumor/data.stan"` after base_data
5. **stan/sf-ssls-lfo.stan** - Add `#include "modules/tumor/data.stan"` after base_data
6. **stan/sf-ssls-lfo-endpoints.stan** - Add `#include "modules/tumor/data.stan"` after base_data
7. **R data preparation** - Update to separate base vs tumor data lists (optional)

## Migration Strategy

1. Create `tumor/data.stan` with tumor-specific variables
2. Update model files to include tumor/data.stan
3. Verify compilation
4. Remove tumor section from base_data.stan
5. Verify compilation again
6. Update R data prep functions (if needed)

## Benefits

- **Cleaner separation**: Base structure vs. biomarker measurements
- **Easy PSA addition**: Just add `psa/data.stan` and `#include` it
- **Joint modeling**: Multiple biomarkers share visit schedule
- **Backward compatible**: No changes to data structure, just organization
