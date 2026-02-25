# _sf_functions.stan Analysis

## Overview
- **Lines**: 1926
- **Functions**: 29 definitions (18 unique names with overloads)
- **Purpose**: Stein-Fojo (SF) tumor growth state-space model functions

## Function Categories

### 1. State-Space Model Core (6 function groups)
**Purpose**: Core dynamics of tumor growth/regression in log-space

- `sf_log_space_transition` (2 overloads) - State transitions
- `sf_log_space_transition_ncp` (2 overloads) - Non-centered parameterization transitions
- `sf_log_space_transition_lpdf` - Transition log probability density
- `sf_log_space_trajectory_ncp` (4 overloads) - Full trajectory generation
- `sf_log_space_trajectory_ncp_vectorized` - Vectorized trajectory computation
- `sf_log_space_obs_lpdf` - Observation log probability density

### 2. State Calculation (5 function groups)
**Purpose**: Calculate patient-specific tumor trajectories

- `calc_states` (3 overloads) - Calculate state trajectories
- `calc_patient_states` (3 overloads) - Patient-specific state calculations
- `calc_patient_states_rect` - For map_rect parallelization
- `calc_patient_process_noise` - Process noise calculation
- `calc_log_sld_mean` - Calculate log SLD means from states

### 3. Outcome Simulation (4 function groups)
**Purpose**: Generate PFS, ORR, and other clinical endpoints

- `generate_patient_states_rng` - Generate patient state trajectories (RNG)
- `generate_patient_states_with_means_rng` - Generate with expected means
- `calculate_all_patients_endpoints_rng` - Calculate PFS, ORR endpoints for all patients
- `aggregate_trial_metrics` - Aggregate outcomes by trial
- `aggregate_conditional_group_metrics` - Aggregate by conditioning groups

### 4. Helper Functions (2 function groups)
**Purpose**: Utility functions for tumor model

- `get_growth_lag_factor` (2 overloads) - Growth lag factor calculation
- `multi_normal_rng` (2 overloads) - Multivariate normal random generation

## Current Issue

The `_sf_functions.stan` file:
1. Has underscore prefix (suggests temporary/internal status)
2. Contains 100% tumor-specific code (nothing general-purpose)
3. Doesn't fit in general utility files like `util.stan`, `pos.stan`, `gp.stan`
4. Is separate from the newly created `modules/tumor/` structure

## Recommendation

**Move to `modules/tumor/functions.stan`**

### Rationale:
1. **All functions are tumor-specific** - They're all about Stein-Fojo tumor dynamics
2. **Consistent with module structure** - `modules/tumor/` is now the home for tumor code
3. **Not general utilities** - These aren't general-purpose like `pos.stan` (position arrays) or `util.stan` (sorting, filtering)
4. **Matches other specialized files** - Similar to how `recist.stanfunctions` is RECIST-specific

### Proposed Structure:
```
modules/tumor/
├── data.stan - Tumor measurement inputs
├── transformed_data.stan - Tumor data processing
└── functions.stan - Tumor state-space model functions (NEW, from _sf_functions.stan)
```

### Alternative Consideration:

If the functions need to be more granular, could split into:
```
modules/tumor/
├── functions.stan - Core SF state-space functions
├── state_functions.stan - State calculation functions
└── outcome_functions.stan - PFS/ORR simulation functions
```

But this seems like over-engineering. A single `functions.stan` keeps it clean.

## Why NOT Move to Existing Utility Files?

### util.stan
- Contains: General sorting, filtering, validation
- **Bad fit**: Tumor state-space dynamics are too specialized

### gp.stan
- Contains: Gaussian Process functions (kernels, covariance matrices)
- **Bad fit**: These are tumor growth model functions, not GP functions

### pfs_functions.stan
- Contains: PFS outcome calculations (Kaplan-Meier, hazard functions)
- **Possible but messy**: Would mix state-space dynamics with outcome metrics

### pos.stan
- Contains: Position array indexing utilities
- **Bad fit**: No relation to position arrays

## Implementation Plan

1. Rename `stan/_sf_functions.stan` → `stan/modules/tumor/functions.stan`
2. Update model includes:
   - Change `#include "_sf_functions.stan"` → `#include "modules/tumor/functions.stan"`
3. Update all three model files:
   - sf-ssm-log-space.stan
   - sf-ssls-lfo.stan
   - sf-ssls-lfo-endpoints.stan
4. Verify compilation

## Files to Modify

- stan/_sf_functions.stan (move to modules/tumor/functions.stan)
- stan/sf-ssm-log-space.stan (update include)
- stan/sf-ssls-lfo.stan (update include)
- stan/sf-ssls-lfo-endpoints.stan (update include)
