# State-Space Module Organization

**Date**: 2026-02-10

## Summary

Organized all Stein-Fojo (SF) state-space infrastructure into a dedicated `modules/state_space/` module, removing underscore prefixes and grouping related files together.

## Motivation

The SF state-space files were scattered at the root level with underscore prefixes (`_sf_*`), making them hard to identify as a cohesive unit. These files represent the **state-space modeling infrastructure** that orchestrates the parameter modules (tr, frac, init, etc.).

## What is the State-Space Module?

The state-space module contains the infrastructure for:
- Two-component log-space dynamics (decrease + growth)
- State trajectory computation
- Outcome predictions (PFS, ORR, etc.)
- Model validation and accuracy metrics

It's **distinct from parameter modules** (tr/frac/init):
- Parameter modules: Define specific model parameters (tumor regression, growth fraction, etc.)
- State-space module: Orchestrates those parameters into trajectories and outcomes

## File Reorganization

| Before | After | Purpose |
|--------|-------|---------|
| `sf_state_space.stanfunctions` | `modules/state_space/functions.stanfunctions` | Core SF dynamics functions |
| `_sf_outcomes_info.stan` | `modules/state_space/data.stan` | Outcome configuration (quantiles, timepoints) |
| `_sf_transformed_data.stan` | `modules/state_space/transformed_data.stan` | Setup and preprocessing |
| `_sf_transformed_parameters.stan` | `modules/state_space/transformed_parameters.stan` | State trajectory computation |
| `_sf_accuracy_generated_quantities.stan` | `modules/state_space/generated_quantities.stan` | RECIST accuracy metrics |
| `_sf-checks.stan` | `modules/state_space/checks.stan` | Validation checks |
| `_sf-ssls-lfo-data.stan` | `modules/state_space/lfo_data.stan` | LFO-specific data |

## New Structure

```
stan/modules/
├── state_space/ - STATE-SPACE INFRASTRUCTURE MODULE
│   ├── functions.stanfunctions - Core SF dynamics (1900+ lines)
│   │   ├── sf_log_space_transition
│   │   ├── sf_log_space_trajectory_ncp
│   │   ├── calc_states, calc_patient_states
│   │   ├── generate_patient_states_rng
│   │   └── calculate_all_patients_endpoints_rng
│   │
│   ├── data.stan - Outcome configuration
│   │   ├── PFS quantiles, timepoints
│   │   └── Conditioning groups
│   │
│   ├── lfo_data.stan - LFO-specific configuration
│   │
│   ├── transformed_data.stan - Setup
│   │   ├── Normalized SLD
│   │   ├── Visit indexing
│   │   ├── RECIST constants
│   │   ├── QR decomposition of covariates
│   │   └── Visit cumsum matrix
│   │
│   ├── transformed_parameters.stan - State computation
│   │   ├── Combine rates from tr/frac modules
│   │   ├── Calculate patient states
│   │   └── Process noise handling
│   │
│   ├── checks.stan - Validation
│   │   └── Assert state dimensions
│   │
│   └── generated_quantities.stan - Outcomes
│       └── RECIST prediction accuracy
│
└── modules/ - PARAMETER MODULES
    ├── tr/ - Tumor regression parameters
    ├── frac/ - Growth fraction parameters
    ├── init/ - Initial state parameters
    ├── measurement/ - Measurement model
    ├── other_events/ - Non-tumor progression
    └── tumor/ - Tumor-specific functions
```

## Updated Includes

### Before:
```stan
functions {
  #include "sf_state_space.stanfunctions"
  ...
}
data {
  #include "_sf_outcomes_info.stan"
  #include "_sf-ssls-lfo-data.stan"  // LFO models only
  ...
}
transformed data {
  #include "_sf_transformed_data.stan"
  #include "_sf-checks.stan"
  ...
}
transformed parameters {
  #include "_sf_transformed_parameters.stan"
  ...
}
generated quantities {
  #include "_sf_accuracy_generated_quantities.stan"
  ...
}
```

### After:
```stan
functions {
  #include "modules/state_space/functions.stanfunctions"
  ...
}
data {
  #include "modules/state_space/data.stan"
  #include "modules/state_space/lfo_data.stan"  // LFO models only
  ...
}
transformed data {
  #include "modules/state_space/transformed_data.stan"
  #include "modules/state_space/checks.stan"
  ...
}
transformed parameters {
  #include "modules/state_space/transformed_parameters.stan"
  ...
}
generated quantities {
  #include "modules/state_space/generated_quantities.stan"  // main model only
  ...
}
```

## Benefits

1. **Clear Module Identity**: `modules/state_space/` directory makes it obvious these files are related
2. **No Underscore Prefixes**: Cleaner filenames without `_sf_` prefix
3. **Logical Grouping**: All state-space infrastructure in one place
4. **Separation of Concerns**: Clear distinction between state-space infrastructure and parameter modules
5. **Easier Navigation**: Single directory to find all state-space code

## Conceptual Architecture

```
┌─────────────────────────────────────────────────────┐
│                 STATE-SPACE MODULE                  │
│          (Orchestration & Computation)              │
│                                                     │
│  - Combines parameters from modules                │
│  - Computes state trajectories                     │
│  - Generates predictions & outcomes                │
└─────────────────────────────────────────────────────┘
                         ▲
                         │ Uses
                         │
    ┌────────────────────┼────────────────────┐
    │                    │                    │
┌───┴───┐          ┌─────┴─────┐       ┌─────┴─────┐
│  tr/  │          │   frac/   │       │   init/   │
│module │          │  module   │       │  module   │
└───────┘          └───────────┘       └───────────┘
   └─── Tumor regression, growth fraction, initial state parameters ───┘
```

## State-Space vs Parameter Modules

| Aspect | State-Space Module | Parameter Modules |
|--------|-------------------|-------------------|
| **Purpose** | Orchestrate computation | Define parameters |
| **Scope** | Full model infrastructure | Specific components |
| **Dependencies** | Uses all parameter modules | Independent |
| **Location** | `modules/state_space/` | `modules/<name>/` |
| **Examples** | Trajectory computation, outcomes | tr, frac, init |

## Verification

All three models compile successfully:
```bash
✅ stan/sf-ssm-log-space.stan
✅ stan/sf-ssls-lfo.stan
✅ stan/sf-ssls-lfo-endpoints.stan
```

## Future: PSA Addition

When adding PSA, the state-space module remains unchanged:
- PSA will use the same `modules/state_space/` infrastructure
- Only need to add PSA-specific functions in `modules/psa/`
- State-space module orchestrates both tumor and PSA states

The separation makes this clear: state-space is the framework, biomarker modules are the content.
