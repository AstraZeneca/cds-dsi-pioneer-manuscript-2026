# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Sclc is a Bayesian hierarchical modeling system for analyzing tumor dynamics and progression-free survival (PFS) in oncology trials. The codebase supports three projects:

1. **SCLC-01**: Longitudinal tumor analysis with state-space modeling (primary focus)
2. **Endometrial to LUNG**: Trial outcome predictions
3. **Breast-01 to -04**: Cross-validation predictions

## Build and Development Commands

### Initial Setup
```bash
# Restore R packages
Rscript -e 'renv::restore()'

# Install CmdStan
Rscript -e 'cmdstanr::install_cmdstan()'
```

### Running Analysis Pipelines
```bash
# Set project via TAR_PROJECT env var (see _targets.yaml for options)
export TAR_PROJECT=sclc
Rscript -e 'targets::tar_make()'
```

### Stan Model Syntax Check (fast)
```bash
~/.cmdstan/cmdstan-2.37.0/bin/stanc --include-paths=stan,stan/ssls stan/ssls/sf-ssm-log-space.stan
```

### Running Tests
```bash
Rscript -e 'testthat::test_dir("tests/testthat")'
```

## Architecture

### Stan Module System

Stan code uses modular `#include` architecture in `stan/ssls/modules/`:
- **tr/** - Tumor regression (decrease) dynamics
- **frac/** - Growth fraction dynamics
- **init/** - Initial state modeling
- **other_events/** - Non-target progression and death events
- **measurement/** - Observation model (measurement error)

Each module follows a 7-file pattern:
1. `flags.stan` - Feature switches
2. `data.stan` - Module-specific data
3. `hyperparams.stan` - Prior hyperparameters
4. `transformed_data.stan` - Data preprocessing
5. `parameters.stan` - Parameter declarations
6. `transformed_parameters.stan` - Derived quantities
7. `priors.stan` - Prior distributions

### Three-Level Hierarchical Parameters

Parameters follow a population → trial → patient hierarchy:
- Population level: `*_pop` (e.g., `tr_intercept_pop`)
- Trial level SD: `*_sd_trial_*` (e.g., `tr_sd_trial_intercept`)
- Trial level raw effects: `*_raw_trial_*` (non-centered parameterization)
- Patient level SD: `*_sd_patient_*`
- Patient level raw effects: `*_raw_patient_*`

### Module Prefixes
- `tr_*` - Tumor regression
- `frac_*` - Growth fraction
- `init_*` - Initial state
- `oe_*` - Other events

### Key Stan Models
- `stan/ssls/sf-ssm-log-space.stan` - Main state-space longitudinal survival model
- `stan/ssls/sf-ssls-lfo.stan` - Leave-future-out cross-validation variant

## Key Files

### R Code
- `r/priors.R` - Prior specifications (shared)
- `r/initializers.R` - Stan model initializers (shared)
- `r/sclc/prepare_analysis_data.R` - SCLC-specific data prep
- `targets/sclc_targets.R` - SCLC pipeline definition

### Configuration
- `_targets.yaml` - Workflow configurations for all projects
- `renv.lock` - Package versions
- `_quarto.yml` - Documentation settings

## Analysis Results Storage

Analysis results are stored in `/mnt/data/analysis-results/karim_naguib/sclc/<run_name>/`:
- **Targets store**: `/mnt/data/analysis-results/karim_naguib/sclc/<run_name>/_targets`
- **Fit CSVs**: `/mnt/data/analysis-results/karim_naguib/sclc/<run_name>/fit`

## Coding Guidelines

### Stan
- Use built-in zero constructors: `zeros_vector()`, `zeros_int_array()`
- Don't pass array/vector sizes as arguments; use `size()` internally
- Ignore linter warnings about code sections (modular `#include` architecture places code fragments across sections)
- Use non-centered parameterization (NCP) for hierarchical parameters

### R
- Use modern pipe operator `|>` (not `%>%`)
- Follow tidyverse style guide
- Prefer `purrr` and `dplyr` over base R loops
- Use `testthat` for unit tests

### Targets
- **NEVER use `tar_config_set(store = ...)`** - it changes global state and causes conflicts
- Always use explicit `store` argument: `tar_read(name, store = "path/_targets")`
- Same applies to all targets functions: `tar_meta()`, `tar_load()`, etc.

### Adding Module Parameters
1. Add feature flag in `modules/<module>/flags.stan`
2. Add hyperparameters in `modules/<module>/hyperparams.stan`
3. Add parameters in `modules/<module>/parameters.stan`
4. Implement transforms in `modules/<module>/transformed_parameters.stan`
5. Add priors in `modules/<module>/priors.stan`
6. Update `r/priors.R` with defaults
7. Update initializer in `r/initializers.R`
8. Update `targets/sclc_targets.R` with flag value

## Documentation

Key technical docs in `docs/`:
- `ARCHITECTURE.md` - System architecture and optimization details
- `ar1_process_noise.md` - AR(1) time-varying process noise implementation
- `OTHER_EVENTS_MODEL.md` - Other events model design
- `CHANGELOG.md` - Major changes and design decisions
