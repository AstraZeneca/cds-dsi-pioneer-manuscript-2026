# Sclc System Architecture

**Last Updated:** November 14, 2025

## Table of Contents

1. [Overview](#overview)
2. [Multi-Level Hierarchical Parameter Design](#multi-level-hierarchical-parameter-design)
3. [Stan Optimization for State-Space Models](#stan-optimization-for-state-space-models)
4. [Module Structure](#module-structure)

---

## Overview

This document describes the architectural design patterns and computational optimization strategies used in the Sclc Bayesian hierarchical modeling system for oncology trial analysis.

Sclc uses a modular Stan architecture with:
- **Three-level hierarchical parameters** (population → trial → patient)
- **Optimized state-space computations** for tumor dynamics
- **Flexible feature flags** for incremental model complexity
- **QR decomposition** for numerical stability

---

## Multi-Level Hierarchical Parameter Design

### Design Philosophy

The Sclc model uses a **three-level hierarchy** for most parameters:

1. **Population level** - Shared across all trials and patients
2. **Trial level** - Shared within a trial, varying between trials
3. **Patient level** - Individual patient effects

This design allows:
- **Information sharing** across trials while respecting trial-specific effects
- **Shrinkage toward population means** when data is sparse
- **Explicit modeling** of between-trial and within-trial variation
- **Interpretable hyperparameters** that control shrinkage strength

### Hierarchical Structure Pattern

#### Standard Three-Level Hierarchy

For a parameter like tumor regression rate intercept:

```stan
// Population level (global mean)
real tr_intercept_pop;

// Trial level (trial-specific deviations)
real<lower=0> tr_sd_trial_intercept;           // Between-trial SD
vector[n_trials] tr_raw_trial_intercept;       // Raw trial effects (std normal)
vector[n_trials] tr_trial_intercept;           // Centered trial effects

// Patient level (patient-specific deviations)
real<lower=0> tr_sd_patient_intercept;         // Within-trial SD
vector[n_patients] tr_raw_patient_intercept;   // Raw patient effects (std normal)
vector[n_patients] tr_patient_intercept;       // Final patient-level parameters
```

#### Centering Computation

```stan
// Trial effects centered on population mean
tr_trial_intercept = tr_intercept_pop + tr_sd_trial_intercept * tr_raw_trial_intercept;

// Patient effects centered on their trial mean
tr_patient_intercept = tr_trial_intercept[patient_trial] + 
                       tr_sd_patient_intercept * tr_raw_patient_intercept;
```

### Module-Specific Naming Conventions

All parameters follow consistent prefixes based on their module:

| Module | Prefix | Example |
|--------|--------|---------|
| Tumor Regression | `tr_` | `tr_intercept_pop`, `tr_sd_trial_intercept` |
| Growth Fraction | `frac_` | `frac_intercept_pop`, `frac_sd_patient_intercept` |
| Initial State | `init_` | `init_intercept_pop`, `init_sd_trial_intercept` |
| Other Events | `oe_` | `oe_baseline_hazard_pop`, `oe_sd_trial_baseline_hazard` |

### Parameter Naming Template

```
<module>_<quantity>_<level>_<detail>
```

**Examples:**
- `tr_intercept_pop` - Tumor regression intercept at population level
- `tr_sd_trial_intercept` - Standard deviation of trial-level tumor regression intercepts
- `tr_raw_trial_intercept` - Raw (non-centered) trial-level effects
- `frac_sd_patient_intercept` - Standard deviation of patient-level growth fraction intercepts

**For comprehensive naming details, see:** [`docs/multi_level_hierarchy_design.md`](multi_level_hierarchy_design.md) (archived reference)

### Covariate Effects with QR Decomposition

For numerical stability, covariate effects use QR decomposition:

```stan
// Population-level coefficients (in QR space)
vector[n_covars] tr_coef_qr_pop;

// Trial-level coefficient SDs (in QR space)
vector<lower=0>[n_covars] tr_sd_trial_coef_qr;

// Trial-level raw effects (in QR space)
matrix[n_trials, n_covars] tr_raw_trial_coef_qr;

// Transform back to original covariate space
matrix[n_trials, n_covars] tr_trial_coef = 
  tr_raw_trial_coef_qr * diag_matrix(tr_sd_trial_coef_qr);
```

### Feature Flags for Optional Parameters

Each module uses feature flags to enable/disable hierarchical levels:

```stan
// flags.stan
int<lower=0, upper=1> enable_trial_intercept_tr;
int<lower=0, upper=1> enable_trial_coef_tr;
```

This allows:
- **Baseline models** with all flags disabled (population-level only)
- **Incremental complexity** by enabling features one at a time
- **Model comparison** to assess value of additional hierarchy

### Hyperparameter Specifications

Hyperparameters control the strength of shrinkage:

```stan
// hyperparams.stan
real tr_mean_intercept_pop;                   // Prior mean for population intercept
real<lower=0> tr_sd_intercept_pop;            // Prior SD for population intercept
real<lower=0> tr_sd_trial_intercept_hp;       // Prior scale for between-trial SD
real<lower=0> tr_sd_patient_intercept_hp;     // Prior scale for within-trial SD
```

**Naming convention:**
- Suffix `_hp` indicates a hyperparameter (prior specification)
- Hyperparameters are passed as data, not estimated

### Priors

Standard priors follow this pattern:

```stan
// priors.stan

// Population level
tr_intercept_pop ~ normal(tr_mean_intercept_pop, tr_sd_intercept_pop);

// Trial level (when enabled)
if (enable_trial_intercept_tr) {
    tr_sd_trial_intercept ~ normal(0, tr_sd_trial_intercept_hp);
    tr_raw_trial_intercept ~ std_normal();
}

// Patient level
tr_sd_patient_intercept ~ normal(0, tr_sd_patient_intercept_hp);
tr_raw_patient_intercept ~ std_normal();
```

---

## Stan Optimization for State-Space Models

### Performance Challenge

State-space models typically have sequential dependencies where state $t$ depends on state $t-1$, preventing vectorization. The Sclc tumor dynamics model faced this bottleneck.

### Mathematical Foundation

The Stein-Fojo tumor growth model has linear state transitions:

$$
\begin{aligned}
x_d[t] &= x_d[t-1] - d \cdot \Delta t + \epsilon_{d,t} \\
x_g[t] &= x_g[t-1] + g \cdot \Delta t + \epsilon_{g,t}
\end{aligned}
$$

This can be rewritten as a cumulative sum:

$$
x[t] = x[0] + \sum_{k=1}^{t} \delta[k]
$$

### Optimization Strategy: Vectorized Cumulative Sum

Replace sequential loop with Stan's `cumulative_sum()`:

```stan
// BEFORE (slow)
for (t in 2:T) {
  states[t, 1] = states[t-1, 1] + increment_d[t];
  states[t, 2] = states[t-1, 2] + increment_g[t];
}

// AFTER (fast)
vector[T-1] d_cumsum = cumulative_sum(d_increments);
vector[T-1] g_cumsum = cumulative_sum(g_increments);

for (t in 2:T) {
  states[t, 1] = initial_state_d + d_cumsum[t-1];
  states[t, 2] = initial_state_g + g_cumsum[t-1];
}
```

**Benefits:**
- **2-10x speedup** from eliminating sequential dependencies
- Better CPU cache utilization
- Enables SIMD vectorization

### GPU-Ready Matrix Formulation

For future GPU acceleration, a matrix-based approach eliminates `map_rect`:

**Step 1:** Create cumulative sum indicator matrix (computed once in transformed data):

```stan
// Lower triangular matrix of 1's for cumulative sums
matrix[n_unique_visits, time_range] visit_cumsum_mat;
```

**Step 2:** Batched computation across all patients:

```stan
// All patients at once
matrix[n_unique_visits, n_patients] R_d = visit_cumsum_mat * D_increments';
matrix[n_unique_visits, n_patients] R_g = visit_cumsum_mat * G_increments';
```

**Expected speedup:** 10-40x with GPU hardware

**Trade-offs:**
- Computes states at all time points (not just visits)
- ~2% utilization in sparse visit scenarios
- But GPU matrix multiply is so fast that this is still faster overall

For complete optimization details, see archived [`docs/stan_state_space_optimization.md`](stan_state_space_optimization.md)

---

## Module Structure

### Stan Module System

The Stan code uses a modular architecture with `#include` directives. Modules are organized by feature in `stan/ssls/modules/`:

- **`tr/`** - Tumor regression (decrease) dynamics
- **`frac/`** - Growth fraction dynamics  
- **`init/`** - Initial state modeling
- **`other_events/`** - Non-target progression and death events
  - Uses mechanistic tumor rates as covariates (see [`docs/other_events_mechanistic_rates.md`](other_events_mechanistic_rates.md))

### Module File Pattern

Each module follows a 7-file structure:

1. **`flags.stan`** - Feature switches (e.g., `enable_trial_intercept_tr`)
2. **`data.stan`** - Data declarations specific to module
3. **`hyperparams.stan`** - Prior hyperparameters
4. **`transformed_data.stan`** - Data preprocessing
5. **`parameters.stan`** - Parameter declarations
6. **`transformed_parameters.stan`** - Derived quantities
7. **`priors.stan`** - Prior distributions

This standardization makes modules:
- **Self-contained** - Each module includes all its dependencies
- **Composable** - Modules can be enabled/disabled independently
- **Maintainable** - Consistent structure across all modules

### Code Organization Best Practices

**Stan Guidelines:**
- Use built-in zero constructors: `zeros_vector()`, `zeros_int_array()`
- Don't pass array/vector sizes as arguments - use `size()` internally
- Ignore linter warnings about code section placement (modular `#include` architecture)

**R Guidelines:**
- Follow tidyverse style guide
- Use modern pipe operator `|>` (not `%>%`)
- Use `testthat` framework for unit tests

---

## References

### Key Documentation

- **Archived detailed design:** [`docs/multi_level_hierarchy_design.md`](multi_level_hierarchy_design.md)
- **Archived optimization:** [`docs/stan_state_space_optimization.md`](stan_state_space_optimization.md)
- **Other events model:** [`docs/OTHER_EVENTS_MODEL.md`](OTHER_EVENTS_MODEL.md)
- **Mechanistic covariates:** [`docs/other_events_mechanistic_rates.md`](other_events_mechanistic_rates.md)
- **State space model:** `sld_state_space_model.md`
- **Code ownership:** [`docs/CODEOWNERS`](CODEOWNERS)

### Configuration Files

- **Quarto config:** `_quarto.yml` for documentation generation
- **Targets config:** `_targets.yaml` for workflow configuration
- **R environment:** `renv.lock` for package versions

---

**Document Version:** 2.0 (Consolidated)  
**Status:** Active Reference  
**Maintainer:** Update when architectural patterns change
