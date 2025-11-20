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

## Time-Varying Process Noise (AR(1) Model)

**Added:** November 19, 2025  
**Status:** Implemented for tumor regression (`tr` module)

### Overview

The model includes **AR(1) time-varying process noise** to capture smooth, systematic deviations in tumor dynamics beyond the mechanistic mean trajectory. This separates:

- **Measurement error** (white noise): Independent observation-level noise
- **Process noise** (colored noise): Smoothly correlated temporal deviations in underlying rates

### Mathematical Formulation

**AR(1) Process:**
```
deviation[1] ~ normal(0, σ)
deviation[t] ~ normal(φ * deviation[t-1], σ * sqrt(1 - φ²))
```

**Applied to Rates (multiplicative in log-space):**
```
patient_rate[t] = exp(log(base_rate) + deviation[t])
```

**Parameters:**
- `φ` (phi): AR(1) coefficient (0 ≤ φ < 1)
  - φ = 0: Independent noise (white noise)
  - φ → 1: High temporal correlation (random walk)
  - Typical: φ ≈ 0.8-0.9 for smooth trends
- `σ` (sigma): Innovation standard deviation
  - Controls magnitude of deviations
  - Stationary variance: σ² / (1 - φ²)

### Hierarchical Structure

**Three-level hierarchy:**
```stan
// Population level
tr_log_sd_pop_process_noise      // log(σ) population mean
tr_phi_pop_process_noise         // φ population mean

// Patient-level variation
tr_sd_patient_log_sd_process_noise   // Between-patient SD in log(σ)
tr_sd_patient_phi_process_noise      // Between-patient SD in φ
tr_raw_patient_log_sd_process_noise  // Raw patient effects (NCP)
tr_raw_patient_phi_process_noise     // Raw patient effects (NCP)

// Time-varying deviations
tr_raw_patient_process_noise[i, t]   // Raw AR(1) innovations (NCP)
```

### Non-Centered Parameterization (NCP)

All process noise uses NCP for computational efficiency:

**Raw Parameters:**
```stan
matrix[n_patients, max_t_width] tr_raw_patient_process_noise;
tr_raw_patient_process_noise ~ std_normal();  // Prior
```

**Transformed to Actual Deviations:**
```stan
sigma = exp(patient_log_sd[i]);
phi = patient_phi[i];
sd_innovation = sigma * sqrt(1 - phi^2);

deviation[i, 1] = sigma * tr_raw_patient_process_noise[i, 1];
for (t in 2:max_t_width) {
  deviation[i, t] = phi * deviation[i, t-1] + 
                    sd_innovation * tr_raw_patient_process_noise[i, t];
}
```

**Benefits of NCP:**
- Reduces posterior correlation between parameters
- Improves HMC sampling efficiency
- Standard practice for weakly-identified hierarchical models

### Prior Design Philosophy

**Conservative priors** ensure process noise enhances rather than dominates the mechanistic model:

```r
# R prior hyperparameters (in priors.R)
tr_log_sd_pop_process_noise_sd     = 1.0   # σ range ~0.01-0.2
tr_phi_pop_process_noise_sd        = 0.2   # φ ≈ 0.8-0.9
tr_log_sd_patient_process_noise_sd = 0.5   # Moderate patient variation
tr_phi_patient_process_noise_sd    = 0.15  # Moderate phi variation
```

**Stan priors:**
```stan
tr_log_sd_pop_process_noise ~ normal(log(0.05), 1.0);  // ≈ 5% deviations
tr_phi_pop_process_noise ~ beta(20, 5);                // Centered at 0.8
```

### Computational Complexity

**Per Patient:**
- **Time:** O(T) where T = number of timepoints
- **Space:** O(T) for deviation vector
- **Efficient:** Sequential AR(1) computation is cache-friendly

**Total Model:**
- Adds 7 parameters per model (population + patient variation)
- Matrix `[n_patients × max_t_width]` for raw innovations
- Minimal overhead when disabled via feature flag

### Feature Flags

```r
# Enable/disable in stan_data
enable_patient_process_noise_tr = 1    # Implemented
enable_patient_process_noise_frac = 0  # Not yet implemented
```

### Implementation Status

**Currently Implemented:**
- ✅ `tr` module (tumor regression/decrease rates)
- ✅ Full hierarchical structure (population + patient variation)
- ✅ R priors and initialization
- ✅ Non-centered parameterization
- ✅ Measurement module separation

**Planned:**
- ⏳ `frac` module (growth fraction rates)
- ⏳ Empirical validation with trial data
- ⏳ Comparison of φ estimates across patients

### Module Organization

**Files Modified/Created:**
```
stan/ssls/modules/
├── tr/
│   ├── hyperparams.stan        # Added 4 process noise hyperparams
│   ├── parameters.stan         # Added 7 process noise params
│   ├── transformed_parameters.stan  # Added AR(1) computation
│   └── priors.stan            # Added process noise priors
└── measurement/               # New module
    ├── flags.stan
    ├── data.stan
    ├── hyperparams.stan       # measure_sd_sd
    ├── parameters.stan        # measure_sd
    ├── priors.stan
    └── transformed_data.stan  # lod, log_lod
```

### Design Rationale

**Why AR(1) instead of GP?**
- **Computational:** O(T) vs O(T³)
- **Stationarity:** Natural with |φ| < 1
- **Interpretability:** Two parameters vs complex kernel
- **Sufficient:** Captures smooth trends without overfitting
- **Priors:** Easier to elicit beliefs about temporal correlation

**Why multiplicative in log-space?**
- **Positivity:** Guarantees rates remain positive
- **Scale invariance:** Deviations proportional to rate magnitude  
- **Additivity:** Simpler mathematics in log-space
- **HMC geometry:** Better posterior geometry

**Why separate process and measurement noise?**
- **Identifiability:** Different roles (dynamics vs observation)
- **Interpretability:** Process = systematic trends, measurement = random error
- **Flexibility:** Can disable one without affecting the other

### Clinical Implications

**Process noise can capture:**
- Gradual treatment adaptation (acquired resistance)
- Tumor heterogeneity evolution over time
- Biological changes not in mechanistic model
- Patient-specific deviations from population trajectory

**Improves:**
- Uncertainty quantification for extrapolation
- Detection of unstable dynamics
- Prediction of future measurements

### References

- **Implementation details:** [`docs/CHANGELOG.md`](CHANGELOG.md) (November 19, 2025 entry)
- **Module structure:** This document, Module Structure section
- **Hierarchical design:** This document, Multi-Level Hierarchical Parameter Design

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
