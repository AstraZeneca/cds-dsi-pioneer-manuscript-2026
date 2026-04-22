# Propensity-Weighted Bayesian Borrowing — Implementation Spec

## Summary

Replace the trial-level random effect in the pioneer combined model with a
joint Bayesian propensity submodel. The propensity model learns
`P(target | covariates)` via logistic regression and uses it to weight each
non-target patient's likelihood contribution. Target patients always contribute
with weight 1.0.

### Design Decisions (from brainstorming)

| Decision | Choice | Rationale |
|----------|--------|-----------|
| Weight function | Raw propensity `w = P(target \| x)` | Simplest; generalizable later |
| Covariate space | Original `covar_design_matrix` (not QR) | Logistic regression doesn't need QR; keeps coefficients interpretable |
| Weighted likelihoods | Both PSA and multistate | Consistent influence control |
| Multistate weighting | Add weight vector to `multistate_lpmf` | Function already loops per-patient internally |
| Target identification | Level-agnostic: `propensity_split_level` + `propensity_target_group` | Mirrors Laplace pattern; works for any hierarchy level |

## New Module: `stan/modules/propensity/`

### `flags.stan`

```stan
int<lower=0, upper=1> enable_propensity_weighting;
int<lower=0> propensity_split_level;   // 0 = disabled, >0 = hierarchy level to split on
int<lower=0> propensity_target_group;  // group ID at that level (weight = 1.0)
```

### `hyperparams.stan`

Prior SDs for the propensity logistic regression coefficients. Set from R to
control the prior predictive distribution of propensity scores.

```stan
real<lower=0> propensity_intercept_sd;  // SD for intercept prior
real<lower=0> propensity_coef_sd;       // SD for coefficient priors
```

With `p` standardized covariates, the prior predictive linear predictor has
`SD ≈ sqrt(intercept_sd² + p * coef_sd²)`. To keep prior probabilities
away from 0/1 (avoiding a U-shape), target a logit-scale SD of ~2.
Default: `intercept_sd = 1.0`, `coef_sd = 2 / sqrt(n_covar)`.

For `n_covar = 16`: `coef_sd ≈ 0.5`, giving `logit SD ≈ sqrt(1 + 16 * 0.25) ≈ 2.2`.

### `parameters.stan`

```stan
// Propensity logistic regression coefficients
// Guarded by enable_propensity_weighting to avoid declaring zero-sized parameters
array[enable_propensity_weighting ? 1 : 0] real beta_propensity_intercept;
vector[enable_propensity_weighting ? n_covar : 0] beta_propensity;
```

Note: `beta_propensity_intercept` is needed because we use the original
`covar_design_matrix` (centered but not QR-transformed), which does not absorb
the intercept. Uses modern Stan array syntax (`array[N] real` not `real[N]`).

### `transformed_data.stan`

Compute target patient range and binary indicator from the hierarchy structure.
Level-agnostic — no dependency on `trial_patient_pos`.

```stan
// Target patient range (contiguous by construction of combined data)
int propensity_target_start = 1;
int propensity_target_end = n_patients;
int propensity_n_target = n_patients;

// Binary indicator for bernoulli likelihood (computed once at init)
array[n_patients] int<lower=0, upper=1> propensity_is_target = rep_array(1, n_patients);

if (enable_propensity_weighting && propensity_split_level > 0) {
  // Find contiguous range of target patients
  propensity_target_start = 0;
  propensity_target_end = 0;
  int n_found = 0;
  for (i in 1:n_patients) {
    if (patient_level_groups[i, propensity_split_level] == propensity_target_group) {
      if (propensity_target_start == 0) propensity_target_start = i;
      propensity_target_end = i;
      n_found += 1;
    }
  }

  propensity_n_target = n_found;

  // Validate contiguity
  if (propensity_target_end - propensity_target_start + 1 != n_found)
    fatal_error("Propensity target group patients are not contiguous: ",
                "found ", n_found, " patients in range [",
                propensity_target_start, ", ", propensity_target_end, "]");

  // Build binary indicator
  propensity_is_target = zeros_int_array(n_patients);
  for (i in propensity_target_start:propensity_target_end) {
    propensity_is_target[i] = 1;
  }

  print("Propensity weighting enabled:");
  print("  split_level = ", propensity_split_level);
  print("  target_group = ", propensity_target_group);
  print("  target patients: [", propensity_target_start, ", ", propensity_target_end,
        "] (n=", propensity_n_target, ")");
  print("  non-target patients: ", n_patients - propensity_n_target);
}
```

### `transformed_parameters.stan`

Vectorized computation of likelihood weights. No per-iteration loops.

```stan
vector<lower=0, upper=1>[n_patients] likelihood_weight = ones_vector(n_patients);

if (enable_propensity_weighting && propensity_split_level > 0) {
  // Compute propensity scores for all patients (vectorized)
  likelihood_weight = inv_logit(
    beta_propensity_intercept[1] + covar_design_matrix * beta_propensity
  );
  // Target patients always contribute fully (range assignment, no loop)
  likelihood_weight[propensity_target_start:propensity_target_end] =
    ones_vector(propensity_n_target);
}
```

### `priors.stan`

Joint propensity submodel: the bernoulli likelihood informs the propensity
coefficients, while outcome likelihoods provide feedback through the shared
`likelihood_weight`.

```stan
if (enable_propensity_weighting && propensity_split_level > 0) {
  // Priors scaled to avoid U-shaped prior predictive probabilities
  beta_propensity_intercept ~ normal(0, propensity_intercept_sd);
  beta_propensity ~ normal(0, propensity_coef_sd);

  // Propensity likelihood: P(target | covariates)
  propensity_is_target ~ bernoulli_logit(
    beta_propensity_intercept[1] + covar_design_matrix * beta_propensity
  );
}
```

## Changes to `pioneer.stan`

### Includes (data, transformed data, parameters, transformed parameters, priors)

Add propensity module includes to each block:

```stan
data {
  // ... existing includes ...
  #include "modules/propensity/flags.stan"
  #include "modules/propensity/hyperparams.stan"
}

transformed data {
  // ... existing includes ...
  #include "modules/propensity/transformed_data.stan"
}

parameters {
  // ... existing includes ...
  #include "modules/propensity/parameters.stan"
}

transformed parameters {
  // ... existing includes ...
  #include "modules/propensity/transformed_parameters.stan"
}

model {
  // ... existing priors ...
  #include "modules/propensity/priors.stan"
  // ... likelihoods (modified below) ...
}
```

### Model Block — PSA Likelihood (lines 91–110)

Convert from `~` sampling statement to `target +=` with weight:

```stan
if (fit_psa_data) {
  profile("psa loglik") {
    for (j in 1:n_hmc_patients) {
      int p = hmc_patient_idx[j];
      int n_screen = n_patient_screening_visits[p];
      int data_start = patient_visit_pos[p];
      int data_end   = patient_visit_pos[p + 1] - 1;
      int state_start = hmc_visit_pos[j];

      target += likelihood_weight[p] * sf_log_space_obs_lpdf(
        normalized_psa[data_start + n_screen : data_end] |
        states[state_start + n_screen : state_start + (data_end - data_start)],
        measure_sd_psa,
        negative_infinity(),
        measure_nu
      );
    }
  }
}
```

### Model Block — Multistate Likelihood (lines 113–142)

Simplify: remove the redundant single-transition branch. `multistate_lpmf`
already has an internal fast path for 0→1-only configurations (line 158 of
`multistate.stanfunctions`), so the outer `if/else` dispatch in the model block
is unnecessary. This simplification also avoids having to maintain two weighted
code paths.

```stan
if (fit_multistate_data) {
  profile("multistate loglik") {
    ms_final_state[hmc_patient_idx] ~ multistate(
      likelihood_weight[hmc_patient_idx],   // NEW: propensity weight
      enable_ms_01, enable_ms_02, enable_ms_12, ms_time_scale_12,
      enable_ms_03, enable_ms_32,
      ms_time_01[hmc_patient_idx], ms_time_02[hmc_patient_idx], ms_time_12[hmc_patient_idx],
      ms_time_03[hmc_patient_idx], ms_time_32[hmc_patient_idx],
      ms_censored_01[hmc_patient_idx], ms_censored_02[hmc_patient_idx], ms_censored_12[hmc_patient_idx],
      ms_censored_32[hmc_patient_idx],
      ms_prog_deterministic[hmc_patient_idx],
      ms_ic_gap_01[hmc_patient_idx],
      t_patient_visits,
      patient_visit_pos,
      log_cond_surv_01,
      log_cond_surv_02,
      log_cond_surv_12_s,
      log_cond_surv_12_t,
      log_cond_surv_03,
      log_cond_surv_32
    );
  }
}
```

## Changes to `multistate.stanfunctions`

### `multistate_lpmf` — add weight parameter

Add `vector weight` as the second argument (first after the required `int[]`
for the `~` syntax). The internal per-patient loop multiplies each patient's
log-likelihood by their weight.

```stan
real multistate_lpmf(
  array[] int final_state,
  vector weight,              // NEW: per-patient likelihood weight
  int enable_01, int enable_02, int enable_12, int ms_time_scale_12,
  int enable_03, int enable_32,
  // ... rest of existing signature unchanged ...
) {
  // Single-transition fast path (also weighted)
  if (enable_01 && !enable_02 && !enable_12 && !enable_03) {
    return dot_product(weight,
      calc_ms_single_transition_loglik(time_01, censored_01, log_cond_surv_01));
  }

  int n_patients = size(final_state);
  real total_ll = 0;

  for (i in 1:n_patients) {
    real patient_ll = 0;

    // ... existing per-patient logic completely unchanged ...

    total_ll += weight[i] * patient_ll;  // CHANGED: was total_ll += patient_ll
  }
  return total_ll;
}
```

Note: the internal single-transition fast path now uses `dot_product(weight, ...)`
instead of `sum(...)`, so the outer model block no longer needs to duplicate this
dispatch — it always calls `multistate_lpmf` and the function handles both cases.

### All other `multistate_lpmf` call sites

Every call site must pass a weight vector to match the updated signature.
Models without propensity weighting pass `ones_vector(n)`.

| File | Call pattern | Weight argument |
|------|-------------|-----------------|
| `stan/psa/ms-standalone.stan` | `ms_final_state[hmc_patient_idx] ~ multistate(...)` | `ones_vector(n_hmc_patients)` |
| `stan/tumor/ms-standalone.stan` | `ms_final_state ~ multistate(...)` | `ones_vector(n_patients)` |
| `stan/tumor/sf-ssm-log-space.stan` | `ms_final_state ~ multistate(...)` | `ones_vector(n_patients)` |
| `stan/tumor/sf-ssls-lfo.stan` | No `multistate_lpmf` call (includes functions only) | N/A |

Also update `calc_ms_single_transition_loglik` call sites in non-pioneer
models to use `sum(...)` (unchanged — they don't need weighting):

| File | Pattern | Change needed |
|------|---------|---------------|
| `stan/psa/ms-standalone.stan` | `target += sum(calc_ms_single_transition_loglik(...))` | None |
| `stan/tumor/sf-ssm-log-space.stan` | Handled inside `multistate_lpmf` | None |
| `stan/tumor/sf-ssls-lfo.stan` | `target += sum(calc_ms_single_transition_loglik(...))` | None |

## Changes to R Code

### `r/pioneer/prepare_analysis_data.R`

In `prepare_pioneer_stan_data()`, no changes needed — `is_trial` is derived
in Stan from `patient_level_groups` and `propensity_split_level/target_group`.

### `r/pioneer/initializers.R`

In `create_pioneer_initializer()`, add propensity parameter initialization.
Guard must match Stan's conditional sizing: intercept exists when
`enable_propensity_weighting == 1`, coefficients exist when additionally
`n_covar > 0`.

```r
# Inside the biomarker_init list:
beta_propensity_intercept = if (enable_propensity_weighting) array(0, dim = 1),
beta_propensity = if (enable_propensity_weighting) {
  if (n_covar > 0) rep(0, n_covar)
},
```

Zero initialization is appropriate: it starts all weights at 0.5, which is a
neutral starting point. The propensity likelihood will quickly move coefficients
toward the MLE during warmup.

### `r/pioneer/priors.R`

Add propensity hyperparameter defaults, scaled by `n_covar`:

```r
# In get_pioneer_priors() or equivalent:
propensity_intercept_sd = 1.0,
propensity_coef_sd = if (stan_data$n_covar > 0) 2.0 / sqrt(stan_data$n_covar) else 1.0,
```

For `n_covar = 16`: `coef_sd = 0.5`. This gives a prior predictive logit SD
of ~2.2, keeping most prior propensity scores between ~0.1 and ~0.9.

### `targets/pioneer_targets.R`

Add a new model variant row in the `tar_map` tribble:

```r
tribble(
  ~model,              ~hist,  ~model_formula,        ~trial_level, ~arm_level, ~propensity,
  "combined",          TRUE,   covar_formula,         TRUE,         TRUE,       FALSE,
  "combined_no_trial", TRUE,   covar_formula,         FALSE,        TRUE,       FALSE,
  "no_hist",           FALSE,  covar_formula,         FALSE,        TRUE,       FALSE,
  "combined_safe_covar", TRUE, safe_covar_formula,    TRUE,         TRUE,       FALSE,
  "propensity",        TRUE,   covar_formula,         FALSE,        TRUE,       TRUE,
),
```

And in `default_settings`, add the propensity flags:

```r
enable_propensity_weighting = as.integer(propensity),
propensity_split_level = if (propensity) 1L else 0L,
propensity_target_group = if (propensity) 1L else 0L,
```

The `propensity` model uses `trial_level = FALSE` (no trial-level RE) because
propensity weighting replaces the trial-level hierarchy. It keeps
`arm_level = TRUE` for within-trial arm structure.

## What Is NOT Affected

- **Priors on PSA/multistate/tr/frac/init parameters** — unweighted. Propensity
  controls how much each patient's DATA influences the posterior, not the
  hierarchical shrinkage structure.
- **Generated quantities** — all endpoint computations (KM, PFS, OS, ORR) use
  HMC patients unweighted. The weights only affect parameter estimation.
  `likelihood_weight` is declared in `transformed parameters` so it IS saved to
  CSV output — this is intentional for posterior diagnostics (inspecting which
  RWD patients get high/low weight).
- **Standalone multistate models** (`psa/ms-standalone.stan`,
  `tumor/ms-standalone.stan`) — pass `ones_vector` (no weighting). Need
  signature update only.
- **SCLC tumor models** (`tumor/sf-ssm-log-space.stan`,
  `tumor/sf-ssls-lfo.stan`) — pass `ones_vector` where they call
  `multistate_lpmf`. No propensity module or weighting logic added.
- **Laplace marginalization** — orthogonal. Propensity replaces the trial-level
  RE; Laplace marginalizes patient-level effects. They can be combined in the
  future.
- **QR decomposition** — unchanged. Propensity uses the original covariate
  matrix; QR is used only by outcome models.

### Edge case: `n_covar = 0`

When no covariates are available, `beta_propensity` is a zero-length vector and
`covar_design_matrix * beta_propensity` produces a zero vector. The propensity
score for all patients becomes `inv_logit(beta_propensity_intercept[1])` — a
single shared weight. This is valid but uninformative; the model degenerates to
uniform downweighting of RWD. The `bernoulli_logit` likelihood still identifies
the intercept (base rate of trial membership).

## New Parameter Count

| Parameter | Size | Notes |
|-----------|------|-------|
| `beta_propensity_intercept` | 1 | Propensity logit intercept |
| `beta_propensity` | `n_covar` (~16) | Propensity logit coefficients |
| **Total** | **~17** | Logistic regression is convex; trivial for HMC |

## Testing Plan

1. **Syntax check**: `stanc` compilation with propensity module included
2. **Disabled mode**: Verify `enable_propensity_weighting = 0` produces identical
   results to `combined_no_trial`
3. **Prior predictive**: Run with `fit_psa_data = 0, fit_multistate_data = 0` to
   verify propensity coefficients recover the data-generating process
4. **Posterior comparison**: Compare `propensity` vs `combined_no_trial` vs
   `combined` posterior predictive KM curves
5. **Weight inspection**: Extract `likelihood_weight` from posterior to verify
   RWD patients near the trial covariate space get higher weights

## File Change Summary

| File | Action | Complexity |
|------|--------|------------|
| `stan/modules/propensity/flags.stan` | Create | Small |
| `stan/modules/propensity/hyperparams.stan` | Create | Small |
| `stan/modules/propensity/parameters.stan` | Create | Small |
| `stan/modules/propensity/transformed_data.stan` | Create | Small |
| `stan/modules/propensity/transformed_parameters.stan` | Create | Small |
| `stan/modules/propensity/priors.stan` | Create | Small |
| `stan/psa/pioneer.stan` | Edit: add includes + weight PSA/MS likelihoods | Medium |
| `stan/multistate.stanfunctions` | Edit: add weight arg to `multistate_lpmf` | Medium |
| `stan/psa/ms-standalone.stan` | Edit: pass `ones_vector` to `multistate_lpmf` | Trivial |
| `stan/tumor/ms-standalone.stan` | Edit: pass `ones_vector` to `multistate_lpmf` | Trivial |
| `stan/tumor/sf-ssm-log-space.stan` | Edit: pass `ones_vector` to `multistate_lpmf` | Trivial |
| `r/pioneer/initializers.R` | Edit: add propensity inits | Trivial |
| `r/pioneer/priors.R` | Edit: add propensity hyperparameter defaults | Trivial |
| `r/pioneer/prepare_analysis_data.R` | Edit: add propensity flags to `prepare_ms_standalone_stan_data()` | Trivial |
| `targets/pioneer_targets.R` | Edit: add propensity model variant + flags | Small |
