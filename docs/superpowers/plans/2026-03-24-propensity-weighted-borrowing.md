# Propensity-Weighted Bayesian Borrowing Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a joint propensity submodel to `pioneer.stan` that weights each RWD patient's PSA and multistate likelihood contribution by their estimated probability of being trial-eligible, replacing the trial-level random effect.

**Architecture:** A new `stan/modules/propensity/` module (5 files) learns `P(target | covariates)` via logistic regression jointly with the outcome model. Per-patient weights flow into both the PSA loop (`target += w * lpdf(...)`) and `multistate_lpmf` (new `weight` argument). The `propensity_split_level` / `propensity_target_group` pattern mirrors Laplace routing — no hardcoded trial assumptions.

**Tech Stack:** Stan 2.38, R tidyverse, `targets`, CmdStanR.

---

## File Map

| File | Action |
|------|--------|
| `stan/multistate.stanfunctions` | Edit: add `vector weight` param to `multistate_lpmf`; use `dot_product` in fast path |
| `tests/testthat/stan/test_multistate_loglik_all.stan` | Edit: add `weight` to wrapper + data block |
| `tests/testthat/test-stan-multistate-loglik.R` | Edit: pass `weight = rep(1, n)` in `make_ms_data`, add weighted test |
| `stan/tumor/sf-ssm-log-space.stan` | Edit: pass `ones_vector(n_patients)` to `multistate_lpmf` |
| `stan/tumor/ms-standalone.stan` | Edit: pass `ones_vector(n_patients)` to `multistate_lpmf` |
| `stan/psa/ms-standalone.stan` | Edit: simplify if/else to single call + `ones_vector(n_hmc_patients)` |
| `tests/testthat/stan/test_propensity_weights.stan` | Create: test harness for propensity module |
| `tests/testthat/test-stan-propensity.R` | Create: range detection + weight value tests |
| `stan/modules/propensity/flags.stan` | Create |
| `stan/modules/propensity/hyperparams.stan` | Create |
| `stan/modules/propensity/parameters.stan` | Create |
| `stan/modules/propensity/transformed_data.stan` | Create |
| `stan/modules/propensity/transformed_parameters.stan` | Create |
| `stan/modules/propensity/priors.stan` | Create |
| `stan/psa/pioneer.stan` | Edit: add includes; convert PSA `~` to `target +=`; simplify multistate block |
| `r/pioneer/priors.R` | Edit: add propensity hyperparameter defaults |
| `r/pioneer/initializers.R` | Edit: add `beta_propensity_intercept` / `beta_propensity` inits |
| `r/pioneer/prepare_analysis_data.R` | Edit: add propensity flags to `prepare_ms_standalone_stan_data()` |
| `targets/pioneer_targets.R` | Edit: add `propensity` model variant + flags |

---

## Task 1: Add `weight` to `multistate_lpmf` and update tests

This is the highest-risk change — it touches shared infrastructure used by all Stan models. Do it first, verify tests pass, then fix the callers.

**Files:**
- Modify: `stan/multistate.stanfunctions:138-160`
- Modify: `tests/testthat/stan/test_multistate_loglik_all.stan`
- Modify: `tests/testthat/test-stan-multistate-loglik.R`

- [ ] **Step 1: Update `multistate_lpmf` signature in `stan/multistate.stanfunctions`**

Change lines 138–160. Add `vector weight` as the second argument (first after `final_state`), update the fast path to use `dot_product`, and update the accumulation line:

```stan
real multistate_lpmf(
  array[] int final_state,
  vector weight,              // NEW: per-patient likelihood weight [0, 1]
  int enable_01, int enable_02, int enable_12, int ms_time_scale_12,
  int enable_03, int enable_32,
  array[] int time_01, array[] int time_02, array[] int time_12,
  array[] int time_03, array[] int time_32,
  array[] int censored_01, array[] int censored_02, array[] int censored_12,
  array[] int censored_32,
  array[] int prog_deterministic,
  array[] int ms_ic_gap_01,
  array[] int t_patient_visits,
  array[] int patient_visit_pos,
  matrix log_cond_surv_01,
  matrix log_cond_surv_02,
  matrix log_cond_surv_12_s,
  matrix log_cond_surv_12_t,
  matrix log_cond_surv_03,
  matrix log_cond_surv_32
) {
  // Single-transition fast path (also weighted)
  if (enable_01 && !enable_02 && !enable_12 && !enable_03) {
    return dot_product(weight,
      calc_ms_single_transition_loglik(time_01, censored_01, log_cond_surv_01));
  }
  // ...
  total_ll += weight[i] * patient_ll;  // only this accumulation line changes
```

The only changes in the body are: (1) fast path `sum(...)` → `dot_product(weight, ...)`, (2) `total_ll += patient_ll` → `total_ll += weight[i] * patient_ll`.

- [ ] **Step 2: Update the test Stan wrapper `tests/testthat/stan/test_multistate_loglik_all.stan`**

Add `weight` to the data block, the wrapper function signature, and every `multistate_loglik(...)` call in generated quantities:

```stan
// In data block, add:
vector[N] weight;

// Wrapper function signature: add vector weight as second arg
real multistate_loglik(
  array[] int final_state,
  vector weight,           // NEW
  int enable_01, ...
) {
  return multistate_lpmf(
    final_state |
    weight,                // NEW
    enable_01, ...
  );
}

// In generated quantities, update both calls:
real ms_lpmf_01only = multistate_loglik(
  final_state,
  weight,                  // NEW
  1, 0, 0, 0, ...
);
real ms_lpmf_01_02 = multistate_loglik(
  final_state,
  weight,                  // NEW
  1, 1, 0, 0, ...
);
```

- [ ] **Step 3: Update `make_ms_data` in `tests/testthat/test-stan-multistate-loglik.R`**

Add `weight = rep(1.0, n)` to the returned list in `make_ms_data()`:

```r
make_ms_data <- function(n, max_t, event_time, censored, final_state,
                         lcs_01, lcs_02 = NULL, time_02 = NULL, censored_02 = NULL) {
  # ... existing code ...
  list(
    # ... existing fields ...
    weight = rep(1.0, n)   # NEW: all-ones for backward-compatible tests
  )
}
```

- [ ] **Step 4: Add a weighted test to `test-stan-multistate-loglik.R`**

Verify that halving one patient's weight halves their contribution:

```r
test_that("multistate_lpmf: weight=0.5 halves patient contribution", {
  n <- 2L; max_t <- 8L
  lcs <- matrix(-0.08, nrow = n, ncol = max_t)
  data_full <- make_ms_data(n, max_t, c(4L, 6L), c(0L, 1L),
                            final_state = c(1L, 0L), lcs_01 = lcs)
  data_half <- modifyList(data_full, list(weight = c(0.5, 1.0)))

  fit_full <- test_stan_function(
    here("tests", "testthat", "stan", "test_multistate_loglik_all.stan"),
    data_full
  )
  fit_half <- test_stan_function(
    here("tests", "testthat", "stan", "test_multistate_loglik_all.stan"),
    data_half
  )

  d_full <- posterior::as_draws_df(fit_full$draws())
  d_half <- posterior::as_draws_df(fit_half$draws())

  # Patient 1 contributes ms_lpmf_01only_full − ms_lpmf_01only_half
  # That difference should equal 0.5 * patient_1_ll
  full_ll  <- get_stan_val(d_full, "ms_lpmf_01only")
  half_ll  <- get_stan_val(d_half, "ms_lpmf_01only")
  p1_ll    <- get_stan_val(d_full, "st_llik", 1)
  expect_equal(full_ll - half_ll, 0.5 * p1_ll, tolerance = 1e-6)
})
```

- [ ] **Step 5: Run tests**

```bash
cd /mnt/code/.worktrees/karim/propensity-borrowing
Rscript -e 'testthat::test_file("tests/testthat/test-stan-multistate-loglik.R")'
```

Expected: all tests pass including the new weighted test.

- [ ] **Step 6: Commit**

```bash
git add stan/multistate.stanfunctions \
        tests/testthat/stan/test_multistate_loglik_all.stan \
        tests/testthat/test-stan-multistate-loglik.R
git commit -m "feat: add per-patient weight to multistate_lpmf"
```

---

## Task 2: Update non-pioneer `multistate_lpmf` callers

Three Stan models call `multistate_lpmf` but don't use propensity weighting. Pass `ones_vector(n)` to match the new signature.

**Files:**
- Modify: `stan/tumor/sf-ssm-log-space.stan:81`
- Modify: `stan/tumor/ms-standalone.stan:114`
- Modify: `stan/psa/ms-standalone.stan:81-106`

- [ ] **Step 1: Update `stan/tumor/sf-ssm-log-space.stan`**

At line 81, insert `ones_vector(n_patients)` as the second argument:

```stan
ms_final_state ~ multistate(
  ones_vector(n_patients),   // no propensity weighting in tumor model
  enable_ms_01, enable_ms_02, enable_ms_12, ms_time_scale_12,
  ...
```

- [ ] **Step 2: Update `stan/tumor/ms-standalone.stan`**

At line 114, insert `ones_vector(n_patients)`:

```stan
ms_final_state ~ multistate(
  ones_vector(n_patients),   // no propensity weighting in standalone
  enable_ms_01, enable_ms_02, enable_ms_12, ms_time_scale_12,
  ...
```

- [ ] **Step 3: Update `stan/psa/ms-standalone.stan`**

At lines 81–106, simplify the existing if/else (same redundant dispatch as in `pioneer.stan`) to a single call, passing `ones_vector(n_hmc_patients)`:

```stan
if (fit_multistate_data) {
  profile("multistate loglik") {
    ms_final_state[hmc_patient_idx] ~ multistate(
      ones_vector(n_hmc_patients),   // no propensity weighting in standalone
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

- [ ] **Step 4: Syntax-check all three models**

```bash
cd /mnt/code/.worktrees/karim/propensity-borrowing
~/.cmdstan/cmdstan-2.38.0/bin/stanc --include-paths=stan --include-paths=stan/tumor \
  stan/tumor/sf-ssm-log-space.stan

~/.cmdstan/cmdstan-2.38.0/bin/stanc --include-paths=stan --include-paths=stan/tumor \
  stan/tumor/ms-standalone.stan

~/.cmdstan/cmdstan-2.38.0/bin/stanc --include-paths=stan --include-paths=stan/psa \
  stan/psa/ms-standalone.stan
```

Expected: `Model name: ...` with no errors for all three.

- [ ] **Step 5: Commit**

```bash
git add stan/tumor/sf-ssm-log-space.stan \
        stan/tumor/ms-standalone.stan \
        stan/psa/ms-standalone.stan
git commit -m "feat: pass ones_vector weight to multistate_lpmf in non-pioneer models"
```

---

## Task 3: Create `stan/modules/propensity/` module

Five small files. Create them all, then do a single syntax check via `pioneer.stan` in Task 4.

**Files:**
- Create: `stan/modules/propensity/flags.stan`
- Create: `stan/modules/propensity/hyperparams.stan`
- Create: `stan/modules/propensity/parameters.stan`
- Create: `stan/modules/propensity/transformed_data.stan`
- Create: `stan/modules/propensity/transformed_parameters.stan`
- Create: `stan/modules/propensity/priors.stan`

- [ ] **Step 1: Create `stan/modules/propensity/flags.stan`**

```stan
// modules/propensity/flags.stan
// Feature flags for propensity-weighted borrowing.
// propensity_split_level = 0 disables the module entirely.
int<lower=0, upper=1> enable_propensity_weighting;
int<lower=0> propensity_split_level;   // >0: which hierarchy level to split on
int<lower=0> propensity_target_group;  // group ID at that level (weight = 1.0)
```

- [ ] **Step 2: Create `stan/modules/propensity/hyperparams.stan`**

```stan
// modules/propensity/hyperparams.stan
// Prior SDs for propensity logistic regression.
// With p standardized covariates, target logit-SD ≈ 2 (avoids U-shaped scores).
// Recommended defaults (set from R): intercept_sd=1.0, coef_sd=2/sqrt(n_covar).
real<lower=0> propensity_intercept_sd;
real<lower=0> propensity_coef_sd;
```

- [ ] **Step 3: Create `stan/modules/propensity/parameters.stan`**

```stan
// modules/propensity/parameters.stan
// Propensity logistic regression coefficients.
// Sized 0 when weighting is disabled to avoid sampling unused parameters.
// Explicit intercept needed because covar_design_matrix is centered (not QR).
array[enable_propensity_weighting ? 1 : 0] real beta_propensity_intercept;
vector[enable_propensity_weighting ? n_covar : 0] beta_propensity;
```

- [ ] **Step 4: Create `stan/modules/propensity/transformed_data.stan`**

```stan
// modules/propensity/transformed_data.stan
// Locate the contiguous range of target patients using the level-agnostic
// patient_level_groups matrix. Runs once at initialization.
int propensity_target_start = 1;
int propensity_target_end = n_patients;
int propensity_n_target = n_patients;

// Binary indicator: 1 = target patient (weight=1), 0 = non-target (weight=propensity score)
array[n_patients] int<lower=0, upper=1> propensity_is_target = rep_array(1, n_patients);

if (enable_propensity_weighting && propensity_split_level > 0) {
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

  // Validate that target patients are contiguous (required for range assignment)
  if (propensity_target_end - propensity_target_start + 1 != n_found)
    fatal_error("Propensity target group patients are not contiguous: ",
                "found ", n_found, " patients but range [",
                propensity_target_start, ", ", propensity_target_end, "] has ",
                propensity_target_end - propensity_target_start + 1, " slots");

  propensity_is_target = zeros_int_array(n_patients);
  for (i in propensity_target_start:propensity_target_end) {
    propensity_is_target[i] = 1;
  }

  print("Propensity weighting enabled:");
  print("  split_level=", propensity_split_level,
        " target_group=", propensity_target_group);
  print("  target: [", propensity_target_start, ", ", propensity_target_end,
        "] n=", propensity_n_target,
        " non-target: ", n_patients - propensity_n_target);
}
```

- [ ] **Step 5: Create `stan/modules/propensity/transformed_parameters.stan`**

```stan
// modules/propensity/transformed_parameters.stan
// Compute per-patient likelihood weights.
// Vectorized: matrix-vector multiply for all patients; range assignment for targets.
//
// NOTE: uses covar_design_matrix (original centered/scaled matrix) NOT Q_covar_design_matrix.
// All other modules use the QR-decomposed Q matrix for numerical stability, but propensity
// logistic regression is convex and doesn't need QR. Using the original matrix keeps
// beta_propensity coefficients interpretable (log-odds per SD change per covariate).
vector<lower=0, upper=1>[n_patients] likelihood_weight = ones_vector(n_patients);

if (enable_propensity_weighting && propensity_split_level > 0) {
  likelihood_weight = inv_logit(
    beta_propensity_intercept[1] + covar_design_matrix * beta_propensity
  );
  // Override target patients to weight 1.0 (range assignment, no loop)
  likelihood_weight[propensity_target_start:propensity_target_end] =
    ones_vector(propensity_n_target);
}
```

- [ ] **Step 6: Create `stan/modules/propensity/priors.stan`**

```stan
// modules/propensity/priors.stan
// Joint propensity submodel: P(target | covariates).
// The bernoulli likelihood identifies beta_propensity from the data source split.
// The outcome likelihoods provide feedback via likelihood_weight.
if (enable_propensity_weighting && propensity_split_level > 0) {
  // Priors scaled to avoid U-shaped prior predictive: coef_sd = 2/sqrt(n_covar)
  beta_propensity_intercept ~ normal(0, propensity_intercept_sd);
  beta_propensity ~ normal(0, propensity_coef_sd);

  // Propensity likelihood
  propensity_is_target ~ bernoulli_logit(
    beta_propensity_intercept[1] + covar_design_matrix * beta_propensity
  );
}
```

- [ ] **Step 7: Commit module files**

```bash
git add stan/modules/propensity/
git commit -m "feat: add propensity module (flags, hyperparams, parameters, transformed_data, transformed_parameters, priors)"
```

---

## Task 4: Test the propensity module

The range-finding logic in `transformed_data.stan` (hierarchy indexing, contiguity
validation) and weight computation in `transformed_parameters.stan` need explicit
tests before wiring into `pioneer.stan`.

**Files:**
- Create: `tests/testthat/stan/test_propensity_weights.stan`
- Create: `tests/testthat/test-stan-propensity.R`

- [ ] **Step 1: Create `tests/testthat/stan/test_propensity_weights.stan`**

The test program declares hierarchy data and propensity betas as data (to fix
values without sampling), includes `transformed_data.stan` in the `transformed
data` block, and inlines the `transformed_parameters.stan` logic in GQ (block
includes cannot go in GQ, so we copy the computation there).

```stan
// Test harness for stan/modules/propensity/ transformed_data and
// transformed_parameters logic.
// betas are passed as data to fix their values without sampling.
functions {
  #include "util.stanfunctions"
  #include "pos.stanfunctions"
}
data {
  // Hierarchy (minimal: n_levels levels, identity mapping at last level)
  int<lower=1> n_patients;
  int<lower=1> n_levels;
  array[n_patients, n_levels] int patient_level_groups;

  int<lower=0> n_covar;
  matrix[n_patients, n_covar] covar_design_matrix;

  // Propensity module flags / hyperparams
  int<lower=0, upper=1> enable_propensity_weighting;
  int<lower=0> propensity_split_level;
  int<lower=0> propensity_target_group;
  real<lower=0> propensity_intercept_sd;
  real<lower=0> propensity_coef_sd;

  // Propensity betas — passed as data so values are fixed for testing
  real beta_propensity_intercept_1;  // will become beta_propensity_intercept[1]
  vector[n_covar] beta_propensity;
}
transformed data {
  // Replicate parameters.stan conditional sizing using data fields
  array[enable_propensity_weighting ? 1 : 0] real beta_propensity_intercept;
  if (enable_propensity_weighting) {
    beta_propensity_intercept[1] = beta_propensity_intercept_1;
  }

  // Module under test: sets propensity_target_start/end/n_target, propensity_is_target
  #include "modules/propensity/transformed_data.stan"
}
generated quantities {
  // Outputs from transformed_data.stan
  int out_target_start           = propensity_target_start;
  int out_target_end             = propensity_target_end;
  int out_n_target               = propensity_n_target;
  array[n_patients] int out_is_target = propensity_is_target;

  // transformed_parameters.stan logic inlined (block includes cannot go in GQ)
  vector[n_patients] out_likelihood_weight = ones_vector(n_patients);
  if (enable_propensity_weighting && propensity_split_level > 0) {
    out_likelihood_weight = inv_logit(
      beta_propensity_intercept[1] + covar_design_matrix * beta_propensity
    );
    out_likelihood_weight[propensity_target_start:propensity_target_end] =
      ones_vector(propensity_n_target);
  }
}
```

- [ ] **Step 2: Create `tests/testthat/test-stan-propensity.R`**

```r
library(testthat)
library(here)

# 5 patients: 3 trial (group 1 at level 1), 2 RWD (group 2).
# Two hierarchy levels: (trial, patient). Patient level is identity.
make_propensity_data <- function(
    n_patients     = 5L,
    n_trial        = 3L,
    n_covar        = 2L,
    enable         = 1L,
    split_level    = 1L,
    target_group   = 1L,
    beta_intercept = 0.0,
    beta           = NULL) {
  n_rwd <- n_patients - n_trial
  # Level 1: trial group (1 = trial, 2 = RWD)
  # Level 2: patient identity (1..n_patients)
  trial_group  <- c(rep(1L, n_trial), rep(2L, n_rwd))
  patient_id   <- seq_len(n_patients)
  level_groups <- cbind(trial_group, patient_id)

  # Simple covariate matrix: each patient gets a unique row
  X <- matrix(seq_len(n_patients * n_covar) / 10, nrow = n_patients, ncol = n_covar)

  list(
    n_patients                  = n_patients,
    n_levels                    = 2L,
    patient_level_groups        = level_groups,
    n_covar                     = n_covar,
    covar_design_matrix         = X,
    enable_propensity_weighting = enable,
    propensity_split_level      = split_level,
    propensity_target_group     = target_group,
    propensity_intercept_sd     = 1.0,
    propensity_coef_sd          = 0.5,
    beta_propensity_intercept_1 = beta_intercept,
    beta_propensity             = if (is.null(beta)) rep(0.0, n_covar) else beta
  )
}

test_stan_propensity <- function(data) {
  test_stan_function(
    here("tests", "testthat", "stan", "test_propensity_weights.stan"),
    data
  )
}

test_that("propensity transformed_data: target range and indicator are correct", {
  data <- make_propensity_data()
  fit  <- test_stan_propensity(data)
  d    <- posterior::as_draws_df(fit$draws())

  expect_equal(get_stan_val(d, "out_target_start"), 1)
  expect_equal(get_stan_val(d, "out_target_end"),   3)
  expect_equal(get_stan_val(d, "out_n_target"),     3)

  # Trial patients (1–3) are target; RWD patients (4–5) are not
  for (i in 1:3) expect_equal(get_stan_val(d, "out_is_target", i), 1,
    label = paste0("out_is_target[", i, "]"))
  for (i in 4:5) expect_equal(get_stan_val(d, "out_is_target", i), 0,
    label = paste0("out_is_target[", i, "]"))
})

test_that("propensity weights: target patients always get weight 1.0", {
  data <- make_propensity_data(beta_intercept = 2.0, beta = c(0.5, -0.5))
  fit  <- test_stan_propensity(data)
  d    <- posterior::as_draws_df(fit$draws())

  # Target patients (1–3): weight must be exactly 1.0 regardless of betas
  for (i in 1:3) expect_equal(get_stan_val(d, "out_likelihood_weight", i), 1.0,
    tolerance = 1e-10, label = paste0("weight[", i, "]"))
})

test_that("propensity weights: non-target patients get inv_logit score", {
  beta_intercept <- 0.5
  beta           <- c(0.3, -0.2)
  data <- make_propensity_data(beta_intercept = beta_intercept, beta = beta)
  fit  <- test_stan_propensity(data)
  d    <- posterior::as_draws_df(fit$draws())

  X <- data$covar_design_matrix

  # RWD patients (4, 5): weight = inv_logit(intercept + X[i,] %*% beta)
  for (i in 4:5) {
    expected <- plogis(beta_intercept + sum(X[i, ] * beta))
    expect_equal(get_stan_val(d, "out_likelihood_weight", i), expected,
      tolerance = 1e-6, label = paste0("weight[", i, "]"))
  }
})

test_that("propensity weights: disabled mode gives all-ones weights", {
  data <- make_propensity_data(enable = 0L, beta_intercept = 99.0, beta = c(5, 5))
  fit  <- test_stan_propensity(data)
  d    <- posterior::as_draws_df(fit$draws())

  for (i in 1:5) expect_equal(get_stan_val(d, "out_likelihood_weight", i), 1.0,
    tolerance = 1e-10, label = paste0("weight[", i, "] when disabled"))
})

test_that("propensity transformed_data: correct when RWD patients come first", {
  # RWD patients first, trial patients last — contiguity still holds
  data <- make_propensity_data()
  # Flip patient ordering: RWD (group 2) first, trial (group 1) last
  data$patient_level_groups <- rbind(
    cbind(rep(2L, 2), 1:2),   # patients 1–2: RWD
    cbind(rep(1L, 3), 3:5)    # patients 3–5: trial
  )
  data$propensity_target_group <- 1L
  fit <- test_stan_propensity(data)
  d   <- posterior::as_draws_df(fit$draws())

  expect_equal(get_stan_val(d, "out_target_start"), 3)
  expect_equal(get_stan_val(d, "out_target_end"),   5)
  expect_equal(get_stan_val(d, "out_n_target"),     3)
  for (i in 1:2) expect_equal(get_stan_val(d, "out_is_target", i), 0,
    label = paste0("out_is_target[", i, "]"))
  for (i in 3:5) expect_equal(get_stan_val(d, "out_is_target", i), 1,
    label = paste0("out_is_target[", i, "]"))
})
```

- [ ] **Step 3: Run the propensity tests**

```bash
cd /mnt/code/.worktrees/karim/propensity-borrowing
Rscript -e 'testthat::test_file("tests/testthat/test-stan-propensity.R")'
```

Expected: 5 tests pass. The first run compiles the Stan test model (~30s); subsequent runs use the cached binary.

- [ ] **Step 4: Commit**

```bash
git add tests/testthat/stan/test_propensity_weights.stan \
        tests/testthat/test-stan-propensity.R
git commit -m "test: add propensity module tests (range detection, weights, disabled mode)"
```

---

## Task 5: Update `pioneer.stan`  <!-- renumbered from 4 -->

Three changes: (1) add module includes, (2) convert PSA `~` to `target +=` with weight, (3) simplify multistate block.

**Files:**
- Modify: `stan/psa/pioneer.stan`

- [ ] **Step 1: Add propensity module includes to each block**

In `data` block (after existing `#include "modules/laplace/data.stan"`):
```stan
  #include "modules/propensity/flags.stan"
  #include "modules/propensity/hyperparams.stan"
```

In `transformed data` block (after `#include "modules/laplace/transformed_data.stan"`):
```stan
  #include "modules/propensity/transformed_data.stan"
```

In `parameters` block (after `#include "modules/multistate/parameters.stan"`):
```stan
  #include "modules/propensity/parameters.stan"
```

In `transformed parameters` block (after `#include "modules/multistate/transformed_parameters.stan"`):
```stan
  #include "modules/propensity/transformed_parameters.stan"
```

In `model` block priors section (after `#include "modules/multistate/priors.stan"`):
```stan
  #include "modules/propensity/priors.stan"
```

- [ ] **Step 2: Convert PSA likelihood to weighted `target +=` (lines 90–110)**

Replace:
```stan
        normalized_psa[data_start + n_screen : data_end] ~ sf_log_space_obs(
          states[state_start + n_screen : state_start + (data_end - data_start)],
          measure_sd_psa,
          negative_infinity(),  // No LOD censoring for PSA
          measure_nu
        );
```

With:
```stan
        target += likelihood_weight[p] * sf_log_space_obs_lpdf(
          normalized_psa[data_start + n_screen : data_end] |
          states[state_start + n_screen : state_start + (data_end - data_start)],
          measure_sd_psa,
          negative_infinity(),  // No LOD censoring for PSA
          measure_nu
        );
```

- [ ] **Step 3: Simplify and weight multistate block (lines 113–142)**

Replace the entire `if (fit_multistate_data)` block with:
```stan
    if (fit_multistate_data) {
      profile("multistate loglik") {
        ms_final_state[hmc_patient_idx] ~ multistate(
          likelihood_weight[hmc_patient_idx],
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

- [ ] **Step 4: Syntax-check `pioneer.stan`**

```bash
cd /mnt/code/.worktrees/karim/propensity-borrowing
~/.cmdstan/cmdstan-2.38.0/bin/stanc --include-paths=stan --include-paths=stan/psa \
  stan/psa/pioneer.stan
```

Expected: `Model name: pioneer_model` with no errors.

- [ ] **Step 5: Commit**

```bash
git add stan/psa/pioneer.stan
git commit -m "feat: add propensity weighting to pioneer.stan PSA and multistate likelihoods"
```

---

## Task 6: R — priors, initializers, and standalone data prep

**Files:**
- Modify: `r/pioneer/priors.R`
- Modify: `r/pioneer/initializers.R`
- Modify: `r/pioneer/prepare_analysis_data.R`

- [ ] **Step 1: Add propensity hyperparameter defaults to `r/pioneer/priors.R`**

In `get_pioneer_priors()`, append to the body of the `lst(...)` call — the line immediately before the closing `) |>` (not inside the `list_assign` pipe chain that follows):

```r
    # =========================================================================
    # Propensity module hyperparams
    # SD scaled to avoid U-shaped prior predictive: logit SD ≈ sqrt(1 + p * coef_sd²)
    # For n_covar=16: coef_sd ≈ 0.5, logit SD ≈ 2.2 (keeps scores in ~0.1–0.9)
    # =========================================================================
    propensity_intercept_sd = 1.0,
    propensity_coef_sd = if (n_covar > 0) 2.0 / sqrt(n_covar) else 1.0,
```

- [ ] **Step 2: Add propensity initializers to `r/pioneer/initializers.R`**

In `create_pioneer_initializer()`, inside the `biomarker_init <- lst(...)` call, append (before the closing `)`):

```r
        # Propensity coefficients — zero init gives weight=0.5 for all RWD patients,
        # a neutral starting point; likelihood quickly moves toward MLE during warmup
        beta_propensity_intercept = if (enable_propensity_weighting) array(0, dim = 1),
        beta_propensity = if (enable_propensity_weighting && n_covar > 0) rep(0, n_covar),
```

- [ ] **Step 3: Add propensity flags to `prepare_ms_standalone_stan_data()` in `r/pioneer/prepare_analysis_data.R`**

The standalone model doesn't use propensity weighting. In `prepare_ms_standalone_stan_data()`, add to the `list_assign(...)` call:

```r
      enable_propensity_weighting = 0L,
      propensity_split_level      = 0L,
      propensity_target_group     = 0L,
      propensity_intercept_sd     = 1.0,
      propensity_coef_sd          = 1.0,
```

- [ ] **Step 4: Run existing tests to confirm no breakage**

```bash
cd /mnt/code/.worktrees/karim/propensity-borrowing
Rscript -e 'testthat::test_dir("tests/testthat")'
```

Expected: all previously passing tests still pass.

- [ ] **Step 5: Commit**

```bash
git add r/pioneer/priors.R \
        r/pioneer/initializers.R \
        r/pioneer/prepare_analysis_data.R
git commit -m "feat: add propensity priors, initializers, and standalone data-prep flags"
```

---

## Task 7: Add `propensity` model variant to targets pipeline

**Files:**
- Modify: `targets/pioneer_targets.R`

- [ ] **Step 1: Add `~propensity` column to the `tar_map` tribble**

Locate the `tribble(` at line ~337. Add a `~propensity` column and a new row:

```r
    tribble(
      ~model,                ~hist,  ~model_formula,      ~trial_level, ~arm_level, ~propensity,
      "combined",            TRUE,   covar_formula,        TRUE,         TRUE,       FALSE,
      "combined_no_trial",   TRUE,   covar_formula,        FALSE,        TRUE,       FALSE,
      "no_hist",             FALSE,  covar_formula,        FALSE,        TRUE,       FALSE,
      "combined_safe_covar", TRUE,   safe_covar_formula,   TRUE,         TRUE,       FALSE,
      "propensity",          TRUE,   covar_formula,        FALSE,        TRUE,       TRUE,
    ),
```

- [ ] **Step 2: Add propensity flags to `default_settings`**

In the `default_settings` target (the `lst(...)` around line ~386), add:

```r
        enable_propensity_weighting = as.integer(propensity),
        propensity_split_level      = if (propensity) 1L else 0L,
        # Mirror the Laplace pattern: read target group dynamically from stan data
        # rather than hardcoding 1L. Trial patients are always the first group at
        # level 1, but this is robust to any future reordering of factor levels.
        propensity_target_group     = if (propensity) {
          base_pioneer_stan_data$patient_level_groups[1, 1]
        } else 0L,
```

The `propensity_target_group` is read dynamically from `patient_level_groups[1, 1]` — the group ID of the first patient (always a trial patient) at level 1. This mirrors the Laplace pattern at line ~434: `base_pioneer_stan_data$patient_level_groups[target_idx[1], 1]`. For non-propensity variants (`propensity = FALSE`), all three fields are `0L`.

- [ ] **Step 3: Verify R parses cleanly**

```bash
cd /mnt/code/.worktrees/karim/propensity-borrowing
Rscript -e 'source("targets/pioneer_targets.R")'
```

Expected: no errors (targets definitions don't execute on source).

- [ ] **Step 4: Commit**

```bash
git add targets/pioneer_targets.R
git commit -m "feat: add propensity model variant to pioneer targets pipeline"
```

---

## Task 8: End-to-end verification

- [ ] **Step 1: Final syntax check of all Stan models**

```bash
cd /mnt/code/.worktrees/karim/propensity-borrowing

~/.cmdstan/cmdstan-2.38.0/bin/stanc --include-paths=stan --include-paths=stan/psa \
  stan/psa/pioneer.stan

~/.cmdstan/cmdstan-2.38.0/bin/stanc --include-paths=stan --include-paths=stan/psa \
  stan/psa/ms-standalone.stan

~/.cmdstan/cmdstan-2.38.0/bin/stanc --include-paths=stan --include-paths=stan/tumor \
  stan/tumor/sf-ssm-log-space.stan

~/.cmdstan/cmdstan-2.38.0/bin/stanc --include-paths=stan --include-paths=stan/tumor \
  stan/tumor/ms-standalone.stan
```

Expected: all four compile cleanly with no errors.

- [ ] **Step 2: Run full test suite**

```bash
cd /mnt/code/.worktrees/karim/propensity-borrowing
Rscript -e 'testthat::test_dir("tests/testthat")'
```

Expected: all tests pass.

- [ ] **Step 3: Verify disabled mode is a no-op**

With `enable_propensity_weighting = 0L`, `likelihood_weight = ones_vector(n_patients)` — the model should behave identically to `combined_no_trial`. Confirm syntactically via `stanc` (full runtime verification is done in the Domino job comparison).

- [ ] **Step 4: Final commit**

```bash
git add -A
git commit -m "feat: propensity-weighted borrowing — full implementation

- New stan/modules/propensity/ module (5 files)
- multistate_lpmf gains per-patient weight vector
- pioneer.stan PSA and multistate likelihoods weighted
- R priors, initializers, targets pipeline updated
- All existing tests pass"
```
