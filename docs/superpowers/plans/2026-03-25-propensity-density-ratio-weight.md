# Propensity Weight: Capped Density Ratio Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the propensity weight formula from raw `P(trial|X)` to a capped density ratio `min(1, P(X|trial)/P(X|RWD))`, and fully disable trial-level hierarchy for the propensity model variant.

**Architecture:** The density ratio `P(X|trial)/P(X|RWD)` equals the stabilized propensity weight: `exp(logit(P(trial|X)) - logit(P(trial)))`. In the identical-twin scenario (same covariates in trial and RWD), this gives weight=1 instead of the current P(trial)~0.09. Capping at 1 prevents any RWD patient from contributing more than a trial patient. The Bernoulli propensity likelihood (priors.stan) is unchanged — only the weight transformation changes.

**Tech Stack:** Stan 2.38, R tidyverse, `targets`, CmdStanR.

---

## Motivation

The current weight `P(trial|X)` conflates covariate similarity with sample-size ratio. When N_rwd >> N_trial, even a perfect covariate match gets weight ~ N_trial/(N_trial+N_rwd) ~ 0.09. The capped density ratio isolates pure covariate similarity:

| Patient profile | Current P(trial\|X) | New min(1, density ratio) |
|----------------|---------------------|---------------------------|
| Perfect match (identical covariates) | ~0.09 | 1.0 |
| Trial-typical but RWD-rare | ~0.5 | 1.0 (capped) |
| RWD-typical but trial-rare | ~0.01 | ~0.1 |
| Completely dissimilar | ~0.001 | ~0.01 |

## File Map

| File | Action | Description |
|------|--------|-------------|
| `stan/modules/propensity/transformed_data.stan` | Modify | Add `propensity_log_marginal_odds` computation |
| `stan/modules/propensity/transformed_parameters.stan` | Modify | Change weight from `inv_logit(...)` to capped density ratio |
| `tests/testthat/stan/test_propensity_weights.stan` | Modify | Update inlined GQ logic to match new formula |
| `tests/testthat/test-stan-propensity.R` | Modify | Update expected values, add density-ratio-specific tests |
| `targets/pioneer_targets.R` | Modify | Use `trial_level` for multistate baseline hazard flag |

---

## Task 1: Add `log_marginal_odds` to propensity transformed_data

The density ratio needs the marginal log-odds `log(N_target / N_nontarget)` computed once at data time.

**Files:**
- Modify: `stan/modules/propensity/transformed_data.stan:36-43`

- [ ] **Step 1: Add `propensity_log_marginal_odds` after existing `propensity_is_target` loop**

After the `propensity_is_target` loop (line 36) and before the print statements, add:

```stan
  // Marginal log-odds: log(N_target / N_nontarget)
  // Used to convert logit(P(trial|X)) → log density ratio
  real propensity_log_marginal_odds = log(propensity_n_target * 1.0)
                                     - log((n_patients - propensity_n_target) * 1.0);
```

Also declare the variable at the top of the file (before the `if` block), with a default value for the disabled case:

```stan
real propensity_log_marginal_odds = 0.0;
```

The full file should read (showing the new lines integrated):

```stan
// modules/propensity/transformed_data.stan
// Locate the contiguous range of target patients using the level-agnostic
// patient_level_groups matrix. Runs once at initialization.
int propensity_target_start = 1;
int propensity_target_end = n_patients;
int propensity_n_target = n_patients;

// Binary indicator: 1 = target patient (weight=1), 0 = non-target (weight=propensity score)
array[n_patients] int<lower=0, upper=1> propensity_is_target = rep_array(1, n_patients);

// Marginal log-odds of being a target patient: log(N_target / N_nontarget).
// Used to stabilize propensity scores into density ratios.
real propensity_log_marginal_odds = 0.0;

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

  // Marginal log-odds: log(N_target / N_nontarget)
  // Stabilizes the propensity logit into a density ratio:
  //   log(P(X|trial)/P(X|RWD)) = logit(P(trial|X)) - logit(P(trial))
  propensity_log_marginal_odds = log(propensity_n_target * 1.0)
                                 - log((n_patients - propensity_n_target) * 1.0);

  print("Propensity weighting enabled:");
  print("  split_level=", propensity_split_level,
        " target_group=", propensity_target_group);
  print("  target: [", propensity_target_start, ", ", propensity_target_end,
        "] n=", propensity_n_target,
        " non-target: ", n_patients - propensity_n_target);
  print("  log_marginal_odds=", propensity_log_marginal_odds);
}
```

- [ ] **Step 2: Verify no syntax errors**

```bash
# Quick check — the module is included via pioneer.stan
~/.cmdstan/cmdstan-2.38.0/bin/stanc --include-paths=stan --include-paths=stan/psa \
  stan/psa/pioneer.stan
```

Expected: compiles cleanly.

- [ ] **Step 3: Commit**

```bash
git add stan/modules/propensity/transformed_data.stan
git commit -m "feat(propensity): compute log_marginal_odds for density ratio stabilization"
```

---

## Task 2: Change weight formula to capped density ratio

Replace `inv_logit(logit_score)` with `min(1, exp(logit_score - log_marginal_odds))`.

**Files:**
- Modify: `stan/modules/propensity/transformed_parameters.stan`

- [ ] **Step 1: Replace the weight computation**

The full file should read:

```stan
// modules/propensity/transformed_parameters.stan
// Compute per-patient likelihood weights as capped density ratios.
//
// The density ratio P(X|trial)/P(X|RWD) measures pure covariate similarity,
// free of the sample-size base rate that makes raw P(trial|X) uniformly low
// when N_rwd >> N_trial.
//
// Math: log(P(X|trial)/P(X|RWD)) = logit(P(trial|X)) - logit(P(trial))
//     = (beta_0 + X*beta) - log(N_target/N_nontarget)
//
// Capping at 1.0 ensures no RWD patient contributes more than a trial patient.
//
// NOTE: uses covar_design_matrix (original centered/scaled matrix) NOT Q_covar_design_matrix.
// All other modules use the QR-decomposed Q matrix for numerical stability, but propensity
// logistic regression is convex and doesn't need QR. Using the original matrix keeps
// beta_propensity coefficients interpretable (log-odds per SD change per covariate).
vector<lower=0, upper=1>[n_patients] likelihood_weight = ones_vector(n_patients);

if (enable_propensity_weighting && propensity_split_level > 0) {
  // Log density ratio for all patients
  vector[n_patients] log_density_ratio =
    beta_propensity_intercept[1] + covar_design_matrix * beta_propensity
    - propensity_log_marginal_odds;

  // Cap at 1.0: min(1, exp(log_dr)) = exp(min(0, log_dr))
  for (i in 1:n_patients) {
    likelihood_weight[i] = exp(fmin(0.0, log_density_ratio[i]));
  }

  // Override target patients to weight 1.0 (range assignment, no loop)
  likelihood_weight[propensity_target_start:propensity_target_end] =
    ones_vector(propensity_n_target);
}
```

- [ ] **Step 2: Verify syntax**

```bash
~/.cmdstan/cmdstan-2.38.0/bin/stanc --include-paths=stan --include-paths=stan/psa \
  stan/psa/pioneer.stan
```

Expected: compiles cleanly.

- [ ] **Step 3: Commit**

```bash
git add stan/modules/propensity/transformed_parameters.stan
git commit -m "feat(propensity): switch weight to capped density ratio min(1, P(X|trial)/P(X|RWD))"
```

---

## Task 3: Update test Stan file and R tests

The test Stan file inlines the `transformed_parameters.stan` logic in GQ. Must update to match the new formula. The R test expected values change accordingly.

**Files:**
- Modify: `tests/testthat/stan/test_propensity_weights.stan:38-54`
- Modify: `tests/testthat/test-stan-propensity.R:74-89`

- [ ] **Step 1: Update test Stan file GQ block**

Replace the `generated quantities` block (lines 38-54) with:

```stan
generated quantities {
  // Outputs from transformed_data.stan
  int out_target_start           = propensity_target_start;
  int out_target_end             = propensity_target_end;
  int out_n_target               = propensity_n_target;
  array[n_patients] int out_is_target = propensity_is_target;
  real out_log_marginal_odds     = propensity_log_marginal_odds;

  // transformed_parameters.stan logic inlined (block includes cannot go in GQ)
  // Capped density ratio: min(1, exp(logit_score - log_marginal_odds))
  vector[n_patients] out_likelihood_weight = ones_vector(n_patients);
  if (enable_propensity_weighting && propensity_split_level > 0) {
    vector[n_patients] log_density_ratio =
      beta_propensity_intercept[1] + covar_design_matrix * beta_propensity
      - propensity_log_marginal_odds;

    for (i in 1:n_patients) {
      out_likelihood_weight[i] = exp(fmin(0.0, log_density_ratio[i]));
    }

    out_likelihood_weight[propensity_target_start:propensity_target_end] =
      ones_vector(propensity_n_target);
  }
}
```

- [ ] **Step 2: Update R test — rename existing test and fix expected values**

In `tests/testthat/test-stan-propensity.R`, replace the test at lines 74-89 ("non-target patients get inv_logit score") with a density-ratio test:

```r
test_that("propensity weights: non-target patients get capped density ratio", {
  beta_intercept <- 0.5
  beta           <- c(0.3, -0.2)
  data <- make_propensity_data(beta_intercept = beta_intercept, beta = beta)
  fit  <- test_stan_propensity(data)
  d    <- posterior::as_draws_df(fit$draws())

  X <- data$covar_design_matrix
  n_trial <- 3L
  n_total <- 5L
  log_marginal_odds <- log(n_trial) - log(n_total - n_trial)

  # RWD patients (4, 5): weight = min(1, exp(logit_score - log_marginal_odds))
  for (i in 4:5) {
    logit_score <- beta_intercept + sum(X[i, ] * beta)
    log_dr <- logit_score - log_marginal_odds
    expected <- min(1.0, exp(log_dr))
    expect_equal(get_stan_val(d, "out_likelihood_weight", i), expected,
      tolerance = 1e-6, label = paste0("weight[", i, "]"))
  }
})
```

- [ ] **Step 3: Add identical-twin test**

This test verifies the key motivating property: when trial and RWD patients have identical covariates, the density ratio is 1.0 (not P(trial)).

```r
test_that("propensity weights: identical covariates give weight=1 (not base rate)", {
  # 3 trial + 2 RWD, all with IDENTICAL covariates
  # Logistic regression with beta=0 should give weight = 1.0 for RWD patients
  # (density ratio = 1 when covariates are non-discriminative)
  data <- make_propensity_data(beta_intercept = 0.0, beta = c(0.0, 0.0))

  # With beta=0: logit_score = intercept = 0.0
  # log_marginal_odds = log(3/2) = 0.405
  # log_density_ratio = 0.0 - 0.405 = -0.405
  # weight = exp(-0.405) = 0.667
  # This is NOT 1.0 because the intercept doesn't match the empirical logit.
  # But if the intercept equals logit(P(trial)) = logit(3/5) = 0.405:
  data2 <- make_propensity_data(
    beta_intercept = log(3/2),  # = logit(P(trial)) when beta=0
    beta = c(0.0, 0.0)
  )
  fit <- test_stan_propensity(data2)
  d   <- posterior::as_draws_df(fit$draws())

  # Now: logit_score = log(3/2), log_marginal_odds = log(3/2)
  # log_density_ratio = 0 → weight = 1.0
  for (i in 4:5) {
    expect_equal(get_stan_val(d, "out_likelihood_weight", i), 1.0,
      tolerance = 1e-6,
      label = paste0("identical-twin weight[", i, "] should be 1.0"))
  }
})
```

- [ ] **Step 4: Add weight-capping test**

Verify that weights > 1 are capped to 1.0:

```r
test_that("propensity weights: density ratio > 1 is capped at 1.0", {
  # Large positive beta makes RWD patients look very trial-like
  # (high logit_score → density ratio >> 1 → capped at 1.0)
  data <- make_propensity_data(beta_intercept = 5.0, beta = c(2.0, 2.0))
  fit  <- test_stan_propensity(data)
  d    <- posterior::as_draws_df(fit$draws())

  # RWD patients should have density ratio >> 1, capped to exactly 1.0
  for (i in 4:5) {
    expect_equal(get_stan_val(d, "out_likelihood_weight", i), 1.0,
      tolerance = 1e-10, label = paste0("capped weight[", i, "]"))
  }
})
```

- [ ] **Step 5: Add `out_log_marginal_odds` check to existing transformed_data test**

In the first test ("propensity transformed_data: target range and indicator are correct"), add:

```r
  # log_marginal_odds = log(3/2)
  expect_equal(get_stan_val(d, "out_log_marginal_odds"), log(3/2), tolerance = 1e-6)
```

- [ ] **Step 6: Run tests**

```bash
Rscript -e 'testthat::test_file("tests/testthat/test-stan-propensity.R")'
```

Expected: 7 tests pass (5 original + 2 new). The first run recompiles the test Stan model (~30s).

- [ ] **Step 7: Commit**

```bash
git add tests/testthat/stan/test_propensity_weights.stan \
        tests/testthat/test-stan-propensity.R
git commit -m "test(propensity): update tests for capped density ratio weight formula"
```

---

## Task 4: Disable trial-level multistate for propensity variant

The propensity variant already disables trial-level intercepts for TR/frac/init modules (`trial_level = FALSE`), but the multistate baseline hazard flag is hardcoded to `trial = TRUE` regardless. This should also respect `trial_level`.

**Files:**
- Modify: `targets/pioneer_targets.R:411`

- [ ] **Step 1: Make `enable_ms_level_baseline_hazard` respect `trial_level`**

Change line 411 from:

```r
enable_ms_level_baseline_hazard = as.array(as.integer(c(trial = TRUE, arm = FALSE, patient = FALSE))),
```

To:

```r
enable_ms_level_baseline_hazard = as.array(as.integer(c(trial = trial_level, arm = FALSE, patient = FALSE))),
```

This means the propensity variant (`trial_level = FALSE`) will have NO trial-level multistate GP — population parameters only, with propensity weights handling covariate alignment.

- [ ] **Step 2: Verify R parses cleanly**

```bash
Rscript -e 'source("targets/pioneer_targets.R")'
```

Expected: no errors.

- [ ] **Step 3: Commit**

```bash
git add targets/pioneer_targets.R
git commit -m "fix(propensity): disable trial-level multistate baseline hazard when trial_level=FALSE"
```

---

## Task 5: End-to-end verification

- [ ] **Step 1: Syntax-check all Stan models**

```bash
~/.cmdstan/cmdstan-2.38.0/bin/stanc --include-paths=stan --include-paths=stan/psa \
  stan/psa/pioneer.stan

~/.cmdstan/cmdstan-2.38.0/bin/stanc --include-paths=stan --include-paths=stan/psa \
  stan/psa/ms-standalone.stan

~/.cmdstan/cmdstan-2.38.0/bin/stanc --include-paths=stan --include-paths=stan/tumor \
  stan/tumor/sf-ssm-log-space.stan
```

Expected: all compile cleanly.

- [ ] **Step 2: Run full test suite**

```bash
Rscript -e 'testthat::test_dir("tests/testthat")'
```

Expected: all tests pass.

- [ ] **Step 3: Verify disabled mode is still a no-op**

When `enable_propensity_weighting = 0`, the if-block is skipped and `likelihood_weight` stays `ones_vector(n_patients)`. The disabled-mode test from Task 3 confirms this. No additional verification needed.

- [ ] **Step 4: Squash or final commit if needed**

```bash
git log --oneline -5  # Review commit history
# Optionally squash if desired
```

---

## Summary of Mathematical Change

**Before (raw propensity):**
```
w_i = inv_logit(beta_0 + X_i * beta)  = P(trial | X_i)
```

**After (capped density ratio):**
```
w_i = min(1, exp(beta_0 + X_i * beta - log(N_trial / N_rwd)))
    = min(1, P(X_i | trial) / P(X_i | RWD))
```

The Bernoulli likelihood `propensity_is_target ~ bernoulli_logit(beta_0 + X*beta)` is unchanged. It still identifies the propensity coefficients. Only the transformation from coefficients to weights changes.

## Impact on Other Model Variants

- **`combined`** (`propensity=FALSE`): No change. `enable_propensity_weighting=0`, all weights are 1.0.
- **`combined_no_trial`** (`propensity=FALSE`): No change.
- **`no_hist`** (`propensity=FALSE`): No change.
- **`combined_safe_covar`** (`propensity=FALSE`): No change.
- **`propensity`** (`propensity=TRUE`): Weight formula changes + multistate trial-level disabled.

The `enable_ms_level_baseline_hazard` change (Task 4) affects other variants only if they set `trial_level=FALSE`: that's `combined_no_trial` and `no_hist`. Both already had `trial_level=FALSE` for TR/frac/init, so making multistate consistent is correct behavior.
