# Visit-Gated 0->1 Transition Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Decouple the PSA state-space model from the multistate 0->1 hazard by making progression detection visit-gated (hazard only fires at assessment visits using observed PSA), while keeping 0->2 optionally coupled through the modeled PSA trajectory.

**Architecture:** A new flag `enable_ms_visit_gated_01` controls whether the 0->1 transition uses visit-gated likelihood (new) or continuous weekly hazard (current). When visit-gated: (1) the 0->1 likelihood accumulates only at assessment weeks (like 0->3 dropout), (2) the time-varying covariate for 0->1 is built from observed PSA in `transformed data` (not modeled PSA in `transformed parameters`), (3) interval censoring for 0->1 is disabled. A second flag `enable_ms_02_time_varying_cov` independently controls whether 0->2 uses the modeled PSA time-varying covariate or GP baseline + time-invariant covariates only.

**Tech Stack:** Stan (modules), R (targets pipeline data prep, priors, initializers), testthat

---

## File Map

### Stan files to modify
- `stan/modules/multistate/flags.stan` — add `enable_ms_visit_gated_01` flag
- `stan/multistate.stanfunctions` — add `sum_at_visits_below()` generic, modify `multistate_lpmf` for visit-gated 0->1 path
- `stan/modules/multistate/transformed_parameters.stan` — guard 0->1 time-varying covariate block, add observed-PSA covariate path, declare `ms_observed_covar_01` placeholder for all models
- `stan/modules/multistate/data.stan` — add per-transition covariate counts
- `stan/modules/multistate/priors.stan` — update guards from `n_time_varying_covar` to per-transition counts

### Stan files to create
- `stan/psa/_psa_observed_covar_transformed_data.stan` — compute observed-PSA covariate matrix at visit times, expanded to weekly grid via carry-forward

### Stan files with minor edits
- `stan/psa/pioneer.stan` — add include for new transformed data file
- `stan/psa/ms-standalone.stan` — set flag to 0 (standalone has no PSA data)
- `stan/tumor/sf-ssm-log-space.stan` — pass new flag to `multistate_lpmf`
- `stan/tumor/ms-standalone.stan` — pass new flag to `multistate_lpmf`

### R files to modify
- `r/pioneer/prepare_analysis_data.R` — pass new flags, set per-transition covariate counts
- `r/priors.R` — update `get_multistate_priors()` to accept per-transition covariate counts
- `r/pioneer/priors.R` — pass per-transition counts from `stan_data`
- `r/initializers_ms.R` — use per-transition covariate counts
- `r/initializers_fixed.R` — use per-transition covariate counts (also references `n_time_varying_covar` for `time_varying_coef_01/02`)
- `r/pioneer/prepare_analysis_data.R` (`prepare_ms_standalone_stan_data`) — pass new flag = 0

### R files to create
- `r/pioneer/build_observed_psa_covar.R` — helper to build observed-PSA covariate arrays from visit data

### Test files to modify
- `tests/testthat/test-pioneer-multistate.R` — add tests for visit-gated flag propagation

### Test files to create
- `tests/testthat/test-observed-psa-covar.R` — unit tests for observed PSA covariate construction

---

## Design Decisions

### Flag semantics

```
enable_ms_visit_gated_01 = 1:
  - 0->1 likelihood: visit-gated (sum at visit weeks only, like 0->3)
  - 0->1 time-varying covariate: observed PSA (data, no gradient coupling)
  - 0->1 interval censoring: disabled (ms_ic_gap_01 ignored)
  - 0->2 time-varying covariate: controlled by separate flag

enable_ms_visit_gated_01 = 0:
  - Current behavior: continuous weekly hazard with modeled PSA trajectory
  - Full backward compatibility

enable_02_time_varying_cov (R-only argument, not a Stan data field):
  TRUE  -> n_time_varying_covar_02 = n_time_varying_covar (modeled PSA for 0->2)
  FALSE -> n_time_varying_covar_02 = 0 (GP baseline + time-invariant covariates only)
  Controls gradient coupling between PSA state-space and 0->2 hazard.
  When FALSE with visit_gated_01=TRUE: zero coupling between PSA and multistate.
```

### Retained data fields

`n_time_varying_covar` is retained as a Stan data field. It represents the full covariate
count (3: PSA level, decrease rate, growth rate). The `_12` transition continues to size
on `n_time_varying_covar` (currently unsupported, stays 0). The new per-transition counts
`n_time_varying_covar_01` and `n_time_varying_covar_02` are independent data fields set by R.

### Observed-PSA covariate construction (transformed data)

For visit-gated mode, build `ms_observed_covar_01[n_time_varying_covar_01][n_patients, max_all_t]` in transformed data:
- At each patient's visit weeks: use `normalized_psa[v]` -> `log(baseline) + log(normalized)` -> standardize
- Between visits: carry forward the last observed value (step function)
- This matrix is only used for the 0->1 hazard at visit weeks, so the between-visit values don't affect the likelihood — they exist only to keep the matrix rectangular

We only provide Feature 1 (PSA level). Features 2-3 (rates) are dropped in visit-gated mode. This means `n_time_varying_covar` may differ between 0->1 and 0->2. To handle this cleanly:
- In visit-gated mode, 0->1 gets its own covariate count: `n_time_varying_covar_01 = 1` (observed PSA only)
- 0->2 continues using `n_time_varying_covar` (0 or 3, depending on `enable_ms_02_time_varying_cov`)
- The Stan parameter `time_varying_coef_01` is sized by `n_time_varying_covar_01`, not `n_time_varying_covar`

### Splitting covariate counts: 0->1 vs 0->2

Currently `time_varying_coef_01` and `time_varying_coef_02` are both sized `vector[n_time_varying_covar]`. In the new design:
- `time_varying_coef_01` is sized by `n_time_varying_covar_01` (1 in visit-gated mode, `n_time_varying_covar` otherwise)
- `time_varying_coef_02` is sized by `n_time_varying_covar_02` (0 or `n_time_varying_covar` depending on `enable_ms_02_time_varying_cov`)
- Hyperparameters (`time_varying_coef_01_mean`, `_sd`) must match their respective sizes

This is the cleanest approach: each transition has its own covariate dimension, and the R data prep sets them independently.

---

## Tasks

### Task 1: Add flags to Stan and R data prep

**Files:**
- Modify: `stan/modules/multistate/flags.stan`
- Modify: `r/pioneer/prepare_analysis_data.R` (around line 520)

- [ ] **Step 1: Add Stan flags**

In `stan/modules/multistate/flags.stan`, add after the existing covariate flags (after line 34):

```stan
// --- Visit-Gated 0->1 Mode ---
// When enabled, 0->1 hazard accumulates only at assessment visit weeks
// (not every calendar week). Time-varying covariates for 0->1 are built
// from observed PSA (data) rather than modeled trajectory (parameters).
// Interval censoring for 0->1 is disabled in this mode.
int<lower=0, upper=1> enable_ms_visit_gated_01;
```

- [ ] **Step 2: Add per-transition covariate counts to Stan data**

In `stan/modules/multistate/data.stan`, find where `n_time_varying_covar` is used. We need per-transition counts. Add new data declarations:

```stan
// Per-transition time-varying covariate counts
// In visit-gated mode, 0->1 uses observed PSA only (n=1),
// while 0->2 may use modeled trajectory (n=0 or n_time_varying_covar).
int<lower=0> n_time_varying_covar_01;
int<lower=0> n_time_varying_covar_02;
```

- [ ] **Step 3: Update R data prep to pass new flags**

In `r/pioneer/prepare_analysis_data.R`, in the stan data list (around line 520), add:

```r
enable_ms_visit_gated_01 = as.integer(visit_gated_01),
n_time_varying_covar_01  = if (visit_gated_01) 1L else n_time_varying_covar,
n_time_varying_covar_02  = if (enable_02_time_varying_cov) n_time_varying_covar else 0L,
```

Where `visit_gated_01` and `enable_02_time_varying_cov` are new arguments to the data prep function (default `FALSE` for backward compatibility).

- [ ] **Step 4: Update standalone data prep**

In `prepare_ms_standalone_stan_data()` (around line 750), add:

```r
enable_ms_visit_gated_01 = 0L,
n_time_varying_covar_01  = 0L,
n_time_varying_covar_02  = 0L,
```

- [ ] **Step 5: Syntax-check the Stan model**

Run: `~/.cmdstan/cmdstan-2.38.0/bin/stanc --include-paths=stan --include-paths=stan/psa stan/psa/pioneer.stan`
Expected: Clean compile (new data fields declared, not yet consumed differently)

- [ ] **Step 6: Commit**

```bash
git add stan/modules/multistate/flags.stan stan/modules/multistate/data.stan \
  r/pioneer/prepare_analysis_data.R
git commit -m "Add visit-gated 0->1 flag and per-transition covariate counts"
```

---

### Task 2: Split parameter sizing by transition

**Files:**
- Modify: `stan/modules/multistate/parameters.stan` (lines 24, 52)
- Modify: `stan/modules/multistate/transformed_parameters.stan` (lines 79-85, 197-205)
- Modify: `stan/modules/multistate/priors.stan`
- Modify: `stan/modules/multistate/hyperparams.stan`
- Modify: `r/priors.R` (`get_multistate_priors`)
- Modify: `r/pioneer/priors.R`
- Modify: `r/pioneer/initializers.R` or `r/initializers_ms.R`

- [ ] **Step 1: Read current parameter declarations**

Read `stan/modules/multistate/parameters.stan` lines 24 and 52 to confirm current sizing:
- `time_varying_coef_01` is `vector[enable_ms_01 && enable_ms_pop_time_varying_cov ? n_time_varying_covar : 0]`
- `time_varying_coef_02` is `vector[enable_ms_02 && enable_ms_pop_time_varying_cov ? n_time_varying_covar : 0]`

- [ ] **Step 2: Update parameter sizing**

In `stan/modules/multistate/parameters.stan`, change:

```stan
// Line 24: was n_time_varying_covar, now per-transition
vector[enable_ms_01 && enable_ms_pop_time_varying_cov ? n_time_varying_covar_01 : 0] time_varying_coef_01;

// Line 52: was n_time_varying_covar, now per-transition
vector[enable_ms_02 && enable_ms_pop_time_varying_cov ? n_time_varying_covar_02 : 0] time_varying_coef_02;
```

Leave `time_varying_coef_12` unchanged (already unsupported, stays at 0).

- [ ] **Step 3: Read and update hyperparams**

In `stan/modules/multistate/hyperparams.stan`, find `time_varying_coef_*_mean` and `time_varying_coef_*_sd` declarations. Update their sizes to use per-transition counts:

```stan
vector[n_time_varying_covar_01] time_varying_coef_01_mean;
vector[n_time_varying_covar_01] time_varying_coef_01_sd;
vector[n_time_varying_covar_02] time_varying_coef_02_mean;
vector[n_time_varying_covar_02] time_varying_coef_02_sd;
```

- [ ] **Step 4: Update priors guard conditions**

In `stan/modules/multistate/priors.stan`, the existing guard `if (enable_ms_pop_time_varying_cov && n_time_varying_covar > 0)` wraps all time-varying covariate priors. This must be split into per-transition guards:

```stan
// 0->1 time-varying covariate priors
if (enable_ms_pop_time_varying_cov && n_time_varying_covar_01 > 0) {
  time_varying_coef_01 ~ normal(time_varying_coef_01_mean, time_varying_coef_01_sd);
}
// 0->2 time-varying covariate priors
if (enable_ms_pop_time_varying_cov && n_time_varying_covar_02 > 0) {
  time_varying_coef_02 ~ normal(time_varying_coef_02_mean, time_varying_coef_02_sd);
}
// 1->2 time-varying covariate priors (currently unsupported, n_time_varying_covar used)
if (enable_ms_pop_time_varying_cov && n_time_varying_covar > 0) {
  time_varying_coef_12 ~ normal(time_varying_coef_12_mean, time_varying_coef_12_sd);
}
```

Without this, when `n_time_varying_covar > 0` but `n_time_varying_covar_02 == 0`, the code would try to set a prior on a zero-length parameter.

- [ ] **Step 5: Update R priors**

In `r/priors.R` `get_multistate_priors()`, the prior vectors are sized by `n_time_varying_covar`. Update to accept per-transition counts:

```r
get_multistate_priors <- function(n_levels,
                                   n_time_varying_covar,
                                   n_time_invariant_covar,
                                   n_time_varying_covar_01 = n_time_varying_covar,
                                   n_time_varying_covar_02 = n_time_varying_covar) {
  # ...
  time_varying_coef_01_mean = rep(0, n_time_varying_covar_01),
  time_varying_coef_01_sd = rep(0.5, n_time_varying_covar_01),
  time_varying_coef_02_mean = rep(0, n_time_varying_covar_02),
  time_varying_coef_02_sd = rep(0.5, n_time_varying_covar_02),
  # ...
}
```

Update the call site in `r/pioneer/priors.R` to pass the new arguments from `stan_data`.

- [ ] **Step 6: Update R initializer**

In `r/initializers_ms.R`, the initializer for `time_varying_coef_01` uses `n_time_varying_covar`. Update to use `n_time_varying_covar_01`:

```r
time_varying_coef_01 = if (enable_ms_01 && enable_ms_pop_time_varying_cov && n_time_varying_covar_01 > 0) {
  rnorm(n_time_varying_covar_01, time_varying_coef_01_mean, time_varying_coef_01_sd * qr_init_scale)
},
time_varying_coef_02 = if (enable_ms_02 && enable_ms_pop_time_varying_cov && n_time_varying_covar_02 > 0) {
  rnorm(n_time_varying_covar_02, time_varying_coef_02_mean, time_varying_coef_02_sd * qr_init_scale)
},
```

- [ ] **Step 7: Update `r/initializers_fixed.R`**

This file also references `n_time_varying_covar` for `time_varying_coef_01/02` initialization. Update to use per-transition counts:

```r
time_varying_coef_01 = if (enable_ms_01 && enable_ms_pop_time_varying_cov && n_time_varying_covar_01 > 0) {
  rnorm(n_time_varying_covar_01, time_varying_coef_01_mean, time_varying_coef_01_sd * qr_init_scale)
},
time_varying_coef_02 = if (enable_ms_02 && enable_ms_pop_time_varying_cov && n_time_varying_covar_02 > 0) {
  rnorm(n_time_varying_covar_02, time_varying_coef_02_mean, time_varying_coef_02_sd * qr_init_scale)
},
```

- [ ] **Step 8: Verify backward compatibility of `get_tumor_priors()`**

`r/priors.R` line ~417: `get_tumor_priors()` calls `get_multistate_priors(n_levels, n_time_varying_covar, n_time_invariant_covar)`. The new default arguments (`n_time_varying_covar_01 = n_time_varying_covar`, `n_time_varying_covar_02 = n_time_varying_covar`) ensure the tumor path uses the old sizing. No code change needed — just verify.

- [ ] **Step 9: Syntax-check**

Run: `~/.cmdstan/cmdstan-2.38.0/bin/stanc --include-paths=stan --include-paths=stan/psa stan/psa/pioneer.stan`
Expected: Clean compile

- [ ] **Step 10: Commit**

```bash
git add stan/modules/multistate/parameters.stan stan/modules/multistate/hyperparams.stan \
  stan/modules/multistate/priors.stan r/priors.R r/pioneer/priors.R \
  r/initializers_ms.R r/initializers_fixed.R
git commit -m "Split time-varying covariate sizing per transition (01 vs 02)"
```

---

### Task 3: Build observed-PSA covariate in transformed data

**Files:**
- Create: `stan/psa/_psa_observed_covar_transformed_data.stan`
- Modify: `stan/psa/pioneer.stan` (add include)

- [ ] **Step 1: Create the observed-PSA covariate file**

Create `stan/psa/_psa_observed_covar_transformed_data.stan`:

```stan
// ============================================================================
// OBSERVED-PSA COVARIATE FOR VISIT-GATED 0->1 TRANSITION
// ============================================================================
// Populates ms_observed_covar_01 (declared in modules/multistate/transformed_data.stan)
// from OBSERVED PSA values (data), not from the modeled state-space trajectory.
//
// Between visits, carries forward the last observed value (step function).
// The between-visit values don't affect the likelihood (visit-gated hazard
// only evaluates at visit weeks), but the matrix must be rectangular.
//
// MUST be included AFTER: modules/psa/transformed_data.stan
// (needs normalized_psa, log_baseline_psa, median_log_psa_obs, iqr_log_psa_obs)
// MUST be included AFTER: modules/multistate/transformed_data.stan
// (which declares ms_observed_covar_01)

if (enable_ms_visit_gated_01 && n_time_varying_covar_01 > 0) {
  // Feature 1: Standardized log(observed PSA) — carry-forward between visits
  for (i in 1:n_patients) {
    int visit_start, visit_end;
    (visit_start, visit_end) = get_pos(patient_visit_pos, i);

    // Initialize with baseline PSA
    real last_log_psa = log_baseline_psa[i];

    for (t in 1:max_all_t) {
      // Check if there's a visit at absolute week t
      for (v in visit_start:visit_end) {
        if (t_patient_visits[v] == t && psa_measured[v]) {
          last_log_psa = log_psa_values[v];
        }
      }
      ms_observed_covar_01[1][i, t] = (last_log_psa - median_log_psa_obs) / iqr_log_psa_obs;
    }
  }
}
```

- [ ] **Step 2: Declare `ms_observed_covar_01` in `modules/multistate/transformed_data.stan`**

Add to `stan/modules/multistate/transformed_data.stan` (at the end, after existing declarations):

```stan
// Visit-gated 0->1 observed covariate matrix.
// Declared here so ALL models (PSA, tumor, standalone) have the symbol in scope.
// Zero-sized when visit-gated mode is off — never accessed.
// Populated by the PSA-specific _psa_observed_covar_transformed_data.stan include.
array[enable_ms_visit_gated_01 && n_time_varying_covar_01 > 0 ? n_time_varying_covar_01 : 0]
  matrix[n_patients, max_all_t] ms_observed_covar_01;
```

This ensures the symbol exists in all four Stan models that include `modules/multistate/transformed_parameters.stan`, even when `enable_ms_visit_gated_01 = 0` (zero-sized array, never accessed). The tumor models and standalones pass `enable_ms_visit_gated_01 = 0`, so the array is empty.

- [ ] **Step 3: Add include in pioneer.stan**

In `stan/psa/pioneer.stan`, in the `transformed data` block (after `modules/multistate/transformed_data.stan`), add:

```stan
#include "_psa_observed_covar_transformed_data.stan"
```

This must come after both `modules/psa/transformed_data.stan` (which computes `log_psa_values`, `log_baseline_psa`, `median_log_psa_obs`, `iqr_log_psa_obs`) and `modules/multistate/transformed_data.stan` (which declares `ms_observed_covar_01`).

- [ ] **Step 3: Syntax-check**

Run: `~/.cmdstan/cmdstan-2.38.0/bin/stanc --include-paths=stan --include-paths=stan/psa stan/psa/pioneer.stan`
Expected: Clean compile

- [ ] **Step 4: Commit**

```bash
git add stan/psa/_psa_observed_covar_transformed_data.stan stan/psa/pioneer.stan
git commit -m "Build observed-PSA covariate matrix in transformed data for visit-gated mode"
```

---

### Task 4: Route 0->1 covariate consumption by flag

**Files:**
- Modify: `stan/modules/multistate/transformed_parameters.stan` (lines 76-85)
- Modify: `stan/psa/_psa_ms_time_varying_covar.stan` (guard with flag)

- [ ] **Step 1: Guard the modeled-PSA covariate construction**

In `stan/psa/_psa_ms_time_varying_covar.stan`, the entire file builds `ms_time_varying_covar_01` from `states_full_grid`. Wrap the 0->1 portion with a guard so it's only built when NOT in visit-gated mode (or when 0->2 still needs it):

Change the outer guard (line 21) from:
```stan
if (enable_ms_01 && enable_ms_pop_time_varying_cov && n_time_varying_covar > 0) {
```
to:
```stan
if (enable_ms_pop_time_varying_cov && n_time_varying_covar > 0 &&
    (!enable_ms_visit_gated_01 || n_time_varying_covar_02 > 0)) {
```

And the matrix declaration (line 18) needs to match: only allocate when needed for the continuous path (either non-visit-gated 0->1, or 0->2 with time-varying covariates).

The matrix size stays `n_time_varying_covar` (the full 3 features) since 0->2 may still need all 3. When visit-gated and 0->2 has no time-varying covars, the matrix is zero-sized.

- [ ] **Step 2: Route 0->1 covariate in transformed_parameters.stan**

In `stan/modules/multistate/transformed_parameters.stan`, replace the 0->1 time-varying block (lines 79-85):

```stan
// Time-varying covariate effects for 0->1
if (enable_ms_pop_time_varying_cov) {
  if (enable_ms_visit_gated_01 && n_time_varying_covar_01 > 0) {
    // Visit-gated mode: use observed-PSA covariate (built in transformed data)
    for (k in 1:n_time_varying_covar_01) {
      log_cond_surv_01 += time_varying_coef_01[k] * ms_observed_covar_01[k][hmc_patient_idx, :];
    }
  } else if (n_time_varying_covar_01 > 0) {
    // Continuous mode: use modeled-PSA covariate (built in transformed parameters)
    for (k in 1:n_time_varying_covar_01) {
      log_cond_surv_01 += time_varying_coef_01[k] * ms_time_varying_covar_01[k];
    }
  }
}
```

Note: `ms_observed_covar_01` is indexed by unified patient `[n_patients, max_all_t]` so we subset with `hmc_patient_idx`. `ms_time_varying_covar_01` is already HMC-local `[n_hmc_patients, max_all_t]`.

- [ ] **Step 3: Guard 0->2 covariate block**

In `stan/modules/multistate/transformed_parameters.stan`, update the 0->2 time-varying block (lines 197-205) to use `n_time_varying_covar_02`:

```stan
if (enable_ms_pop_time_varying_cov && n_time_varying_covar_02 > 0) {
  for (k in 1:n_time_varying_covar_02) {
    log_cond_surv_02 += time_varying_coef_02[k] * ms_time_varying_covar_01[k];
  }
}
```

Remove the `&& enable_ms_01` guard — 0->2 covariate availability is now controlled by `n_time_varying_covar_02`.

- [ ] **Step 4: Syntax-check**

Run: `~/.cmdstan/cmdstan-2.38.0/bin/stanc --include-paths=stan --include-paths=stan/psa stan/psa/pioneer.stan`
Expected: Clean compile

- [ ] **Step 5: Commit**

```bash
git add stan/modules/multistate/transformed_parameters.stan \
  stan/psa/_psa_ms_time_varying_covar.stan
git commit -m "Route 0->1 covariate by visit-gated flag: observed vs modeled PSA"
```

---

### Task 5: Visit-gated 0->1 likelihood

**Files:**
- Modify: `stan/multistate.stanfunctions` — add generic `sum_at_visits_below()`, modify `multistate_lpmf`

- [ ] **Step 1: Read `sum_03_at_visits_below` for reference**

The existing function (lines 82-99) sums `log_cond_surv_03[i, wk]` at visit weeks below a limit. We need the same pattern for 0->1.

- [ ] **Step 2: Rename to generic `sum_at_visits_below()`**

Replace `sum_03_at_visits_below` with a generic version that takes the matrix as an argument:

```stan
/**
 * Sum log conditional survival at assessment visit weeks below a limit.
 *
 * Used for visit-gated transitions where hazard accumulates only at
 * assessment visits (e.g., 0->3 dropout, visit-gated 0->1 progression).
 */
real sum_at_visits_below(
    matrix log_cond_surv,
    int i,
    array[] int t_patient_visits,
    array[] int patient_visit_pos,
    int limit
) {
  real ll = 0;
  int v_start; int v_end;
  (v_start, v_end) = get_pos(patient_visit_pos, i);
  for (v in v_start:v_end) {
    int wk = t_patient_visits[v];
    if (wk < 1) continue;
    if (wk >= limit) break;
    ll += log_cond_surv[i, wk];
  }
  return ll;
}
```

Keep `sum_03_at_visits_below` as a thin wrapper that calls the generic version (to avoid breaking all call sites at once), or do a find-and-replace of existing calls.

- [ ] **Step 3: Add `enable_ms_visit_gated_01` parameter to `multistate_lpmf`**

Add `int enable_ms_visit_gated_01` to the function signature (after `enable_ms_32`).

- [ ] **Step 4: Modify state 0 (censored) handling**

In `multistate_lpmf`, for `final_state[i] == 0` (lines 172-194), add a branch:

```stan
// 0->1 survival contribution
if (enable_01 && time_01[i] > 0) {
  if (enable_ms_visit_gated_01) {
    patient_ll += sum_at_visits_below(
      log_cond_surv_01, i, t_patient_visits, patient_visit_pos,
      time_01[i] + 1);
  } else {
    patient_ll += sum(log_cond_surv_01[i, 1:time_01[i]]);
  }
}
```

- [ ] **Step 5: Modify state 1 (progressed) handling**

For `final_state[i] == 1`, the 0->1 survival and event contributions need visit-gating. In the no-IC path (gap == 0 or deterministic), replace:

```stan
// Survival up to T_c
if (enable_01 && T_c > 0) {
  if (enable_ms_visit_gated_01) {
    patient_ll += sum_at_visits_below(
      log_cond_surv_01, i, t_patient_visits, patient_visit_pos, T_c + 1);
  } else {
    patient_ll += sum(log_cond_surv_01[i, 1:T_c]);
  }
}

// Survival from T_c+1 to T_d-1
if (enable_01 && T_d > T_c + 1) {
  if (enable_ms_visit_gated_01) {
    // In visit-gated mode, visits between T_c and T_d are already covered above
    // (sum_at_visits_below with limit = T_d covers them)
    // Nothing extra needed here
  } else {
    patient_ll += sum(log_cond_surv_01[i, (T_c + 1):(T_d - 1)]);
  }
}
```

Actually, in visit-gated mode the survival and event logic simplifies greatly. For a progressed patient:
- Survival: sum log_cond_surv_01 at all visits before T_d
- Event: `log1m_exp(log_cond_surv_01[i, T_d])` at the detection visit

So for visit-gated mode, the entire state-1 0->1 path collapses to:

```stan
if (enable_ms_visit_gated_01) {
  // Visit-gated: no IC, no gap logic. Survival at visits before detection,
  // event at detection visit.
  if (enable_01) {
    patient_ll += sum_at_visits_below(
      log_cond_surv_01, i, t_patient_visits, patient_visit_pos, T_d);
    if (!prog_deterministic[i] && T_d > 0)
      patient_ll += log1m_exp(log_cond_surv_01[i, T_d]);
  }
  // 0->2 still continuous
  if (enable_02 && T_d > 0)
    patient_ll += sum(log_cond_surv_02[i, 1:(T_d - 1)]);
  // 0->3 visit-gated (existing)
  if (enable_03)
    patient_ll += sum_at_visits_below(
      log_cond_surv_03, i, t_patient_visits, patient_visit_pos, T_d);
  // Post-progression sojourn (unchanged)
  // ... existing 1->2 code ...
}
```

This should be a clean early-return branch before the existing IC/no-IC paths.

- [ ] **Step 6: Handle states 2, 3 similarly**

For state 2 (died on trial) and state 3 (dropout), the 0->1 survival contributions follow the same pattern — replace `sum(log_cond_surv_01[i, 1:T])` with `sum_at_visits_below(...)` when visit-gated.

- [ ] **Step 7: Update all call sites**

Update the calls to `multistate_lpmf` in:
- `stan/psa/pioneer.stan` (line 122)
- `stan/psa/ms-standalone.stan` (line 88)
- `stan/tumor/sf-ssm-log-space.stan` (line ~80)
- `stan/tumor/ms-standalone.stan` (line ~115)

Pass `enable_ms_visit_gated_01` as the new argument. For tumor models and standalones, pass `0` (current behavior).

- [ ] **Step 8: Syntax-check both models**

Run:
```bash
~/.cmdstan/cmdstan-2.38.0/bin/stanc --include-paths=stan --include-paths=stan/psa stan/psa/pioneer.stan
~/.cmdstan/cmdstan-2.38.0/bin/stanc --include-paths=stan --include-paths=stan/psa stan/psa/ms-standalone.stan
~/.cmdstan/cmdstan-2.38.0/bin/stanc --include-paths=stan --include-paths=stan/tumor stan/tumor/sf-ssm-log-space.stan
```
Expected: All clean

- [ ] **Step 9: Commit**

```bash
git add stan/multistate.stanfunctions stan/psa/pioneer.stan \
  stan/psa/ms-standalone.stan stan/tumor/sf-ssm-log-space.stan \
  stan/tumor/ms-standalone.stan
git commit -m "Visit-gated 0->1 likelihood: hazard fires only at assessment visits"
```

---

### Task 6: R helper for observed-PSA covariate validation

**Files:**
- Create: `r/pioneer/build_observed_psa_covar.R`
- Create: `tests/testthat/test-observed-psa-covar.R`

- [ ] **Step 1: Write failing tests**

Create `tests/testthat/test-observed-psa-covar.R`:

```r
test_that("build_observed_psa_covar returns correct structure", {
  # Minimal mock data: 2 patients, known PSA values at known visit weeks
  visit_data <- list(
    tibble(week = c(-2, 0, 4, 8), psa = c(10, 10, 5, 20), measured = c(1, 1, 1, 1)),
    tibble(week = c(-1, 0, 6), psa = c(50, 50, 25), measured = c(1, 1, 1))
  )
  result <- build_observed_psa_covar(
    visit_data = visit_data,
    baseline_psa = c(10, 50),
    max_all_t = 10,
    median_log_psa = log(30),
    iqr_log_psa = 1.0
  )
  expect_equal(dim(result), c(2, 10))
  # At week 4 for patient 1: log(5) standardized
  expect_equal(result[1, 4], (log(5) - log(30)) / 1.0)
  # At week 5 for patient 1: carry-forward from week 4
  expect_equal(result[1, 5], result[1, 4])
  # At week 8 for patient 1: new observation
  expect_equal(result[1, 8], (log(20) - log(30)) / 1.0)
})

test_that("build_observed_psa_covar handles missing measurements", {
  visit_data <- list(
    tibble(week = c(0, 4, 8), psa = c(10, NA, 20), measured = c(1, 0, 1))
  )
  result <- build_observed_psa_covar(
    visit_data = visit_data,
    baseline_psa = c(10),
    max_all_t = 10,
    median_log_psa = log(10),
    iqr_log_psa = 1.0
  )
  # Week 4 has no measurement — carry forward baseline
  expect_equal(result[1, 4], (log(10) - log(10)) / 1.0)
  # Week 8 has new observation
  expect_equal(result[1, 8], (log(20) - log(10)) / 1.0)
})
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `Rscript -e 'testthat::test_file("tests/testthat/test-observed-psa-covar.R")'`
Expected: FAIL (function not found)

- [ ] **Step 3: Implement the helper**

Create `r/pioneer/build_observed_psa_covar.R`:

```r
#' Build observed-PSA covariate matrix for visit-gated 0->1 mode
#'
#' Constructs a [n_patients x max_all_t] matrix of standardized log(PSA)
#' using observed values at visit times, with carry-forward between visits.
#'
#' @param visit_data List of per-patient tibbles with columns: week, psa, measured
#' @param baseline_psa Numeric vector of baseline PSA per patient
#' @param max_all_t Integer, maximum calendar week
#' @param median_log_psa Scalar, median of observed log(PSA) for standardization
#' @param iqr_log_psa Scalar, IQR of observed log(PSA) for standardization
#' @return Matrix [n_patients, max_all_t]
build_observed_psa_covar <- function(visit_data, baseline_psa, max_all_t,
                                      median_log_psa, iqr_log_psa) {
  n_patients <- length(visit_data)
  mat <- matrix(NA_real_, nrow = n_patients, ncol = max_all_t)

  for (i in seq_len(n_patients)) {
    vd <- visit_data[[i]]
    last_log_psa <- log(baseline_psa[i])

    for (t in seq_len(max_all_t)) {
      # Check if there's a measured visit at this week
      visit_match <- vd |>
        dplyr::filter(week == t, measured == 1)
      if (nrow(visit_match) > 0) {
        last_log_psa <- log(visit_match$psa[1])
      }
      mat[i, t] <- (last_log_psa - median_log_psa) / iqr_log_psa
    }
  }
  mat
}
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `Rscript -e 'testthat::test_file("tests/testthat/test-observed-psa-covar.R")'`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add r/pioneer/build_observed_psa_covar.R tests/testthat/test-observed-psa-covar.R
git commit -m "Add R helper for observed-PSA covariate construction with tests"
```

---

### Task 7: Update R data prep to pass flag from pipeline

**Files:**
- Modify: `r/pioneer/prepare_analysis_data.R`
- Modify: `targets/pioneer_targets.R` (if flags are set in pipeline)

- [ ] **Step 1: Add function arguments**

In `prepare_psa_stan_data()` (or whichever function assembles the full stan data list), add arguments:

```r
visit_gated_01 = FALSE,
enable_02_time_varying_cov = TRUE
```

Wire them through to the stan data list as designed in Task 1.

- [ ] **Step 2: Update targets pipeline**

In `targets/pioneer_targets.R`, find where the stan data targets are defined. Add the new flags as configurable options that flow from the `tar_map` dimensions or a configuration list. Default to `FALSE` (current behavior) so existing pipeline runs are unaffected.

- [ ] **Step 3: Run existing tests**

Run: `Rscript -e 'testthat::test_dir("tests/testthat")'`
Expected: All existing tests PASS (backward-compatible defaults)

- [ ] **Step 4: Commit**

```bash
git add r/pioneer/prepare_analysis_data.R targets/pioneer_targets.R
git commit -m "Wire visit-gated flags through R data prep and targets pipeline"
```

---

### Task 8: GQ forecasting compatibility

**Files:**
- Modify: `stan/psa/_psa_endpoints_generated_quantities.stan` or `stan/modules/state_space/burden_endpoints.stan`
- Review: `stan/pfs.stanfunctions` (forecast functions)

The GQ block generates future PFS/OS endpoints by sampling from the hazard. For visit-gated 0->1:
- **Observed patients** (not right-censored): already have their progression time — no forecasting needed for 0->1.
- **Censored patients**: need to forecast when they'd progress. The forecast should use the **modeled** PSA trajectory (from `states_full_grid`) at future assessment times.

- [ ] **Step 1: Review existing forecast pathway**

Read `stan/pfs.stanfunctions` functions `assessment_gated_survival_time_rng` and `visit_only_survival_time_rng`. These already implement visit-gated sampling for GQ — they accumulate hazard only at visit weeks. Verify they can be reused for the 0->1 forecast in visit-gated mode.

- [ ] **Step 2: Check if forecast pathway needs modification**

The existing `assessment_gated_survival_time_rng` uses `log_cond_surv_01` which, in the GQ block, is still computed from the modeled trajectory (transformed parameters). In visit-gated mode during fitting, we use observed PSA for the covariate. But in GQ, we need to:
1. Keep using the modeled trajectory for `log_cond_surv_01` computation (it's in transformed parameters, always computed)
2. The forecast functions sample from this — they should work as-is

The key insight: `log_cond_surv_01` in transformed parameters is always computed for all weeks (needed for GQ). The visit-gated change only affects the **likelihood** (which weeks contribute to the log-probability). The GQ block can use the full `log_cond_surv_01` for forecasting.

However, in visit-gated mode, `log_cond_surv_01` uses observed PSA for the covariate (from `ms_observed_covar_01`). For forecasting beyond the last visit, there's no observed PSA — the carry-forward value from the last observation would be used. This is acceptable for now but could be improved later by switching to the modeled trajectory for forecast-only weeks.

- [ ] **Step 3: Document the forecasting behavior**

Add a comment in the GQ section of `pioneer.stan` noting:
```stan
// In visit-gated mode, log_cond_surv_01 uses observed PSA (carry-forward)
// for the covariate. For forecasting beyond last visit, the carry-forward
// value from the patient's last observed PSA is used. This is conservative
// (assumes PSA stays at last observed level).
```

- [ ] **Step 4: Commit**

```bash
git add stan/psa/pioneer.stan
git commit -m "Document GQ forecasting behavior in visit-gated mode"
```

---

### Task 9: End-to-end validation

- [ ] **Step 1: Syntax-check all Stan models**

```bash
~/.cmdstan/cmdstan-2.38.0/bin/stanc --include-paths=stan --include-paths=stan/psa stan/psa/pioneer.stan
~/.cmdstan/cmdstan-2.38.0/bin/stanc --include-paths=stan --include-paths=stan/psa stan/psa/ms-standalone.stan
~/.cmdstan/cmdstan-2.38.0/bin/stanc --include-paths=stan --include-paths=stan/tumor stan/tumor/sf-ssm-log-space.stan
~/.cmdstan/cmdstan-2.38.0/bin/stanc --include-paths=stan --include-paths=stan/tumor stan/tumor/ms-standalone.stan
```
Expected: All clean

- [ ] **Step 2: Run all R tests**

```bash
Rscript -e 'testthat::test_dir("tests/testthat")'
```
Expected: All PASS

- [ ] **Step 3: Test backward compatibility**

Build stan data with `visit_gated_01 = FALSE` (default) and verify it produces identical stan data to current main branch. The model should compile and sample identically.

- [ ] **Step 4: Test visit-gated mode**

Build stan data with `visit_gated_01 = TRUE, enable_02_time_varying_cov = FALSE` and verify:
- `n_time_varying_covar_01 == 1`
- `n_time_varying_covar_02 == 0`
- `enable_ms_visit_gated_01 == 1`
- Stan model compiles
- Prior sampling runs without errors (short chain, few iterations)

- [ ] **Step 5: Commit final state**

```bash
git add -A
git commit -m "Visit-gated 0->1 transition: end-to-end validation complete"
```
