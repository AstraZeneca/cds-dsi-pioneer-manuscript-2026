# Per-level Centered Parameterization (`LEVEL_MODE_RE_CP`) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Extend the existing `LEVEL_MODE_*` enum with a new `RE_CP` value so every hierarchy level in `tr`, `frac`, `init`, and `multistate` modules can independently switch between non-centered (NCP) and centered (CP) parameterization, without breaking `FE` or `RE_GP` behavior.

**Architecture:** Currently every hierarchical slice (intercept or slope) declares one `_raw_` parameter vector sampled `std_normal()` / `student_t(ν,0,1)` and scaled `sd * raw` in `transformed_parameters`. Stan cannot mix parameterizations inside one declared vector. The fix: for each slice, split enabled groups into two buckets by level mode — a `_raw_` vector (sized by groups at `FE`, `RE`, `RE_GP` levels, kept exactly as today) and a new `_cp_` vector (sized by groups at `RE_CP` levels) sampled directly `~normal(0, sd)` / `student_t(ν,0,sd)`. `transformed_parameters` glues both back into the existing `_scaled_` vector; downstream gather code is untouched. When all levels stay at their existing modes (no `RE_CP` anywhere), the model is bit-exact to `main`.

**Tech Stack:** Stan (cmdstan 2.38), R (targets, testthat), cmdstanr.

**Key decisions (locked in by user):**
- Name for centered vectors: `_cp_` (parallels existing `_raw_`)
- Granularity: per-level — the existing mode enum drives both intercept and slope at each level. `enable_level_cov_*` stays as a 0/1 gate.
- Scope of PR 1: scaffolding only, all defaults bit-exact vs `main`. A separate PR flips `enable_level_intercept_init[patient] = RE_CP` for the pioneer combined variant.
- `RE_GP` does NOT get a CP variant. Only plain `RE` → `RE_CP` mapping.

---

## File Structure

**Stan — new constant + shared helper:**
- `stan/_hierarchy_transformed_data.stan` — add `LEVEL_MODE_RE_CP = 4` constant alongside the existing three.
- `stan/hierarchy.stanfunctions` — relax the range check (`0-3` → `0-4`) and add a new helper `split_cp_ncp_pos()` that returns `(n_raw, raw_pos, n_cp, cp_pos)` given a mode vector.
- `stan/multistate.stanfunctions` — relax the range check in `compute_ms_level_baseline_flags` (`0-3` → `0-4`) and update `any_re_level` / `is_gp` / `is_enabled` masks to treat `RE_CP` correctly.

**Stan — modules (scaffolding, one slice pair per module):**
- `stan/modules/{tr,frac,init}/flags.stan` — widen the upper bound of `enable_level_intercept_*` from `upper=3` to `upper=4`.
- `stan/modules/{tr,frac,init}/transformed_data.stan` — compute `n_raw_groups_*_intercept`, `n_cp_groups_*_intercept`, `raw_level_pos_*_intercept`, `cp_level_pos_*_intercept` (and the slope analogues). Keep existing `n_enabled_groups_*` for backward-compatible sizing of `_scaled_`.
- `stan/modules/{tr,frac,init}/parameters.stan` — rename raw vectors to be sized by `n_raw_groups_*` and add new `_cp_` vectors sized by `n_cp_groups_*`. Same pattern for slopes.
- `stan/modules/{tr,frac,init}/priors.stan` — add CP branch that samples `_cp_` directly.
- `stan/modules/{tr,frac,init}/transformed_parameters.stan` — glue path: if `mode == RE_CP`, `_scaled_[e_lo:e_hi] = _cp_[c_lo:c_hi]`; else `_scaled_[e_lo:e_hi] = sd[lv] * _raw_[r_lo:r_hi]`. Existing downstream gather (`patient_*_flat_idx`) is unchanged.
- `stan/modules/multistate/flags.stan` — widen `enable_ms_level_baseline_hazard` upper bound `3` → `4`. (`enable_ms_level_cov` stays `0/1`.)
- `stan/modules/multistate/transformed_data.stan` — per transition (01, 02, 12_s, 12_t, 03, 32), compute `n_raw_groups_ms_baseline_<tr>`, `n_cp_groups_ms_baseline_<tr>`, `raw_level_pos_ms_baseline_<tr>`, `cp_level_pos_ms_baseline_<tr>`. For slopes, follow the mode from `enable_ms_level_baseline_hazard` since the user-facing mode is per-level, not per-transition — slope RE_CP routing uses the same mode vector (see implementation notes below).
- `stan/modules/multistate/parameters.stan` — split each `raw_log_lambda_gp_<tr>_level_intercept` into a `_raw_` and `_cp_` vector. Split each `raw_level_slope_<tr>`.
- `stan/modules/multistate/priors.stan` — add CP branches per transition (mirror the intercept+slope pattern).
- `stan/modules/multistate/transformed_parameters.stan` — glue per transition.

**R:**
- `targets/pioneer_targets.R:70-75` — add `re_cp = 4L` to the `level_intercept_mode` lookup.
- `targets/sclc_targets.R` — add `re_cp = 4L` if this file defines a similar lookup (Task 8 checks).
- `r/pioneer/priors.R:28-34` — leave `fe_mode` logic untouched. Document that `RE_CP` levels use the same RE SD hyperparameter as `RE` (no new prior fields).
- `r/pioneer/initializers.R` — add CP routing: groups at `RE_CP` levels get initialized into `_cp_` with `rnorm(n, 0, sd_level[lv])`; other enabled groups stay in `_raw_` with current `rnorm(n, 0, 0.2)` / zeros pattern.
- `r/sclc/initializers.R` and `r/sclc/initializers_fixed.R` — same CP routing additions.

**Tests:**
- `tests/testthat/stan/test_sd_expansion.stan` — extend to include `RE_CP` in the mode enum (mode value `4`); `RE_CP` uses the same `sd_expanded[lv]` as `RE` (free SD param), so the existing expansion logic is just widened.
- `tests/testthat/test-stan-sd-expansion.R` — add test cases for `RE_CP` mode (single-level, mixed).
- `tests/testthat/stan/test_split_cp_ncp_pos.stan` — NEW: unit test for the position-splitter helper.
- `tests/testthat/test-stan-split-cp-ncp-pos.R` — NEW: R oracle + test harness for the helper.
- `tests/testthat/test-pioneer-bit-exact.R` — NEW: compile the full pioneer model with all-default (no `RE_CP`) config and diff generated `model.hpp` (or a small-data `draws()` hash) against a saved baseline.

---

## Implementation Notes (read before starting)

1. **Why `enable_level_cov_*` stays 0/1:** The user chose per-level granularity — one mode per level drives both the intercept and its slopes. When `enable_level_intercept_*[lv] == RE_CP`, the slope slice at that level also flips to CP *if* `enable_level_cov_*[lv] == 1`. The slope mode is read from `enable_level_intercept_*`, not from `enable_level_cov_*`. This keeps the flag surface small and matches the user's "per-level" answer.

2. **FE preservation:** FE groups never enter the `_cp_` bucket. The SD assembly block continues to dispatch on `mode == LEVEL_MODE_FE` (uses `*_fe_sd_level_intercept[lv]`); `mode == LEVEL_MODE_RE || mode == LEVEL_MODE_RE_CP` takes the next slot from `_sd_level_intercept_raw`. FE × CP is a non-combination by design — there's no funnel to break there.

3. **RE_GP is untouched:** `RE_GP` levels continue to route through the `_raw_` bucket (GP residual is a separate pooling mechanism). The only change `RE_GP` sees is that `LEVEL_MODE_RE_CP = 4` now exists above it in the enum.

4. **Prior statement ordering:** Prior statements in `priors.stan` reference `init_sd_level_intercept[lv]` (assembled in `transformed_parameters.stan`). This works because transformed_parameters executes before the model block each iteration — already relied on elsewhere.

5. **Bit-exact means identical posterior samples** when all levels stay at their existing `NONE/FE/RE/RE_GP` values. Any `RE_CP` in config would produce mathematically equivalent posterior means but different geometry (different MCMC trace). Task 14 locks this in with a diff test.

6. **Propensity module:** The `propensity` module has its own `beta_propensity` parameters (not level-indexed). Out of scope.

7. **Process-noise NCP blocks** (`tr_raw_patient_log_sd_process_noise`, `tr_raw_patient_phi_process_noise`, `tr_raw_patient_process_noise`, `tr_raw_pop_process_noise`) are patient-wise single-level, not a level hierarchy. Out of scope for this PR. Revisit separately if diagnostics flag them.

---

## Task 1: Add `LEVEL_MODE_RE_CP` constant

**Files:**
- Modify: `stan/_hierarchy_transformed_data.stan:8-11`

- [ ] **Step 1: Read current constants**

Run: `cat stan/_hierarchy_transformed_data.stan | sed -n '7,12p'`
Expected output:
```stan
// Level intercept mode constants (consistent across tr, frac, init, ms modules)
int LEVEL_MODE_NONE  = 0;  // No intercept at this level
int LEVEL_MODE_FE    = 1;  // Fixed effect: SD is a data hyperparameter (no pooling)
int LEVEL_MODE_RE    = 2;  // Random effect: SD estimated via NCP (hierarchical pooling)
int LEVEL_MODE_RE_GP = 3;  // Random effect + full GP residual (ms module only, for now)
```

- [ ] **Step 2: Add the new constant**

Edit `stan/_hierarchy_transformed_data.stan`:
```stan
// Level intercept mode constants (consistent across tr, frac, init, ms modules)
int LEVEL_MODE_NONE  = 0;  // No intercept at this level
int LEVEL_MODE_FE    = 1;  // Fixed effect: SD is a data hyperparameter (no pooling)
int LEVEL_MODE_RE    = 2;  // Random effect, non-centered (NCP): raw ~ std_normal, scaled = sd * raw
int LEVEL_MODE_RE_GP = 3;  // Random effect + full GP residual (ms module only, for now)
int LEVEL_MODE_RE_CP = 4;  // Random effect, centered (CP): centered ~ normal(0, sd) directly
```

- [ ] **Step 3: Commit**

```bash
git add stan/_hierarchy_transformed_data.stan
git commit -m "stan: add LEVEL_MODE_RE_CP=4 constant for centered parameterization"
```

---

## Task 2: Relax range checks in hierarchy.stanfunctions

**Files:**
- Modify: `stan/hierarchy.stanfunctions:47-52`

- [ ] **Step 1: Read current validator**

Run: `sed -n '47,52p' stan/hierarchy.stanfunctions`

- [ ] **Step 2: Widen the upper bound**

Edit `stan/hierarchy.stanfunctions` — change both validator lines:
```stan
for (lv in 1:n_levels) {
  if (enable_intercept[lv] < 0 || enable_intercept[lv] > 4)
    fatal_error("enable_intercept[", lv, "] must be 0-4; got ", enable_intercept[lv]);
  if (enable_slope[lv] < 0 || enable_slope[lv] > 4)
    fatal_error("enable_slope[", lv, "] must be 0-4; got ", enable_slope[lv]);
}
```

- [ ] **Step 3: Commit**

```bash
git add stan/hierarchy.stanfunctions
git commit -m "stan: widen hierarchy validator range to 0-4 for RE_CP"
```

---

## Task 3: Add `split_cp_ncp_pos` helper

Every module needs to know, per slice, (a) how many groups go into `_raw_`, (b) how many into `_cp_`, and (c) the position arrays into both. Rather than rewriting this in 11 slices, add a single helper.

**Files:**
- Modify: `stan/hierarchy.stanfunctions` (append new function)
- Create: `tests/testthat/stan/test_split_cp_ncp_pos.stan`
- Create: `tests/testthat/test-stan-split-cp-ncp-pos.R`

- [ ] **Step 1: Write the failing Stan test harness**

Create `tests/testthat/stan/test_split_cp_ncp_pos.stan`:
```stan
functions {
  #include pos.stanfunctions
  #include hierarchy.stanfunctions
}

data {
  int<lower=1> n_levels;
  array[n_levels] int<lower=0,upper=4> mode;
  array[n_levels] int<lower=0> n_groups;
}

generated quantities {
  int n_raw;
  int n_cp;
  array[n_levels + 1] int raw_pos;
  array[n_levels + 1] int cp_pos;
  (n_raw, raw_pos, n_cp, cp_pos) = split_cp_ncp_pos(n_levels, n_groups, mode);
}
```

- [ ] **Step 2: Write the failing R test**

Create `tests/testthat/test-stan-split-cp-ncp-pos.R`:
```r
library(testthat)
library(here)
library(posterior)

# R oracle: route a group to _raw_ if mode in {1 FE, 2 RE, 3 RE_GP}; route to
# _cp_ if mode == 4 RE_CP. Skip levels with mode == 0 NONE.
r_split_cp_ncp_pos <- function(mode, n_groups) {
  n_levels <- length(mode)
  raw_counts <- integer(n_levels)
  cp_counts  <- integer(n_levels)
  for (lv in seq_len(n_levels)) {
    if (mode[lv] == 4L) {
      cp_counts[lv] <- n_groups[lv]
    } else if (mode[lv] %in% c(1L, 2L, 3L)) {
      raw_counts[lv] <- n_groups[lv]
    }
  }
  list(
    n_raw   = sum(raw_counts),
    n_cp    = sum(cp_counts),
    raw_pos = as.integer(c(1L, cumsum(raw_counts) + 1L)),
    cp_pos  = as.integer(c(1L, cumsum(cp_counts)  + 1L))
  )
}

run_split_test <- function(mode, n_groups) {
  n_levels <- length(mode)
  stan_data <- list(
    n_levels = n_levels,
    mode     = as.array(as.integer(mode)),
    n_groups = as.array(as.integer(n_groups))
  )
  fit <- test_stan_function(
    here("tests", "testthat", "stan", "test_split_cp_ncp_pos.stan"),
    data = stan_data
  )
  d <- as_draws_df(fit$draws())
  exp <- r_split_cp_ncp_pos(mode, n_groups)
  expect_equal(as.integer(d$n_raw[1]), exp$n_raw)
  expect_equal(as.integer(d$n_cp[1]),  exp$n_cp)
  for (i in seq_len(n_levels + 1)) {
    expect_equal(as.integer(d[[paste0("raw_pos[", i, "]")]][1]), exp$raw_pos[i])
    expect_equal(as.integer(d[[paste0("cp_pos[",  i, "]")]][1]), exp$cp_pos[i])
  }
}

test_that("split: all NONE", {
  run_split_test(c(0L, 0L, 0L), c(5L, 10L, 100L))
})

test_that("split: all RE (legacy)", {
  run_split_test(c(2L, 2L, 2L), c(5L, 10L, 100L))
})

test_that("split: all RE_CP", {
  run_split_test(c(4L, 4L, 4L), c(5L, 10L, 100L))
})

test_that("split: mixed NONE/FE/RE/RE_CP 4-level", {
  run_split_test(c(0L, 1L, 2L, 4L), c(3L, 5L, 10L, 100L))
})

test_that("split: FE goes to raw bucket", {
  run_split_test(c(1L, 4L), c(5L, 50L))
})

test_that("split: RE_GP goes to raw bucket", {
  run_split_test(c(3L, 4L), c(5L, 50L))
})
```

- [ ] **Step 3: Run test to verify it fails**

Run: `Rscript -e 'testthat::test_file("tests/testthat/test-stan-split-cp-ncp-pos.R")'`
Expected: FAIL with message about `split_cp_ncp_pos` not found.

- [ ] **Step 4: Implement the helper**

Append to `stan/hierarchy.stanfunctions`:
```stan
/**
 * Split level groups into NCP (raw) and CP (cp) buckets by mode.
 *
 * Groups at levels with mode in {FE=1, RE=2, RE_GP=3} go into the `_raw_` bucket
 * (sampled std_normal or student_t(nu,0,1), scaled by sd). Groups at levels with
 * mode == RE_CP (4) go into the `_cp_` bucket (sampled normal(0, sd) directly).
 * Groups at levels with mode == NONE (0) are skipped.
 *
 * @param n_levels  Number of hierarchy levels
 * @param n_groups  [n_levels] group counts per level (use n_forecast_groups_per_level)
 * @param mode      [n_levels] level mode: 0=NONE, 1=FE, 2=RE, 3=RE_GP, 4=RE_CP
 * @return Tuple:
 *   (n_raw, raw_pos[n_levels+1], n_cp, cp_pos[n_levels+1])
 */
tuple(int, array[] int, int, array[] int) split_cp_ncp_pos(
  int n_levels,
  array[] int n_groups,
  array[] int mode
) {
  for (lv in 1:n_levels) {
    if (mode[lv] < 0 || mode[lv] > 4)
      fatal_error("mode[", lv, "] must be 0-4; got ", mode[lv]);
  }
  array[n_levels] int raw_mask;
  array[n_levels] int cp_mask;
  for (lv in 1:n_levels) {
    raw_mask[lv] = (mode[lv] == 1 || mode[lv] == 2 || mode[lv] == 3) ? 1 : 0;
    cp_mask[lv]  = (mode[lv] == 4) ? 1 : 0;
  }
  int n_raw = compute_n_enabled_groups(n_groups, raw_mask);
  int n_cp  = compute_n_enabled_groups(n_groups, cp_mask);
  array[n_levels + 1] int raw_pos = create_enabled_pos(n_groups, raw_mask);
  array[n_levels + 1] int cp_pos  = create_enabled_pos(n_groups, cp_mask);
  return (n_raw, raw_pos, n_cp, cp_pos);
}
```

- [ ] **Step 5: Run test to verify it passes**

Run: `Rscript -e 'testthat::test_file("tests/testthat/test-stan-split-cp-ncp-pos.R")'`
Expected: all tests PASS.

- [ ] **Step 6: Commit**

```bash
git add stan/hierarchy.stanfunctions tests/testthat/stan/test_split_cp_ncp_pos.stan tests/testthat/test-stan-split-cp-ncp-pos.R
git commit -m "stan: add split_cp_ncp_pos helper for RE_CP bucket routing"
```

---

## Task 4: Extend SD expansion test harness to RE_CP

`RE_CP` uses the same free-SD parameter as `RE` (both take the next slot from `_sd_level_intercept_raw`). The SD assembly logic just needs to treat mode `4` identically to mode `2` for SD lookup.

**Files:**
- Modify: `tests/testthat/stan/test_sd_expansion.stan:9,23-33`
- Modify: `tests/testthat/test-stan-sd-expansion.R`

- [ ] **Step 1: Widen the Stan test data upper bound**

Edit `tests/testthat/stan/test_sd_expansion.stan`:
```stan
data {
  int<lower=1> n_levels;
  // Mode per level: 0=none, 1=fe, 2=re, 4=re_cp (3=re_gp not tested here; same sd path as re)
  array[n_levels] int<lower=0,upper=4> mode;
  array[n_levels] real<lower=0> fe_sd;
  int<lower=0> n_re_levels;  // counts mode == 2 OR mode == 4
  array[n_re_levels] real<lower=0> re_sd_values;
}

generated quantities {
  array[n_levels] real<lower=0> sd_expanded;
  {
    int sd_idx = 0;
    for (lv in 1:n_levels) {
      if (mode[lv] == 1) {           // LEVEL_MODE_FE
        sd_expanded[lv] = fe_sd[lv];
      } else if (mode[lv] == 2 || mode[lv] == 4) {  // LEVEL_MODE_RE or LEVEL_MODE_RE_CP
        sd_idx += 1;
        sd_expanded[lv] = re_sd_values[sd_idx];
      } else {
        sd_expanded[lv] = 0.0;
      }
    }
  }
}
```

- [ ] **Step 2: Update the R oracle**

Edit `tests/testthat/test-stan-sd-expansion.R` — change the `r_expand_sd` function:
```r
r_expand_sd <- function(mode, fe_sd, re_sd_values) {
  n_levels <- length(mode)
  sd_out <- numeric(n_levels)
  sd_idx <- 0L
  for (lv in seq_len(n_levels)) {
    if (mode[lv] == 1L) {
      sd_out[lv] <- fe_sd[lv]
    } else if (mode[lv] == 2L || mode[lv] == 4L) {
      sd_idx <- sd_idx + 1L
      sd_out[lv] <- re_sd_values[sd_idx]
    } else {
      sd_out[lv] <- 0.0
    }
  }
  sd_out
}

run_sd_expansion_test <- function(mode, fe_sd, re_sd_values) {
  n_levels <- length(mode)
  n_re_levels <- sum(mode == 2L | mode == 4L)
  stopifnot(length(re_sd_values) == n_re_levels)
  # ... rest unchanged
```

- [ ] **Step 3: Add new test cases**

Append to `tests/testthat/test-stan-sd-expansion.R`:
```r
test_that("SD expansion: all RE_CP (mode=4)", {
  run_sd_expansion_test(
    mode         = c(4L, 4L, 4L),
    fe_sd        = c(0.1, 0.2, 0.3),
    re_sd_values = c(0.4, 0.5, 0.6)
  )
})

test_that("SD expansion: mixed NONE/FE/RE/RE_CP 4-level", {
  run_sd_expansion_test(
    mode         = c(0L, 1L, 2L, 4L),
    fe_sd        = c(0.1, 0.05, 0.2, 0.9),
    re_sd_values = c(0.35, 0.80)
  )
})

test_that("SD expansion: pioneer level-3 CP use case", {
  # trial=none, arm=re, patient=re_cp (the funnel fix)
  run_sd_expansion_test(
    mode         = c(0L, 2L, 4L),
    fe_sd        = c(0.1, 0.35, 0.5),
    re_sd_values = c(0.35, 0.42)
  )
})
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `Rscript -e 'testthat::test_file("tests/testthat/test-stan-sd-expansion.R")'`
Expected: all tests PASS (including the new RE_CP cases).

- [ ] **Step 5: Commit**

```bash
git add tests/testthat/stan/test_sd_expansion.stan tests/testthat/test-stan-sd-expansion.R
git commit -m "test: extend SD expansion tests to cover LEVEL_MODE_RE_CP"
```

---

## Task 5: Extend `init` module — intercept slice

Split `init_raw_level_intercept` into `_raw_` + `_cp_` buckets, add the CP prior branch, and wire the glue in `transformed_parameters`.

**Files:**
- Modify: `stan/modules/init/flags.stan:8`
- Modify: `stan/modules/init/transformed_data.stan`
- Modify: `stan/modules/init/parameters.stan:20`
- Modify: `stan/modules/init/priors.stan:13-47`
- Modify: `stan/modules/init/transformed_parameters.stan:28-46`

- [ ] **Step 1: Widen the flag upper bound**

Edit `stan/modules/init/flags.stan`:
```stan
array[n_levels] int<lower=0,upper=4> enable_level_intercept_init;
```

- [ ] **Step 2: Add bucket counts + positions to transformed_data**

Edit `stan/modules/init/transformed_data.stan` — append after the existing `compute_level_module_flags` call:
```stan
// Split intercept groups into NCP (_raw_) and CP (_cp_) buckets by level mode.
// Groups at NONE levels are skipped; FE/RE/RE_GP go to raw, RE_CP goes to cp.
int n_raw_groups_init_intercept;
int n_cp_groups_init_intercept;
array[n_levels + 1] int raw_level_pos_init_intercept;
array[n_levels + 1] int cp_level_pos_init_intercept;
(n_raw_groups_init_intercept, raw_level_pos_init_intercept,
 n_cp_groups_init_intercept,  cp_level_pos_init_intercept) =
  split_cp_ncp_pos(n_levels, n_forecast_groups_per_level, enable_level_intercept_init);

// Count RE-only levels for right-sizing the SD parameter array — RE_CP also needs
// a free SD hyperparameter (same as RE), so both count here.
int n_re_levels_init_intercept = 0;
for (lv in 1:n_levels)
  if (enable_level_intercept_init[lv] == LEVEL_MODE_RE ||
      enable_level_intercept_init[lv] == LEVEL_MODE_RE_CP) n_re_levels_init_intercept += 1;
```

Replace the existing `n_re_levels_init_intercept` loop with the version above.

- [ ] **Step 3: Split the raw parameter, add the cp parameter**

Edit `stan/modules/init/parameters.stan`:
```stan
// Raw standard normal draws for intercepts at FE/RE/RE_GP levels
vector[n_raw_groups_init_intercept] init_raw_level_intercept;

// Centered draws for intercepts at RE_CP levels — sampled ~normal(0, sd) directly
vector[n_cp_groups_init_intercept] init_cp_level_intercept;
```

(Replaces the single `vector[n_enabled_groups_init_intercept] init_raw_level_intercept;` line.)

- [ ] **Step 4: Update priors — SD prior now fires for RE OR RE_CP; add CP branch**

Edit `stan/modules/init/priors.stan` inside the `for (lv in 1:n_levels)` block:
```stan
for (lv in 1:n_levels) {
  int mode = enable_level_intercept_init[lv];

  // SD prior for any level that samples its SD as a free parameter (RE or RE_CP)
  if (mode == LEVEL_MODE_RE || mode == LEVEL_MODE_RE_CP) {
    sd_idx += 1;
    init_sd_level_intercept_raw[sd_idx] ~ normal(0, init_sd_level_intercept_sd[lv]);
  }
  if (n_covar > 0) {
    init_sd_level_slope[lv] ~ normal(0, init_sd_level_slope_sd[lv]);
  }

  // RAW path: FE, RE, RE_GP all sample from std_normal / student_t(nu, 0, 1)
  if (mode == LEVEL_MODE_FE || mode == LEVEL_MODE_RE || mode == LEVEL_MODE_RE_GP) {
    int r_lo = raw_level_pos_init_intercept[lv];
    int r_hi = raw_level_pos_init_intercept[lv + 1] - 1;
    if (r_hi >= r_lo) {
      if (enable_student_t_hierarchy)
        init_raw_level_intercept[r_lo:r_hi] ~ student_t(init_nu_level[lv], 0, 1);
      else
        init_raw_level_intercept[r_lo:r_hi] ~ std_normal();
    }
  }

  // CP path: RE_CP samples directly at scale sd
  if (mode == LEVEL_MODE_RE_CP) {
    int c_lo = cp_level_pos_init_intercept[lv];
    int c_hi = cp_level_pos_init_intercept[lv + 1] - 1;
    if (c_hi >= c_lo) {
      if (enable_student_t_hierarchy)
        init_cp_level_intercept[c_lo:c_hi]
          ~ student_t(init_nu_level[lv], 0, init_sd_level_intercept[lv]);
      else
        init_cp_level_intercept[c_lo:c_hi]
          ~ normal(0, init_sd_level_intercept[lv]);
    }
  }

  // Slope raw effects — unchanged in this task (slopes handled in Task 6)
  if (enable_level_cov_init[lv] && n_covar > 0) {
    int s_lo = enabled_level_pos_init_slope[lv];
    int s_hi = enabled_level_pos_init_slope[lv + 1] - 1;
    if (enable_student_t_hierarchy)
      to_vector(init_raw_level_slope[s_lo:s_hi, :]) ~ student_t(init_nu_level[lv], 0, 1);
    else
      to_vector(init_raw_level_slope[s_lo:s_hi, :]) ~ std_normal();
  }
}
```

(Replaces the existing body of the `for (lv in 1:n_levels)` block.)

- [ ] **Step 5: Wire the glue in transformed_parameters**

Edit `stan/modules/init/transformed_parameters.stan` — replace the existing INTERCEPT EFFECTS Step 1 block:
```stan
// ===== INTERCEPT EFFECTS =====
// Step 1: Assemble scaled intercepts from either raw (NCP) or cp (CP) bucket.
vector[n_enabled_groups_init_intercept] init_scaled_level_intercept;
for (lv in 1:n_levels) {
  int mode = enable_level_intercept_init[lv];
  if (mode == LEVEL_MODE_NONE) continue;
  int e_lo, e_hi;
  (e_lo, e_hi) = get_pos(enabled_level_pos_init_intercept, lv);

  if (mode == LEVEL_MODE_RE_CP) {
    // CP: scaled = centered (identity — the _cp_ vector is already at natural scale)
    int c_lo = cp_level_pos_init_intercept[lv];
    int c_hi = cp_level_pos_init_intercept[lv + 1] - 1;
    init_scaled_level_intercept[e_lo:e_hi] = init_cp_level_intercept[c_lo:c_hi];
  } else {
    // NCP (FE, RE, RE_GP): scaled = sd[lv] * raw
    int r_lo = raw_level_pos_init_intercept[lv];
    int r_hi = raw_level_pos_init_intercept[lv + 1] - 1;
    init_scaled_level_intercept[e_lo:e_hi] =
      init_sd_level_intercept[lv] * init_raw_level_intercept[r_lo:r_hi];
  }
}

// Step 2 (gather) — unchanged
vector[n_forecast_patients] init_linpred_level_intercepts = zeros_vector(n_forecast_patients);
for (lv in 1:n_levels) {
  if (enable_level_intercept_init[lv]) {
    init_linpred_level_intercepts += init_scaled_level_intercept[patient_init_intercept_flat_idx[forecast_patient_idx, lv]];
  }
}
```

- [ ] **Step 6: Stan syntax check**

Run:
```bash
~/.cmdstan/cmdstan-2.38.0/bin/stanc \
  --include-paths=stan --include-paths=stan/psa \
  stan/psa/pioneer.stan
```
Expected: exits with code 0 (no parse errors).

- [ ] **Step 7: Commit**

```bash
git add stan/modules/init/flags.stan stan/modules/init/transformed_data.stan \
        stan/modules/init/parameters.stan stan/modules/init/priors.stan \
        stan/modules/init/transformed_parameters.stan
git commit -m "stan(init): split intercept into _raw_+_cp_ buckets for RE_CP"
```

---

## Task 6: Extend `init` module — slope slice

The slope slice is controlled by `enable_level_cov_init` (a 0/1 gate), but the *mode* (NCP vs CP) is read from `enable_level_intercept_init` per the per-level decision. Slope follows its level's intercept mode.

**Files:**
- Modify: `stan/modules/init/transformed_data.stan` (add slope bucket counts)
- Modify: `stan/modules/init/parameters.stan:26` (split matrix)
- Modify: `stan/modules/init/priors.stan` (CP branch for slopes)
- Modify: `stan/modules/init/transformed_parameters.stan:49-62` (glue path for slopes)

- [ ] **Step 1: Add slope bucket counts to transformed_data**

Append to `stan/modules/init/transformed_data.stan`:
```stan
// Split slope groups the same way, but routed by a combined mask: slope is only
// enabled when enable_level_cov_init[lv] == 1 AND the level has an intercept.
// When routed, the slope parameterization follows the level's intercept mode.
array[n_levels] int init_slope_mode;
for (lv in 1:n_levels) {
  init_slope_mode[lv] = enable_level_cov_init[lv] ? enable_level_intercept_init[lv] : 0;
}
int n_raw_groups_init_slope;
int n_cp_groups_init_slope;
array[n_levels + 1] int raw_level_pos_init_slope;
array[n_levels + 1] int cp_level_pos_init_slope;
(n_raw_groups_init_slope, raw_level_pos_init_slope,
 n_cp_groups_init_slope,  cp_level_pos_init_slope) =
  split_cp_ncp_pos(n_levels, n_forecast_groups_per_level, init_slope_mode);
```

- [ ] **Step 2: Split the slope matrix**

Edit `stan/modules/init/parameters.stan`:
```stan
// Raw slope effects at FE/RE/RE_GP levels
matrix[n_raw_groups_init_slope, n_covar] init_raw_level_slope;

// Centered slope effects at RE_CP levels
matrix[n_cp_groups_init_slope, n_covar] init_cp_level_slope;
```

- [ ] **Step 3: Update the slope prior branch**

In `stan/modules/init/priors.stan`, replace the slope block inside the level loop:
```stan
// Slope effects — parameterization follows the level's intercept mode
if (enable_level_cov_init[lv] && n_covar > 0) {
  int mode = enable_level_intercept_init[lv];

  // RAW path (FE, RE, RE_GP)
  if (mode == LEVEL_MODE_FE || mode == LEVEL_MODE_RE || mode == LEVEL_MODE_RE_GP) {
    int r_lo = raw_level_pos_init_slope[lv];
    int r_hi = raw_level_pos_init_slope[lv + 1] - 1;
    if (r_hi >= r_lo) {
      if (enable_student_t_hierarchy)
        to_vector(init_raw_level_slope[r_lo:r_hi, :]) ~ student_t(init_nu_level[lv], 0, 1);
      else
        to_vector(init_raw_level_slope[r_lo:r_hi, :]) ~ std_normal();
    }
  }

  // CP path (RE_CP): sample at scale sd_level_slope[lv]
  if (mode == LEVEL_MODE_RE_CP) {
    int c_lo = cp_level_pos_init_slope[lv];
    int c_hi = cp_level_pos_init_slope[lv + 1] - 1;
    if (c_hi >= c_lo) {
      for (k in 1:n_covar) {
        if (enable_student_t_hierarchy)
          init_cp_level_slope[c_lo:c_hi, k]
            ~ student_t(init_nu_level[lv], 0, init_sd_level_slope[lv, k]);
        else
          init_cp_level_slope[c_lo:c_hi, k]
            ~ normal(0, init_sd_level_slope[lv, k]);
      }
    }
  }
}
```

- [ ] **Step 4: Wire the slope glue in transformed_parameters**

Edit `stan/modules/init/transformed_parameters.stan` — replace the SLOPE Step 1 block:
```stan
// ===== COVARIATE SLOPE EFFECTS =====
matrix[n_enabled_groups_init_slope, n_covar] init_scaled_level_slope;
if (n_covar > 0 && n_enabled_groups_init_slope > 0) {
  for (lv in 1:n_levels) {
    if (!enable_level_cov_init[lv]) continue;
    int mode = enable_level_intercept_init[lv];
    int e_lo, e_hi;
    (e_lo, e_hi) = get_pos(enabled_level_pos_init_slope, lv);

    if (mode == LEVEL_MODE_RE_CP) {
      int c_lo = cp_level_pos_init_slope[lv];
      int c_hi = cp_level_pos_init_slope[lv + 1] - 1;
      init_scaled_level_slope[e_lo:e_hi, :] = init_cp_level_slope[c_lo:c_hi, :];
    } else {
      int r_lo = raw_level_pos_init_slope[lv];
      int r_hi = raw_level_pos_init_slope[lv + 1] - 1;
      init_scaled_level_slope[e_lo:e_hi, :] =
        init_raw_level_slope[r_lo:r_hi, :] .*
        rep_matrix(init_sd_level_slope[lv]', r_hi - r_lo + 1);
    }
  }
}
// Step 2 (gather) — unchanged
```

- [ ] **Step 5: Stan syntax check**

Run:
```bash
~/.cmdstan/cmdstan-2.38.0/bin/stanc \
  --include-paths=stan --include-paths=stan/psa \
  stan/psa/pioneer.stan
```
Expected: exits with code 0.

- [ ] **Step 6: Commit**

```bash
git add stan/modules/init/
git commit -m "stan(init): split slope into _raw_+_cp_ buckets for RE_CP"
```

---

## Task 7: Extend `tr` module — intercept + slope slices

Mechanical mirror of Tasks 5 + 6, applied to `tr`. Same 6 file edits, same pattern. The only differences are prefix (`tr_` not `init_`), and `tr`'s extra process-noise NCP blocks which are out of scope.

**Files:**
- Modify: `stan/modules/tr/{flags,transformed_data,parameters,priors,transformed_parameters}.stan`

- [ ] **Step 1: Widen flag bound**

Edit `stan/modules/tr/flags.stan` — change `upper=3` to `upper=4` on `enable_level_intercept_tr`.

- [ ] **Step 2: Add bucket counts to transformed_data**

Edit `stan/modules/tr/transformed_data.stan` — append the intercept bucket + slope bucket blocks (copy from Tasks 5.2 and 6.1, replacing `init_` with `tr_`, and use `enable_level_cov_tr` for the slope mask).

Update the `n_re_levels_tr_intercept` loop to include both `LEVEL_MODE_RE` and `LEVEL_MODE_RE_CP`.

- [ ] **Step 3: Split parameter vectors**

Edit `stan/modules/tr/parameters.stan` — replace:
```stan
vector[n_enabled_groups_tr_intercept] tr_raw_level_intercept;
```
with:
```stan
vector[n_raw_groups_tr_intercept] tr_raw_level_intercept;
vector[n_cp_groups_tr_intercept]  tr_cp_level_intercept;
```

And replace:
```stan
matrix[n_enabled_groups_tr_slope, n_covar] tr_raw_level_slope;
```
with:
```stan
matrix[n_raw_groups_tr_slope, n_covar] tr_raw_level_slope;
matrix[n_cp_groups_tr_slope,  n_covar] tr_cp_level_slope;
```

- [ ] **Step 4: Update priors block**

Edit `stan/modules/tr/priors.stan` — apply the same SD-prior relaxation, raw-path condition, CP-branch-added, slope-split pattern as in Tasks 5.4 and 6.3.

- [ ] **Step 5: Wire the glue in transformed_parameters**

Edit `stan/modules/tr/transformed_parameters.stan` — apply the same intercept glue (Task 5.5) and slope glue (Task 6.4) patterns.

- [ ] **Step 6: Stan syntax check**

Run:
```bash
~/.cmdstan/cmdstan-2.38.0/bin/stanc \
  --include-paths=stan --include-paths=stan/tumor \
  stan/tumor/sf-ssm-log-space.stan
```
Expected: exits with code 0.

- [ ] **Step 7: Commit**

```bash
git add stan/modules/tr/
git commit -m "stan(tr): split intercept and slope into _raw_+_cp_ buckets for RE_CP"
```

---

## Task 8: Extend `frac` module — intercept + slope slices

Same mechanical mirror, applied to `frac`.

**Files:**
- Modify: `stan/modules/frac/{flags,transformed_data,parameters,priors,transformed_parameters}.stan`

- [ ] **Step 1: Widen flag bound**

Edit `stan/modules/frac/flags.stan`:
```stan
array[n_levels] int<lower=0,upper=4> enable_level_intercept_frac;
```

- [ ] **Step 2: Add bucket counts to transformed_data**

Edit `stan/modules/frac/transformed_data.stan` — append intercept + slope bucket blocks (copy from Task 7.2, replacing `tr_` with `frac_`, and use `enable_level_cov_frac`).

- [ ] **Step 3: Split parameter vectors**

Edit `stan/modules/frac/parameters.stan`:
```stan
vector[n_raw_groups_frac_intercept] frac_raw_level_intercept;
vector[n_cp_groups_frac_intercept]  frac_cp_level_intercept;
matrix[n_raw_groups_frac_slope, n_covar] frac_raw_level_slope;
matrix[n_cp_groups_frac_slope,  n_covar] frac_cp_level_slope;
```

- [ ] **Step 4: Update priors block**

Apply the SD-prior + raw-path + CP-branch + slope-split pattern from Task 7.4.

- [ ] **Step 5: Wire the glue in transformed_parameters**

Apply the intercept and slope glue pattern from Task 7.5.

- [ ] **Step 6: Stan syntax check**

Run:
```bash
~/.cmdstan/cmdstan-2.38.0/bin/stanc \
  --include-paths=stan --include-paths=stan/psa \
  stan/psa/pioneer.stan
```
Expected: exits with code 0.

- [ ] **Step 7: Commit**

```bash
git add stan/modules/frac/
git commit -m "stan(frac): split intercept and slope into _raw_+_cp_ buckets for RE_CP"
```

---

## Task 9: Extend `multistate` module — baseline intercept per transition

The multistate module has 6 baseline slices (01, 02, 12_s, 12_t, 03, 32), each with its own `raw_log_lambda_gp_<tr>_level_intercept` vector. All are controlled by the single `enable_ms_level_baseline_hazard` mode vector. We split each raw vector into `_raw_` + `_cp_`.

**Files:**
- Modify: `stan/multistate.stanfunctions:73-74` (widen range check)
- Modify: `stan/modules/multistate/flags.stan:22` (widen upper bound)
- Modify: `stan/modules/multistate/transformed_data.stan` (per-transition bucket counts)
- Modify: `stan/modules/multistate/parameters.stan` (split baseline vectors)
- Modify: `stan/modules/multistate/priors.stan` (CP branches for each transition)
- Modify: `stan/modules/multistate/transformed_parameters.stan` (glue per transition)

- [ ] **Step 1: Widen ms validator**

Edit `stan/multistate.stanfunctions:73-74`:
```stan
if (enable_ms_level_baseline_hazard[lv] < 0 || enable_ms_level_baseline_hazard[lv] > 4) {
  if (strict) fatal_error("enable_ms_level_baseline_hazard[", lv, "] must be 0-4; got ", enable_ms_level_baseline_hazard[lv]);
  return (0, 0, zeros_int_array(n_levels), 0, zeros_int_array(n_levels + 1), zeros_int_array(n_levels + 1));
}
```

Also update the `any_re` computation to recognize RE_CP as a "has free SD" case:
```stan
int any_re = 0;
for (lv in 1:n_levels) {
  if (enable_ms_level_baseline_hazard[lv] == 2 || enable_ms_level_baseline_hazard[lv] == 4) {
    any_re = 1;
  }
}
```

(`is_gp` is unchanged — `RE_GP == 3` is still the only GP mode. `is_enabled` is unchanged — any mode ≥ 1 is enabled.)

- [ ] **Step 2: Widen flag bound in multistate/flags.stan**

Edit `stan/modules/multistate/flags.stan:22`:
```stan
array[n_levels] int<lower=0, upper=4> enable_ms_level_baseline_hazard;
```

Update the comment:
```stan
// 5-state flag: 0 = no level effect, 1 = FE intercept, 2 = RE intercept (NCP),
// 3 = RE + GP residual, 4 = RE intercept (CP)
```

- [ ] **Step 3: Add per-transition bucket counts to transformed_data**

Edit `stan/modules/multistate/transformed_data.stan` — after the existing B3 block (around line 83), add:
```stan
// --- Per-transition RAW vs CP bucket counts ---
// Each transition uses the shared enable_ms_level_baseline_hazard mode vector,
// but is gated on whether the transition is enabled. When the transition is
// disabled, all counts are 0.
int n_raw_groups_ms_baseline_shared;
int n_cp_groups_ms_baseline_shared;
array[n_levels + 1] int raw_level_pos_ms_baseline_shared;
array[n_levels + 1] int cp_level_pos_ms_baseline_shared;
(n_raw_groups_ms_baseline_shared, raw_level_pos_ms_baseline_shared,
 n_cp_groups_ms_baseline_shared,  cp_level_pos_ms_baseline_shared) =
  split_cp_ncp_pos(n_levels, n_forecast_groups_per_level, enable_ms_level_baseline_hazard);

int n_raw_groups_ms_baseline_01 = enable_ms_01 ? n_raw_groups_ms_baseline_shared : 0;
int n_cp_groups_ms_baseline_01  = enable_ms_01 ? n_cp_groups_ms_baseline_shared  : 0;
int n_raw_groups_ms_baseline_02 = enable_ms_02 ? n_raw_groups_ms_baseline_shared : 0;
int n_cp_groups_ms_baseline_02  = enable_ms_02 ? n_cp_groups_ms_baseline_shared  : 0;
int n_raw_groups_ms_baseline_12_s = need_12_s_gp ? n_raw_groups_ms_baseline_shared : 0;
int n_cp_groups_ms_baseline_12_s  = need_12_s_gp ? n_cp_groups_ms_baseline_shared  : 0;
int n_raw_groups_ms_baseline_12_t = need_12_t_gp ? n_raw_groups_ms_baseline_shared : 0;
int n_cp_groups_ms_baseline_12_t  = need_12_t_gp ? n_cp_groups_ms_baseline_shared  : 0;
int n_raw_groups_ms_baseline_03 = enable_ms_03 ? n_raw_groups_ms_baseline_shared : 0;
int n_cp_groups_ms_baseline_03  = enable_ms_03 ? n_cp_groups_ms_baseline_shared  : 0;
int n_raw_groups_ms_baseline_32 = enable_ms_32 ? n_raw_groups_ms_baseline_shared : 0;
int n_cp_groups_ms_baseline_32  = enable_ms_32 ? n_cp_groups_ms_baseline_shared  : 0;
// Sanity: each transition's raw+cp sum must equal the pre-existing enabled count
// (diagnostic — remove after validation if desired)
if (n_raw_groups_ms_baseline_01 + n_cp_groups_ms_baseline_01 != n_enabled_groups_ms_baseline_01)
  fatal_error("RE_CP bucket split mismatch for ms_baseline_01");
// ... repeat for each transition (optional; remove in final version if noisy)
```

- [ ] **Step 4: Split baseline parameter vectors**

Edit `stan/modules/multistate/parameters.stan` — for each transition, replace:
```stan
vector[n_enabled_groups_ms_baseline_<tr>] raw_log_lambda_gp_<tr>_level_intercept;
```
with:
```stan
vector[n_raw_groups_ms_baseline_<tr>] raw_log_lambda_gp_<tr>_level_intercept;
vector[n_cp_groups_ms_baseline_<tr>]  cp_log_lambda_gp_<tr>_level_intercept;
```

Apply this to all 6 transitions: `01`, `02`, `12_s`, `12_t`, `03`, `32`.

- [ ] **Step 5: Update priors — add CP branches per transition**

Edit `stan/modules/multistate/priors.stan` — for each transition block, replace:
```stan
if (n_enabled_groups_ms_baseline_<tr> > 0) {
  if (enable_student_t_hierarchy) {
    for (lv in 1:n_levels) {
      if (enable_ms_level_baseline_hazard[lv]) {
        int lv_start = enabled_level_pos_ms_baseline[lv];
        int lv_end = enabled_level_pos_ms_baseline[lv + 1] - 1;
        raw_log_lambda_gp_<tr>_level_intercept[lv_start:lv_end]
            ~ student_t(ms_nu_baseline_level[lv], 0, 1);
      }
    }
  } else {
    raw_log_lambda_gp_<tr>_level_intercept ~ std_normal();
  }
}
```
with:
```stan
for (lv in 1:n_levels) {
  int mode = enable_ms_level_baseline_hazard[lv];
  if (!mode) continue;

  // RAW path (FE, RE, RE_GP)
  if (mode == LEVEL_MODE_FE || mode == LEVEL_MODE_RE || mode == LEVEL_MODE_RE_GP) {
    int r_lo = raw_level_pos_ms_baseline_shared[lv];
    int r_hi = raw_level_pos_ms_baseline_shared[lv + 1] - 1;
    if (r_hi >= r_lo) {
      if (enable_student_t_hierarchy)
        raw_log_lambda_gp_<tr>_level_intercept[r_lo:r_hi]
          ~ student_t(ms_nu_baseline_level[lv], 0, 1);
      else
        raw_log_lambda_gp_<tr>_level_intercept[r_lo:r_hi] ~ std_normal();
    }
  }

  // CP path (RE_CP)
  if (mode == LEVEL_MODE_RE_CP) {
    int c_lo = cp_level_pos_ms_baseline_shared[lv];
    int c_hi = cp_level_pos_ms_baseline_shared[lv + 1] - 1;
    if (c_hi >= c_lo) {
      if (enable_student_t_hierarchy)
        cp_log_lambda_gp_<tr>_level_intercept[c_lo:c_hi]
          ~ student_t(ms_nu_baseline_level[lv], 0, log_lambda_gp_<tr>_level_intercept_sd[lv]);
      else
        cp_log_lambda_gp_<tr>_level_intercept[c_lo:c_hi]
          ~ normal(0, log_lambda_gp_<tr>_level_intercept_sd[lv]);
    }
  }
}
```

Apply to all 6 transitions. The positions `raw_level_pos_ms_baseline_shared` / `cp_level_pos_ms_baseline_shared` are shared across transitions because the enable mode vector is shared.

- [ ] **Step 6: Wire the glue in transformed_parameters**

Edit `stan/modules/multistate/transformed_parameters.stan` — at each transition's scaling block (lines 44-54 for 01, similar blocks for 02/12_s/12_t/03/32), replace:
```stan
if (enable_ms_level_baseline_hazard[lv] == 1) {
  log_lambda_gp_01_level_intercept[lv_start:lv_end] =
    raw_log_lambda_gp_01_level_intercept[lv_start:lv_end] *
    fe_log_lambda_gp_01_level_intercept_sd[lv];
} else {
  log_lambda_gp_01_level_intercept[lv_start:lv_end] =
    raw_log_lambda_gp_01_level_intercept[lv_start:lv_end] *
    log_lambda_gp_01_level_intercept_sd[lv];
}
```
with:
```stan
int mode = enable_ms_level_baseline_hazard[lv];
if (mode == LEVEL_MODE_RE_CP) {
  int c_lo = cp_level_pos_ms_baseline_shared[lv];
  int c_hi = cp_level_pos_ms_baseline_shared[lv + 1] - 1;
  log_lambda_gp_01_level_intercept[lv_start:lv_end] =
    cp_log_lambda_gp_01_level_intercept[c_lo:c_hi];
} else {
  int r_lo = raw_level_pos_ms_baseline_shared[lv];
  int r_hi = raw_level_pos_ms_baseline_shared[lv + 1] - 1;
  real scale = (mode == LEVEL_MODE_FE)
    ? fe_log_lambda_gp_01_level_intercept_sd[lv]
    : log_lambda_gp_01_level_intercept_sd[lv];
  log_lambda_gp_01_level_intercept[lv_start:lv_end] =
    raw_log_lambda_gp_01_level_intercept[r_lo:r_hi] * scale;
}
```

Apply to all 6 transitions. **Note:** `lv_start`/`lv_end` here are from `enabled_level_pos_ms_baseline` (the combined enabled-groups positions) — these stay unchanged on the left-hand side. The right-hand side reads from the per-bucket positions.

- [ ] **Step 7: Stan syntax check**

Run:
```bash
~/.cmdstan/cmdstan-2.38.0/bin/stanc \
  --include-paths=stan --include-paths=stan/psa \
  stan/psa/pioneer.stan
```
Expected: exits with code 0.

- [ ] **Step 8: Commit**

```bash
git add stan/multistate.stanfunctions stan/modules/multistate/
git commit -m "stan(ms): split baseline intercept per transition into _raw_+_cp_ buckets for RE_CP"
```

---

## Task 10: Extend `multistate` module — slope per transition (01, 02, 12)

Slope slices in `multistate` are smaller (only 3 transitions have slopes: 01, 02, 12). They're controlled by `enable_ms_level_cov` (the 0/1 gate), with parameterization following the intercept mode.

**Files:**
- Modify: `stan/modules/multistate/transformed_data.stan` (slope bucket counts)
- Modify: `stan/modules/multistate/parameters.stan:33,63,101` (split slope matrices)
- Modify: `stan/modules/multistate/priors.stan` (slope CP branches)
- Modify: `stan/modules/multistate/transformed_parameters.stan` (slope glue per transition)

- [ ] **Step 1: Add slope bucket counts to transformed_data**

Append to `stan/modules/multistate/transformed_data.stan`:
```stan
// --- Slope bucket routing (same pattern as intercepts) ---
// Slope parameterization follows the intercept mode at that level; slope is
// only included when enable_ms_level_cov[lv] == 1.
array[n_levels] int ms_slope_mode;
for (lv in 1:n_levels) {
  ms_slope_mode[lv] = enable_ms_level_cov[lv] ? enable_ms_level_baseline_hazard[lv] : 0;
}
int n_raw_groups_ms_slope_shared;
int n_cp_groups_ms_slope_shared;
array[n_levels + 1] int raw_level_pos_ms_slope_shared;
array[n_levels + 1] int cp_level_pos_ms_slope_shared;
(n_raw_groups_ms_slope_shared, raw_level_pos_ms_slope_shared,
 n_cp_groups_ms_slope_shared,  cp_level_pos_ms_slope_shared) =
  split_cp_ncp_pos(n_levels, n_forecast_groups_per_level, ms_slope_mode);
```

- [ ] **Step 2: Split slope matrices**

Edit `stan/modules/multistate/parameters.stan` — for transitions 01, 02, 12:
```stan
matrix[enable_ms_01 ? n_raw_groups_ms_slope_shared : 0, n_time_invariant_covar] raw_level_slope_01;
matrix[enable_ms_01 ? n_cp_groups_ms_slope_shared  : 0, n_time_invariant_covar] cp_level_slope_01;
// similarly for _02, _12
```

- [ ] **Step 3: Update slope priors — add CP branches**

Edit `stan/modules/multistate/priors.stan` — for each slope block (01, 02, 12), replace the existing slope prior with the same pattern used in Task 6.3 (read mode from `enable_ms_level_baseline_hazard`, route to raw or cp bucket).

- [ ] **Step 4: Wire slope glue in transformed_parameters**

Edit `stan/modules/multistate/transformed_parameters.stan` — for each transition's slope block, apply the same glue pattern from Task 6.4.

- [ ] **Step 5: Stan syntax check**

Run:
```bash
~/.cmdstan/cmdstan-2.38.0/bin/stanc \
  --include-paths=stan --include-paths=stan/psa \
  stan/psa/pioneer.stan
```
Expected: exits with code 0.

- [ ] **Step 6: Commit**

```bash
git add stan/modules/multistate/
git commit -m "stan(ms): split slope per transition into _raw_+_cp_ buckets for RE_CP"
```

---

## Task 11: Add `re_cp` to targets lookup table

**Files:**
- Modify: `targets/pioneer_targets.R:70-75`
- Modify: `targets/sclc_targets.R` (if it has a similar table)

- [ ] **Step 1: Add re_cp = 4L to pioneer lookup**

Edit `targets/pioneer_targets.R`:
```r
level_intercept_mode <- c(
  none  = 0L,  # No intercept at this level
  fe    = 1L,  # Fixed effect: SD is a data hyperparameter (no pooling)
  re    = 2L,  # Random effect, non-centered: SD estimated via NCP
  re_gp = 3L,  # Random effect + full GP residual (ms module only, for now)
  re_cp = 4L   # Random effect, centered: direct ~normal(0, sd) sampling
)
```

- [ ] **Step 2: Check sclc targets for a similar lookup**

Run: `grep -n "level_intercept_mode\b" targets/sclc_targets.R`

If a similar `level_intercept_mode` list exists, add `re_cp = 4L` to it. If not (sclc uses inline integer modes like `trial_level` / `patient_level`), skip — the integers will accept `4L` without changes since the Stan flag is now `upper=4`.

- [ ] **Step 3: Commit**

```bash
git add targets/pioneer_targets.R
# and sclc_targets.R if modified
git commit -m "targets: add re_cp=4L to level_intercept_mode lookup"
```

---

## Task 12: Update R initializers — route RE_CP draws to _cp_

Initializers build the per-chain starting-value list. Groups at `RE_CP` levels must be initialized into the new `_cp_` fields, not `_raw_`.

**Files:**
- Modify: `r/pioneer/initializers.R:30-35,49-63,94-99`
- Modify: `r/sclc/initializers.R:28-33,49-80,112-116`
- Modify: `r/sclc/initializers_fixed.R:193-197,222-226`

- [ ] **Step 1: Pioneer initializer — split draws**

Edit `r/pioneer/initializers.R` — replace lines 30-35 (enabled group counts) with bucket counts:
```r
# Route each level to raw (FE/RE/RE_GP) or cp (RE_CP) bucket.
route_level <- function(enable_vec) {
  raw_mask <- enable_vec %in% c(1L, 2L, 3L)
  cp_mask  <- enable_vec == 4L
  list(raw_mask = raw_mask, cp_mask = cp_mask)
}
tr_intercept_route   <- route_level(enable_level_intercept_tr)
frac_intercept_route <- route_level(enable_level_intercept_frac)
init_intercept_route <- route_level(enable_level_intercept_init)

n_raw_groups_tr_intercept   <- sum(n_forecast_groups_per_level[tr_intercept_route$raw_mask])
n_cp_groups_tr_intercept    <- sum(n_forecast_groups_per_level[tr_intercept_route$cp_mask])
n_raw_groups_frac_intercept <- sum(n_forecast_groups_per_level[frac_intercept_route$raw_mask])
n_cp_groups_frac_intercept  <- sum(n_forecast_groups_per_level[frac_intercept_route$cp_mask])
n_raw_groups_init_intercept <- sum(n_forecast_groups_per_level[init_intercept_route$raw_mask])
n_cp_groups_init_intercept  <- sum(n_forecast_groups_per_level[init_intercept_route$cp_mask])

# slope: mode follows intercept when enable_level_cov == 1, else 0
tr_slope_route   <- route_level(enable_level_intercept_tr   * (enable_level_cov_tr   == 1L))
frac_slope_route <- route_level(enable_level_intercept_frac * (enable_level_cov_frac == 1L))
init_slope_route <- route_level(enable_level_intercept_init * (enable_level_cov_init == 1L))

n_raw_groups_tr_slope   <- sum(n_forecast_groups_per_level[tr_slope_route$raw_mask])
n_cp_groups_tr_slope    <- sum(n_forecast_groups_per_level[tr_slope_route$cp_mask])
# ... similarly for frac_slope and init_slope
```

- [ ] **Step 2: Pioneer initializer — build raw + cp draws**

Replace the existing `tr_raw_level` / `frac_raw_level` / `init_raw_level` blocks with a helper that builds both buckets:
```r
build_level_draws <- function(enable_vec, sd_level, n_forecast_groups_per_level, n_levels) {
  raw_draws <- list()
  cp_draws  <- list()
  for (lv in seq_len(n_levels)) {
    mode <- enable_vec[lv]
    if (mode == 0L) next
    n <- n_forecast_groups_per_level[lv]
    draws <- if (lv == n_levels) rep(0, n) else rnorm(n, sd = 0.2)
    if (mode == 4L) {
      # CP: init on natural scale (multiply raw by sd_level)
      cp_draws[[length(cp_draws) + 1]] <- draws * sd_level[lv]
    } else {
      raw_draws[[length(raw_draws) + 1]] <- draws
    }
  }
  list(raw = unlist(raw_draws) %||% numeric(0),
       cp  = unlist(cp_draws)  %||% numeric(0))
}

tr_intercept_draws   <- build_level_draws(enable_level_intercept_tr,   tr_sd_level,   n_forecast_groups_per_level, n_levels)
frac_intercept_draws <- build_level_draws(enable_level_intercept_frac, frac_sd_level, n_forecast_groups_per_level, n_levels)
init_intercept_draws <- build_level_draws(enable_level_intercept_init, init_sd_level, n_forecast_groups_per_level, n_levels)
```

- [ ] **Step 3: Pioneer initializer — emit new field names**

Replace lines 94-99 (the `biomarker_init` list assignments for `_raw_level_intercept`) with:
```r
tr_sd_level_intercept_raw   = as.array(tr_sd_level[enable_level_intercept_tr   %in% c(2L, 4L)]),
tr_raw_level_intercept      = as.array(tr_intercept_draws$raw),
tr_cp_level_intercept       = as.array(tr_intercept_draws$cp),
frac_sd_level_intercept_raw = as.array(frac_sd_level[enable_level_intercept_frac %in% c(2L, 4L)]),
frac_raw_level_intercept    = as.array(frac_intercept_draws$raw),
frac_cp_level_intercept     = as.array(frac_intercept_draws$cp),
init_sd_level_intercept_raw = as.array(init_sd_level[enable_level_intercept_init %in% c(2L, 4L)]),
init_raw_level_intercept    = as.array(init_intercept_draws$raw),
init_cp_level_intercept     = as.array(init_intercept_draws$cp),
```

- [ ] **Step 4: Same mechanical refactor for slope matrices**

Replace the `tr_raw_level_slope` / `frac_raw_level_slope` / `init_raw_level_slope` matrix construction (lines 77-82) with per-bucket builds:
```r
tr_raw_level_slope <- matrix(rnorm(n_raw_groups_tr_slope * n_covar, sd = 0.5),
                              nrow = n_raw_groups_tr_slope, ncol = n_covar)
tr_cp_level_slope  <- matrix(rnorm(n_cp_groups_tr_slope * n_covar, sd = 0.1),
                              nrow = n_cp_groups_tr_slope, ncol = n_covar)
# similarly for frac and init
```

And emit them in `biomarker_init`:
```r
tr_raw_level_slope  = if (n_covar > 0) tr_raw_level_slope,
tr_cp_level_slope   = if (n_covar > 0) tr_cp_level_slope,
# ...
```

- [ ] **Step 5: Same for sclc initializers**

Apply the analogous changes to `r/sclc/initializers.R` (single-bucket) and `r/sclc/initializers_fixed.R` (which uses `enable_level_intercept_tr != 0` masking — extend to recognize mode 4).

- [ ] **Step 6: Commit**

```bash
git add r/pioneer/initializers.R r/sclc/initializers.R r/sclc/initializers_fixed.R
git commit -m "R: initializers route RE_CP groups to _cp_ fields, draws on natural scale"
```

---

## Task 13: Update R priors helper for RE_CP

Ensure `get_pioneer_priors()` in `r/pioneer/priors.R` leaves FE SD logic intact but is aware that `RE_CP` uses the same free-SD prior as `RE` (no change needed, but add a defensive comment).

**Files:**
- Modify: `r/pioneer/priors.R:28-34`

- [ ] **Step 1: Add a comment clarifying RE_CP uses RE path**

Edit `r/pioneer/priors.R` — around line 28, update the FE comment:
```r
# FE SD hyperparameters: non-zero only for levels using FE mode (mode == 1L).
# RE (mode 2), RE_GP (mode 3), and RE_CP (mode 4) all use the free SD parameter
# `*_sd_level_intercept_raw` with prior `~normal(0, *_sd_level_intercept_sd[lv])`.
# No new hyperparameter fields are needed for RE_CP.
fe_mode <- 1L
tr_fe_sd   <- rep(0, n_levels)
frac_fe_sd <- rep(0, n_levels)
init_fe_sd <- rep(0, n_levels)
tr_fe_sd[stan_data$enable_level_intercept_tr     == fe_mode] <- 0.50
frac_fe_sd[stan_data$enable_level_intercept_frac == fe_mode] <- 0.35
init_fe_sd[stan_data$enable_level_intercept_init == fe_mode] <- 0.50
```

- [ ] **Step 2: Commit**

```bash
git add r/pioneer/priors.R
git commit -m "R(pioneer): document RE_CP reuses RE's free SD prior"
```

---

## Task 14: Bit-exact regression test against main

The most important guarantee of this PR: when no level is set to `RE_CP`, the model must produce identical posteriors. Use a small fixture and hash the first 10 draws.

**Files:**
- Create: `tests/testthat/test-pioneer-bit-exact.R`

- [ ] **Step 1: Build a minimal fixture**

Run locally (before the PR changes):
```r
# Use a tiny pioneer fixture (first 5 trial patients + 5 RWD patients, 1 chain, 50 warmup, 50 sampling, seed=42)
# Save the draws array to tests/testthat/fixtures/pioneer-pre-re-cp-baseline.rds
```

If the fixture doesn't exist yet, create it on `main` (before the PR) and commit it as part of Task 14.

- [ ] **Step 2: Write the regression test**

Create `tests/testthat/test-pioneer-bit-exact.R`:
```r
library(testthat)
library(here)
library(posterior)

test_that("RE_CP scaffolding is bit-exact when no level uses RE_CP", {
  skip_on_ci()  # Needs cmdstan; run locally only
  skip_if_not(file.exists(here("tests", "testthat", "fixtures", "pioneer-pre-re-cp-baseline.rds")),
              "Baseline fixture missing — regenerate on main")

  baseline <- readRDS(here("tests", "testthat", "fixtures", "pioneer-pre-re-cp-baseline.rds"))

  # Run the same fit with the new code path (RE_CP added but not used)
  # Expect: all draws identical to within 1e-10 (seeded MCMC is deterministic modulo
  # floating-point associativity; the reordered loop should not perturb anything)
  current <- run_tiny_pioneer_fit(seed = 42L)
  for (param in intersect(variables(baseline), variables(current))) {
    expect_equal(
      as.numeric(subset_draws(current, variable = param)),
      as.numeric(subset_draws(baseline, variable = param)),
      tolerance = 1e-10,
      info = paste("Parameter drift:", param)
    )
  }
})
```

- [ ] **Step 3: Run locally on the new branch**

Run: `Rscript -e 'testthat::test_file("tests/testthat/test-pioneer-bit-exact.R")'`
Expected: PASS. If it fails, one of the `transformed_parameters` edits reordered a computation — debug before merging.

- [ ] **Step 4: Commit**

```bash
git add tests/testthat/test-pioneer-bit-exact.R tests/testthat/fixtures/pioneer-pre-re-cp-baseline.rds
git commit -m "test: bit-exact regression baseline for RE_CP scaffolding PR"
```

---

## Task 15: Stan compilation check across all entry-point models

Every top-level Stan file that includes the touched modules must compile.

**Files:** (verification only, no edits)

- [ ] **Step 1: Run stanc on every entry-point model**

```bash
for model in \
  stan/psa/pioneer.stan \
  stan/tumor/sf-ssm-log-space.stan \
  stan/tumor/sf-ssls-lfo.stan \
  stan/ms-standalone.stan; do
  echo "=== $model ==="
  ~/.cmdstan/cmdstan-2.38.0/bin/stanc \
    --include-paths=stan \
    --include-paths=stan/psa \
    --include-paths=stan/tumor \
    "$model" || exit 1
done
```

Expected: every model parses cleanly (empty stdout, exit code 0).

- [ ] **Step 2: Run full R test suite**

```bash
Rscript -e 'testthat::test_dir("tests/testthat")'
```
Expected: all pass (new tests + existing tests).

---

## Task 16: Run a small Domino job to verify end-to-end

**Files:** (verification only)

- [ ] **Step 1: Push the branch**

```bash
git push origin karim/speedup
```

- [ ] **Step 2: Launch a small test job via /start-job**

Invoke the `/start-job` skill:
- target: `pioneer_fit_posterior_combined_markov_fe_psa_re_ms`
- store: `main` (or a scratch store if available)
- hardware: small CPU tier (whatever the previous #638 job used)
- goal: verify the new Stan compiles and fits with all-NCP defaults — no RE_CP anywhere yet.

- [ ] **Step 3: Check completion**

Monitor via `pipeline-monitor` agent. Expect posterior means within MC error of the pre-PR run #638.

---

## Self-Review Checklist

After completing all tasks, verify:

1. **Spec coverage:** Every module in scope (tr, frac, init, ms) has both intercept and slope slices extended. Check by grepping `_cp_level_intercept` and `_cp_level_slope` across `stan/modules/` — expect 4 × 2 = 8 matches minimum (some ms slices are per-transition).

2. **Bit-exact default:** No task changes the model's output when all levels remain at their existing modes. If Task 14 fails, a transformed_parameters edit accidentally changed evaluation order — fix before merging.

3. **FE preservation:** Task 7, 8, 9, 10 each include an `if (mode == LEVEL_MODE_FE)` dispatch in the prior block (raw-path branch includes FE) and transformed_parameters (scale uses `*_fe_sd_level_intercept`). FE groups never appear in `_cp_` bucket.

4. **Type consistency:** Every `_cp_` vector declared in `parameters.stan` is referenced in `priors.stan` and `transformed_parameters.stan`. Every bucket-count scalar declared in `transformed_data.stan` is used in parameter sizing.

5. **Initializer fields match parameter names:** `tr_cp_level_intercept` in R initializer ↔ `tr_cp_level_intercept` in Stan. Any mismatch → cmdstan error at init time.

---

## Subsequent PR (separate) — flip level-3 for the pioneer combined variant

Not part of this plan. In a follow-up:

1. In `targets/pioneer_targets.R`, change `enable_level_intercept_init` for the `combined_markov_fe_psa_re_ms` variant:
   ```r
   enable_level_intercept_init = as.array(c(
     trial   = level_intercept_mode[["none"]],
     arm     = arm_level,
     patient = level_intercept_mode[["re_cp"]]   # was "re"
   )),
   ```
2. Refit via `/start-job`.
3. Re-run diagnostic targets (`pioneer_energy_corr`, `pioneer_group_ess`).
4. Compare E-BFMI against the pre-PR run #638: target E-BFMI ≥ 0.3 across all 4 chains, `|cor(E, init_sd_level_intercept[3])|` < 0.3.
5. Write the ebfmi-diagnosis report (existing pending Task #4 from memory).
