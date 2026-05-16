# Multi-level SD Hierarchy — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Extend the pioneer multilevel parameter framework so the SD scaling the `tr` module's location-level raw draws can itself vary across coarser levels via a 2D sub-hierarchy. Ship the primitive with bit-exactness on the inactive path and a minimal synthetic exercise of the active path.

**Architecture:** Add a 2D mode matrix (`enable_sd_level_intercept_mode_tr[L, ℓ]`) and a self-contained sub-hierarchy block (`stan/modules/tr/`) that assembles a per-group `tr_sd_intercept_pergroup[lv]` vector. Existing location-side assembly in `stan/modules/tr/transformed_parameters.stan:47-48`, `stan/modules/tr/priors.stan:47-51`, and `stan/psa/pioneer.stan:156` is touched with a single substitution (`scalar → per-group lookup`). All new code paths are gated so that an empty mode tibble produces bit-identical draws to `main`.

**Tech Stack:** Stan (cmdstan 2.38 with `#include` modules), R (tidyverse, targets, cmdstanr), testthat.

**Spec:** `docs/superpowers/specs/2026-04-24-multilevel-sd-hierarchy-design.md`.

**Out of scope for this plan** (tracked separately): activating the feature on any existing variant, running Domino fits, comparing convergence to the #109 baseline, or interpreting posterior summaries of per-group log-SDs. Those belong to diagnostic follow-up work once the primitive is in place. The plan's acceptance gate is: bit-exact on inactive path, passing unit tests, passing Stan compile, and one synthetic end-to-end test that proves the active path executes and returns per-group SDs.

---

## File Structure

### Files to create
- `tests/testthat/test-validate-sd-modes.R` — R-side validator unit tests
- `tests/testthat/test-stan-split-sd-cp-ncp-pos.R` — Stan helper test harness (standalone model approach, mirrors `test-stan-split-cp-ncp-pos.R`)
- `tests/testthat/test-sd-hierarchy-bit-exact.R` — Phase 1's bit-exact regression gate
- `stan/modules/tr/_sd_subhierarchy_parameters.stan` — new parameters for the sub-hierarchy
- `stan/modules/tr/_sd_subhierarchy_transformed_parameters.stan` — assembly of `tr_sd_intercept_pergroup[lv]`
- `stan/modules/tr/_sd_subhierarchy_priors.stan` — priors for new sub-hierarchy parameters

### Files to modify
- `stan/hierarchy.stanfunctions` — add `split_sd_cp_ncp_pos` helper (appended)
- `stan/modules/tr/flags.stan` — declare `enable_sd_level_intercept_mode_tr[n_levels, n_levels]` and hyperparams
- `stan/modules/tr/hyperparams.stan` — declare new hyperparameter scalars (`tr_log_sd_level_intercept_pop_mean`, etc.)
- `stan/modules/tr/transformed_data.stan` — compute `has_sd_subhierarchy_tr_intercept[lv]`, position arrays, and the `subgroup_idx_tr_L_to_ell[lv, sub_lv]` flat index via `split_sd_cp_ncp_pos`
- `stan/modules/tr/parameters.stan` — `#include` the new `_sd_subhierarchy_parameters.stan`
- `stan/modules/tr/transformed_parameters.stan` — insert per-group SD assembly at the top; substitute `tr_sd_level_intercept[lv]` → `tr_sd_intercept_pergroup[lv]` at the existing scaling line (line 48)
- `stan/modules/tr/priors.stan` — `#include` the new `_sd_subhierarchy_priors.stan`; substitute the scalar in the RE_CP branch (lines 47, 50)
- `stan/psa/pioneer.stan` — substitute the scalar in the analogous scaling line (line 156)
- `r/multi_level_hierarchy.R` — add `validate_sd_modes()`
- `r/priors.R` — add three new defaults (`tr_log_sd_level_intercept_pop_mean`, `tr_log_sd_level_intercept_pop_sd`, `tr_sd_hyperscale_level_intercept_sd`) and pass-through for `tr_sd_intercept_modes`
- `r/pioneer/initializers.R` — route initial draws for new parameters
- `r/pioneer/prepare_analysis_data.R` — accept `tr_sd_intercept_modes` tibble, validate it, serialize the 2D mode matrix into the Stan data list (both in `prepare_*_stan_data()` paths)
- `targets/pioneer_targets.R` — add `~tr_sd_intercept_modes` column to the outer tribble at line 375; add a new variant row for Phase 2

### Files NOT touched (Phase 3 scope)
- `stan/modules/frac/`, `stan/modules/init/`, `stan/modules/multistate/` — SD hierarchy added here in Phase 3
- Any tumor-only model (`stan/tumor/sf-ssm-log-space.stan`) — not exercised in pioneer

---

## Task 1: `split_sd_cp_ncp_pos` Stan helper — write test harness

**Files:**
- Create: `tests/testthat/test-stan-split-sd-cp-ncp-pos.R`

The pattern follows the existing `tests/testthat/test-stan-split-cp-ncp-pos.R` — build a tiny Stan model that exposes the helper's outputs via a `generated quantities` block, compile it once per test, and assert on the generated values.

- [ ] **Step 1: Read the existing test file for pattern reuse**

Run: `Read tests/testthat/test-stan-split-cp-ncp-pos.R` (full file). Extract the cmdstanr compilation pattern and the per-test `tar_tempfile` convention used there.

- [ ] **Step 2: Write the failing test file**

```r
# tests/testthat/test-stan-split-sd-cp-ncp-pos.R
# Tests the split_sd_cp_ncp_pos helper in stan/hierarchy.stanfunctions.
# Pattern mirrors test-stan-split-cp-ncp-pos.R (the 1D split helper).

library(testthat)
library(cmdstanr)

# A tiny Stan model that loads hierarchy.stanfunctions and exposes
# the helper's outputs via generated quantities.
build_harness_stan_code <- function() '
functions {
  #include "hierarchy.stanfunctions"
}
data {
  int<lower=1> n_levels;
  array[n_levels] int<lower=0> n_groups_per_level;
  array[n_levels, n_levels] int sd_mode;
}
generated quantities {
  int raw_total;
  array[n_levels, n_levels + 1] int raw_pos;
  int cp_total;
  array[n_levels, n_levels + 1] int cp_pos;
  (raw_total, raw_pos, cp_total, cp_pos) =
    split_sd_cp_ncp_pos(n_levels, n_groups_per_level, sd_mode);
}
'

# Compile once and reuse across tests.
compile_harness <- function() {
  f <- tempfile(fileext = ".stan")
  writeLines(build_harness_stan_code(), f)
  cmdstan_model(f, include_paths = "stan")
}

run_harness <- function(model, n_levels, n_groups_per_level, sd_mode) {
  fit <- model$generate_quantities(
    data = list(
      n_levels = n_levels,
      n_groups_per_level = n_groups_per_level,
      sd_mode = sd_mode
    ),
    fitted_params = posterior::as_draws(data.frame(lp__ = 0)),
    parallel_chains = 1
  )
  list(
    raw_total = fit$draws("raw_total", format = "draws_df")$raw_total[1],
    cp_total  = fit$draws("cp_total",  format = "draws_df")$cp_total[1],
    raw_pos   = fit$draws("raw_pos",   format = "draws_matrix"),
    cp_pos    = fit$draws("cp_pos",    format = "draws_matrix")
  )
}

test_that("all-zero mode matrix yields zero totals and identity-position arrays", {
  model <- compile_harness()
  sd_mode <- matrix(0L, 3, 3)
  res <- run_harness(model, 3, c(2L, 3L, 10L), sd_mode)
  expect_equal(res$raw_total, 0)
  expect_equal(res$cp_total, 0)
})

test_that("single RE_CP entry goes into cp bucket only", {
  model <- compile_harness()
  sd_mode <- matrix(0L, 3, 3)
  sd_mode[3, 2] <- 4L  # (L=patient, ℓ=arm) = RE_CP
  res <- run_harness(model, 3, c(2L, 3L, 10L), sd_mode)
  expect_equal(res$raw_total, 0)
  expect_equal(res$cp_total, 3)  # n_groups at ℓ=arm = 3
})

test_that("mixed FE/RE/RE_CP partitions groups correctly across buckets", {
  model <- compile_harness()
  sd_mode <- matrix(0L, 3, 3)
  sd_mode[3, 1] <- 2L  # (patient, trial) = RE  → raw bucket, 2 groups
  sd_mode[3, 2] <- 4L  # (patient, arm)   = RE_CP → cp bucket, 3 groups
  sd_mode[2, 1] <- 1L  # (arm, trial)     = FE   → raw bucket, 2 groups
  res <- run_harness(model, 3, c(2L, 3L, 10L), sd_mode)
  expect_equal(res$raw_total, 4)  # 2 (RE at patient/trial) + 2 (FE at arm/trial)
  expect_equal(res$cp_total, 3)   # 3 (RE_CP at patient/arm)
})

test_that("RE_GP (mode=3) triggers fatal_error", {
  model <- compile_harness()
  sd_mode <- matrix(0L, 3, 3)
  sd_mode[3, 2] <- 3L  # RE_GP not allowed for sub-hierarchy
  expect_error(
    run_harness(model, 3, c(2L, 3L, 10L), sd_mode),
    regexp = "RE_GP|mode 3|fatal"
  )
})

test_that("upper-triangular entries (ℓ >= L) that are non-zero trigger fatal_error", {
  model <- compile_harness()
  sd_mode <- matrix(0L, 3, 3)
  sd_mode[1, 2] <- 4L  # L=trial cannot have sub-level ℓ=arm
  expect_error(
    run_harness(model, 3, c(2L, 3L, 10L), sd_mode),
    regexp = "sub_lv|ℓ|fatal"
  )
})
```

- [ ] **Step 3: Run the tests — expect them to FAIL (function not defined yet)**

Run: `Rscript -e 'testthat::test_file("tests/testthat/test-stan-split-sd-cp-ncp-pos.R")'`
Expected: Compilation fails with `Identifier 'split_sd_cp_ncp_pos' not in scope` (or equivalent Stan syntax error).

- [ ] **Step 4: Commit**

```bash
git add tests/testthat/test-stan-split-sd-cp-ncp-pos.R
git commit -m "test(stan): failing harness for split_sd_cp_ncp_pos helper"
```

---

## Task 2: `split_sd_cp_ncp_pos` Stan helper — implement

**Files:**
- Modify: `stan/hierarchy.stanfunctions` (append; current file is 125 lines, ends at the `}` of `split_cp_ncp_pos`)

- [ ] **Step 1: Append the helper**

Append to `stan/hierarchy.stanfunctions`:

```stan

/**
 * Split sub-hierarchy (L, ℓ) pairs into raw and cp buckets by mode.
 * For each location level L and sub-level ℓ ∈ 1..L-1:
 *   - modes in {FE=1, RE=2} → raw bucket
 *   - mode == RE_CP (4)      → cp bucket
 *   - mode == NONE (0)       → skipped (no storage)
 * RE_GP (3) is rejected via fatal_error.
 * Any non-zero entry at (L, ℓ) with ℓ >= L is also rejected.
 *
 * @param n_levels            number of levels in the stack
 * @param n_groups_per_level  group counts per level (used to size buckets)
 * @param sd_mode             [n_levels, n_levels] matrix; entries (L, ℓ) with ℓ >= L must be 0
 *
 * @return tuple(
 *   int raw_total,
 *   array[n_levels, n_levels + 1] int raw_pos,
 *   int cp_total,
 *   array[n_levels, n_levels + 1] int cp_pos
 * )
 *   raw_pos[L, ℓ+1] - raw_pos[L, ℓ] = count of groups for location-level-L, sub-level-ℓ in raw bucket
 *   (likewise for cp_pos).
 */
tuple(int, array[,] int, int, array[,] int) split_sd_cp_ncp_pos(
  int n_levels,
  array[] int n_groups_per_level,
  array[,] int sd_mode
) {
  // Validate entries
  for (L in 1:n_levels) {
    for (sub_lv in 1:n_levels) {
      int m = sd_mode[L, sub_lv];
      if (m < 0 || m > 4)
        fatal_error("sd_mode[", L, ",", sub_lv, "] must be 0..4 (NONE/FE/RE/RE_GP/RE_CP); got ", m);
      if (m == 3)
        fatal_error("sd_mode[", L, ",", sub_lv, "] = RE_GP (3) is not allowed on the sub-hierarchy");
      if (m != 0 && sub_lv >= L)
        fatal_error(
          "sd_mode[", L, ",", sub_lv, "] non-zero but sub_lv (", sub_lv,
          ") >= L (", L, "); sub-levels must be positionally before L"
        );
    }
  }

  // Build raw_pos and cp_pos as cumulative-offset arrays:
  // raw_pos[L, ℓ]  = starting flat index for (L, ℓ) block in raw bucket
  // raw_pos[L, n_levels+1] = one past end of L's contribution in raw bucket
  array[n_levels, n_levels + 1] int raw_pos;
  array[n_levels, n_levels + 1] int cp_pos;

  int raw_cursor = 1;
  int cp_cursor  = 1;
  for (L in 1:n_levels) {
    for (sub_lv in 1:n_levels) {
      raw_pos[L, sub_lv] = raw_cursor;
      cp_pos[L, sub_lv]  = cp_cursor;
      int m = sd_mode[L, sub_lv];
      if (m == 1 || m == 2) {
        raw_cursor += n_groups_per_level[sub_lv];
      } else if (m == 4) {
        cp_cursor += n_groups_per_level[sub_lv];
      }
    }
    raw_pos[L, n_levels + 1] = raw_cursor;
    cp_pos[L, n_levels + 1]  = cp_cursor;
  }

  int raw_total = raw_cursor - 1;
  int cp_total  = cp_cursor - 1;
  return (raw_total, raw_pos, cp_total, cp_pos);
}
```

- [ ] **Step 2: Run the tests — expect PASS**

Run: `Rscript -e 'testthat::test_file("tests/testthat/test-stan-split-sd-cp-ncp-pos.R")'`
Expected: All 5 tests pass.

- [ ] **Step 3: Commit**

```bash
git add stan/hierarchy.stanfunctions
git commit -m "feat(stan): split_sd_cp_ncp_pos helper for SD sub-hierarchy bucketing"
```

---

## Task 3: `validate_sd_modes` R validator — write test

**Files:**
- Create: `tests/testthat/test-validate-sd-modes.R`

- [ ] **Step 1: Write the failing test**

```r
# tests/testthat/test-validate-sd-modes.R
library(testthat)

source("../../r/multi_level_hierarchy.R")

test_that("empty tibble passes", {
  sd_modes <- tibble::tibble(location_level = character(), sub_level = character(), mode = character())
  expect_silent(validate_sd_modes(sd_modes, level_stack = c("trial", "arm", "patient")))
})

test_that("valid rows pass", {
  sd_modes <- tibble::tribble(
    ~location_level, ~sub_level, ~mode,
    "patient",       "arm",      "re_cp",
    "patient",       "trial",    "re",
    "arm",           "trial",    "fe"
  )
  expect_silent(validate_sd_modes(sd_modes, level_stack = c("trial", "arm", "patient")))
})

test_that("unknown location_level errors", {
  sd_modes <- tibble::tibble(location_level = "galaxy", sub_level = "arm", mode = "re_cp")
  expect_error(
    validate_sd_modes(sd_modes, level_stack = c("trial", "arm", "patient")),
    regexp = "galaxy|unknown|level_stack"
  )
})

test_that("unknown sub_level errors", {
  sd_modes <- tibble::tibble(location_level = "patient", sub_level = "region", mode = "re_cp")
  expect_error(
    validate_sd_modes(sd_modes, level_stack = c("trial", "arm", "patient")),
    regexp = "region|unknown|level_stack"
  )
})

test_that("sub_level not positionally before location_level errors", {
  sd_modes <- tibble::tibble(location_level = "trial", sub_level = "arm", mode = "re")
  expect_error(
    validate_sd_modes(sd_modes, level_stack = c("trial", "arm", "patient")),
    regexp = "position|before|nested"
  )
})

test_that("sub_level equal to location_level errors", {
  sd_modes <- tibble::tibble(location_level = "patient", sub_level = "patient", mode = "re_cp")
  expect_error(
    validate_sd_modes(sd_modes, level_stack = c("trial", "arm", "patient")),
    regexp = "position|before"
  )
})

test_that("disallowed mode 'gp' errors", {
  sd_modes <- tibble::tibble(location_level = "patient", sub_level = "arm", mode = "gp")
  expect_error(
    validate_sd_modes(sd_modes, level_stack = c("trial", "arm", "patient")),
    regexp = "mode|gp|allowed"
  )
})

test_that("typo mode 'random' errors", {
  sd_modes <- tibble::tibble(location_level = "patient", sub_level = "arm", mode = "random")
  expect_error(
    validate_sd_modes(sd_modes, level_stack = c("trial", "arm", "patient")),
    regexp = "mode|random|allowed"
  )
})

test_that("error message names the offending row", {
  sd_modes <- tibble::tibble(
    location_level = c("patient", "trial"),
    sub_level      = c("arm",     "arm"),
    mode           = c("re_cp",   "re")
  )
  expect_error(
    validate_sd_modes(sd_modes, level_stack = c("trial", "arm", "patient")),
    regexp = "row 2|\\(trial, arm\\)"
  )
})
```

- [ ] **Step 2: Run — expect FAIL (function not defined)**

Run: `Rscript -e 'testthat::test_file("tests/testthat/test-validate-sd-modes.R")'`
Expected: FAIL with `could not find function "validate_sd_modes"`.

- [ ] **Step 3: Commit**

```bash
git add tests/testthat/test-validate-sd-modes.R
git commit -m "test(r): failing tests for validate_sd_modes"
```

---

## Task 4: `validate_sd_modes` R validator — implement

**Files:**
- Modify: `r/multi_level_hierarchy.R` (append after `create_hierarchy_structure`, keeping the file's documentation style)

- [ ] **Step 1: Append the validator**

Add to the end of `r/multi_level_hierarchy.R`:

```r

#' Validate an SD sub-hierarchy mode tibble against a level stack
#'
#' Checks each row of `sd_modes` against the current implementation's positional
#' (nested) rule: `sub_level` must appear positionally before `location_level` in
#' the `level_stack`. Empty tibbles pass trivially.
#'
#' @param sd_modes Tibble with columns `location_level`, `sub_level`, `mode`.
#'   An empty tibble (zero rows) is valid and represents the all-NONE configuration.
#' @param level_stack Character vector naming the levels in positional order,
#'   e.g. `c("trial", "arm", "patient")`.
#' @param allowed_modes Character vector of permitted mode strings.
#'   Defaults to `c("none", "fe", "re", "re_cp")`. `"gp"` is NOT in this default
#'   because the sub-hierarchy does not support RE_GP.
#'
#' @return Invisibly returns `sd_modes` if all rows are valid; stops otherwise.
validate_sd_modes <- function(sd_modes,
                              level_stack,
                              allowed_modes = c("none", "fe", "re", "re_cp")) {
  stopifnot(is.data.frame(sd_modes))
  if (nrow(sd_modes) == 0) return(invisible(sd_modes))

  required_cols <- c("location_level", "sub_level", "mode")
  missing <- setdiff(required_cols, names(sd_modes))
  if (length(missing) > 0) {
    stop("validate_sd_modes: missing columns ", paste(missing, collapse = ", "))
  }

  level_pos <- setNames(seq_along(level_stack), level_stack)

  for (i in seq_len(nrow(sd_modes))) {
    row <- sd_modes[i, ]
    loc <- row$location_level
    sub <- row$sub_level
    mode <- row$mode

    if (!(loc %in% level_stack)) {
      stop(sprintf(
        "validate_sd_modes: row %d (%s, %s): unknown location_level %s; not in level_stack = (%s)",
        i, loc, sub, loc, paste(level_stack, collapse = ", ")
      ))
    }
    if (!(sub %in% level_stack)) {
      stop(sprintf(
        "validate_sd_modes: row %d (%s, %s): unknown sub_level %s; not in level_stack = (%s)",
        i, loc, sub, sub, paste(level_stack, collapse = ", ")
      ))
    }
    if (level_pos[sub] >= level_pos[loc]) {
      stop(sprintf(
        "validate_sd_modes: row %d (%s, %s): sub_level must appear positionally before location_level in the (nested) level_stack",
        i, loc, sub
      ))
    }
    if (!(mode %in% allowed_modes)) {
      stop(sprintf(
        "validate_sd_modes: row %d (%s, %s): mode '%s' not in allowed_modes = (%s)",
        i, loc, sub, mode, paste(allowed_modes, collapse = ", ")
      ))
    }
  }

  invisible(sd_modes)
}
```

- [ ] **Step 2: Run — expect PASS**

Run: `Rscript -e 'testthat::test_file("tests/testthat/test-validate-sd-modes.R")'`
Expected: All 9 tests pass.

- [ ] **Step 3: Commit**

```bash
git add r/multi_level_hierarchy.R
git commit -m "feat(r): validate_sd_modes for SD sub-hierarchy config"
```

---

## Task 5: Stan declarations for `tr` module sub-hierarchy

**Files:**
- Modify: `stan/modules/tr/flags.stan`
- Modify: `stan/modules/tr/hyperparams.stan`
- Create: `stan/modules/tr/_sd_subhierarchy_parameters.stan`

This task adds all declarations. Nothing is computed or used yet — so syntax must still parse with the existing main model.

- [ ] **Step 1: Add mode matrix flag**

Add to the end of `stan/modules/tr/flags.stan`:

```stan

// ===== SD sub-hierarchy (issue #110) =====
// Mode matrix: [location-level L, sub-level ℓ]
// Entries with ℓ >= L must be 0 (validated in transformed_data via split_sd_cp_ncp_pos).
// 0=NONE, 1=FE, 2=RE, 4=RE_CP. (RE_GP=3 is rejected.)
array[n_levels, n_levels] int<lower=0, upper=4> enable_sd_level_intercept_mode_tr;
```

- [ ] **Step 2: Add hyperparameter scalars**

Add to the end of `stan/modules/tr/hyperparams.stan`:

```stan

// ===== SD sub-hierarchy hyperparameters (issue #110) =====
// Population log-SD prior (per location level L)
array[n_levels] real tr_log_sd_level_intercept_pop_mean;
array[n_levels] real<lower=0> tr_log_sd_level_intercept_pop_sd;
// Hyperscale prior (half-normal via <lower=0>)
array[n_levels, n_levels] real<lower=0> tr_sd_hyperscale_level_intercept_sd;
// FE hyperscale values used when mode == FE (data-supplied)
array[n_levels, n_levels] real<lower=0> tr_fe_sd_hyperscale_level_intercept;
```

- [ ] **Step 3: Create the new parameters block**

Create `stan/modules/tr/_sd_subhierarchy_parameters.stan`:

```stan
// tr/_sd_subhierarchy_parameters.stan
// Parameters for the SD sub-hierarchy (issue #110).
// All sizes are zero when the mode matrix is all-NONE → no new params allocated → bit-exact.

// Population log-SD per location level (one entry per level with an active sub-hierarchy)
vector[n_subhier_active_tr_intercept] tr_log_sd_level_intercept_pop;

// Free hyperscales for RE / RE_CP sub-levels
vector<lower=0>[n_sd_hyperscales_tr_intercept] tr_sd_hyperscale_level_intercept_raw;

// Flat raw bucket for FE/RE sub-levels (indexed by raw_pos_tr_log_sd_intercept)
vector[n_raw_groups_tr_log_sd_intercept] tr_raw_log_sd_level_intercept;

// Flat cp bucket for RE_CP sub-levels (indexed by cp_pos_tr_log_sd_intercept)
vector[n_cp_groups_tr_log_sd_intercept] tr_cp_log_sd_level_intercept;
```

- [ ] **Step 4: `#include` it in tr/parameters.stan**

Append one line to `stan/modules/tr/parameters.stan` (after line 49):

```stan

#include "modules/tr/_sd_subhierarchy_parameters.stan"
```

- [ ] **Step 5: Verify Stan still compiles (all new code still inert — sizes not yet defined)**

Expected: compilation WILL fail with unresolved identifiers (`n_subhier_active_tr_intercept`, `n_sd_hyperscales_tr_intercept`, etc.). This is expected — the sizes are computed in the next task. Do not commit yet; proceed to Task 6.

---

## Task 6: Stan transformed_data — wire size scalars, position arrays, and subgroup index

**Files:**
- Modify: `stan/modules/tr/transformed_data.stan` (current file size: check with `wc -l`)

- [ ] **Step 1: Append sub-hierarchy bookkeeping to transformed_data**

Append to `stan/modules/tr/transformed_data.stan`:

```stan

// ===== SD sub-hierarchy sizing (issue #110) =====
// All derived from enable_sd_level_intercept_mode_tr (declared in flags.stan).

// Per-level active flag: sum over row L > 0.
array[n_levels] int has_sd_subhierarchy_tr_intercept;
int n_subhier_active_tr_intercept = 0;
for (lv in 1:n_levels) {
  int row_sum = 0;
  for (sub_lv in 1:n_levels) {
    row_sum += enable_sd_level_intercept_mode_tr[lv, sub_lv];
  }
  has_sd_subhierarchy_tr_intercept[lv] = (row_sum > 0) ? 1 : 0;
  n_subhier_active_tr_intercept += has_sd_subhierarchy_tr_intercept[lv];
}

// Position arrays via split_sd_cp_ncp_pos. Total counts size the raw/cp buckets.
int n_raw_groups_tr_log_sd_intercept;
array[n_levels, n_levels + 1] int raw_pos_tr_log_sd_intercept;
int n_cp_groups_tr_log_sd_intercept;
array[n_levels, n_levels + 1] int cp_pos_tr_log_sd_intercept;
(
  n_raw_groups_tr_log_sd_intercept,
  raw_pos_tr_log_sd_intercept,
  n_cp_groups_tr_log_sd_intercept,
  cp_pos_tr_log_sd_intercept
) = split_sd_cp_ncp_pos(
  n_levels,
  n_forecast_groups_per_level,
  enable_sd_level_intercept_mode_tr
);

// Hyperscale count: one per (L, ℓ) with mode ∈ {RE, RE_CP} (not FE — FE uses data value).
int n_sd_hyperscales_tr_intercept = 0;
for (lv in 1:n_levels) {
  for (sub_lv in 1:(lv - 1)) {
    int m = enable_sd_level_intercept_mode_tr[lv, sub_lv];
    if (m == LEVEL_MODE_RE || m == LEVEL_MODE_RE_CP) {
      n_sd_hyperscales_tr_intercept += 1;
    }
  }
}

// Flat hyperscale index: hyperscale_idx_tr_intercept[L, ℓ] = 1-based position in
// tr_sd_hyperscale_level_intercept_raw, or 0 if not a RE/RE_CP entry.
array[n_levels, n_levels] int hyperscale_idx_tr_intercept;
{
  int idx = 0;
  for (lv in 1:n_levels) {
    for (sub_lv in 1:n_levels) {
      int m = enable_sd_level_intercept_mode_tr[lv, sub_lv];
      if (m == LEVEL_MODE_RE || m == LEVEL_MODE_RE_CP) {
        idx += 1;
        hyperscale_idx_tr_intercept[lv, sub_lv] = idx;
      } else {
        hyperscale_idx_tr_intercept[lv, sub_lv] = 0;
      }
    }
  }
}

// Pop-vector index: pop_idx_tr_intercept[L] = 1-based position in
// tr_log_sd_level_intercept_pop, or 0 if L has no active sub-hierarchy.
array[n_levels] int pop_idx_tr_intercept;
{
  int idx = 0;
  for (lv in 1:n_levels) {
    if (has_sd_subhierarchy_tr_intercept[lv]) {
      idx += 1;
      pop_idx_tr_intercept[lv] = idx;
    } else {
      pop_idx_tr_intercept[lv] = 0;
    }
  }
}

// Subgroup flat index: for each location-level-lv group g, which sub-level-ℓ group
// does it map to? Derived from patient_level_groups by picking a representative
// patient per lv-group (well-defined because the implementation is nested).
// Flat layout: subgroup_idx_tr_intercept_flat[lv][sub_lv][g_lv] but stored as a
// jagged array via enabled_level_pos-style indexing.
//
// Simplification for Phase 1: store only the lookup we need — for each location
// level lv with active sub-hierarchy and each active sub-level ℓ, a vector of
// length n_forecast_groups_per_level[lv] giving the ℓ-group id for each lv-group.
array[n_levels, n_levels] vector[max(n_forecast_groups_per_level)] subgroup_idx_tr_intercept;
for (lv in 1:n_levels) {
  if (!has_sd_subhierarchy_tr_intercept[lv]) continue;
  int ngroups_lv = n_forecast_groups_per_level[lv];
  for (sub_lv in 1:(lv - 1)) {
    if (enable_sd_level_intercept_mode_tr[lv, sub_lv] == LEVEL_MODE_NONE) continue;
    // Build a representative-patient map from lv-group → ℓ-group.
    // patient_level_groups[i, lv] = group id at level lv for patient i.
    // For each lv-group g, find the first patient with that group id and read its ℓ-group.
    array[ngroups_lv] int seen = zeros_int_array(ngroups_lv);
    for (i in 1:n_patients) {
      int g_lv = patient_level_groups[i, lv];
      if (g_lv >= 1 && g_lv <= ngroups_lv && seen[g_lv] == 0) {
        subgroup_idx_tr_intercept[lv, sub_lv][g_lv] = patient_level_groups[i, sub_lv];
        seen[g_lv] = 1;
      }
    }
  }
}
```

- [ ] **Step 2: Verify Stan compiles**

Run: `~/.cmdstan/cmdstan-2.38.0/bin/stanc --include-paths=stan --include-paths=stan/psa stan/psa/pioneer.stan`
Expected: Exit 0 with no errors. The `tr_log_sd_level_intercept_pop` / `raw`/`cp` parameters are declared but not yet used — Stan permits that.

- [ ] **Step 3: Commit**

```bash
git add stan/modules/tr/flags.stan stan/modules/tr/hyperparams.stan \
        stan/modules/tr/parameters.stan stan/modules/tr/_sd_subhierarchy_parameters.stan \
        stan/modules/tr/transformed_data.stan
git commit -m "stan(tr): declarations + transformed_data bookkeeping for SD sub-hierarchy"
```

---

## Task 7: Stan transformed_parameters — build `tr_sd_intercept_pergroup` and substitute

**Files:**
- Create: `stan/modules/tr/_sd_subhierarchy_transformed_parameters.stan`
- Modify: `stan/modules/tr/transformed_parameters.stan` (insert `#include` + one line substitution at line 48)

- [ ] **Step 1: Create the self-contained sub-hierarchy transformed_parameters block**

Create `stan/modules/tr/_sd_subhierarchy_transformed_parameters.stan`:

```stan
// tr/_sd_subhierarchy_transformed_parameters.stan
// Build per-group SD vector tr_sd_intercept_pergroup[lv] for each level.
// When sub-hierarchy inactive at lv: constant fill with tr_sd_level_intercept[lv] → bit-exact.
// When active: sum log-deviations additively, then exp once.

array[n_levels] vector[max(n_forecast_groups_per_level)] tr_sd_intercept_pergroup;
for (lv in 1:n_levels) {
  int ngroups_lv = n_forecast_groups_per_level[lv];
  if (!has_sd_subhierarchy_tr_intercept[lv]) {
    // Inactive: constant fill with existing scalar → bit-exact with pre-feature code.
    tr_sd_intercept_pergroup[lv][1:ngroups_lv] =
      rep_vector(tr_sd_level_intercept[lv], ngroups_lv);
  } else {
    // Active: additive log-scale assembly.
    vector[ngroups_lv] log_sd_g = rep_vector(
      tr_log_sd_level_intercept_pop[pop_idx_tr_intercept[lv]],
      ngroups_lv
    );
    for (sub_lv in 1:(lv - 1)) {
      int m = enable_sd_level_intercept_mode_tr[lv, sub_lv];
      if (m == LEVEL_MODE_NONE) continue;

      int n_sub_groups = n_forecast_groups_per_level[sub_lv];
      vector[n_sub_groups] delta;

      if (m == LEVEL_MODE_RE_CP) {
        int c_lo = cp_pos_tr_log_sd_intercept[lv, sub_lv];
        int c_hi = cp_pos_tr_log_sd_intercept[lv, sub_lv + 1] - 1;
        delta = tr_cp_log_sd_level_intercept[c_lo:c_hi];
      } else {
        // FE or RE
        int r_lo = raw_pos_tr_log_sd_intercept[lv, sub_lv];
        int r_hi = raw_pos_tr_log_sd_intercept[lv, sub_lv + 1] - 1;
        real hyperscale = (m == LEVEL_MODE_FE)
          ? tr_fe_sd_hyperscale_level_intercept[lv, sub_lv]
          : tr_sd_hyperscale_level_intercept_raw[hyperscale_idx_tr_intercept[lv, sub_lv]];
        delta = hyperscale * tr_raw_log_sd_level_intercept[r_lo:r_hi];
      }

      // Scatter delta from sub-level-ℓ groups to location-level-L groups via subgroup_idx.
      for (g_lv in 1:ngroups_lv) {
        int g_sub = to_int(subgroup_idx_tr_intercept[lv, sub_lv][g_lv]);
        log_sd_g[g_lv] += delta[g_sub];
      }
    }
    tr_sd_intercept_pergroup[lv][1:ngroups_lv] = exp(log_sd_g);
  }
}
```

- [ ] **Step 2: Insert `#include` at the top of tr/transformed_parameters.stan**

Modify `stan/modules/tr/transformed_parameters.stan`. At line 28 (immediately after the closing `}` of the SD expansion block at line 27), insert:

```stan

// SD sub-hierarchy: builds tr_sd_intercept_pergroup[lv] (bit-exact fill when inactive)
#include "modules/tr/_sd_subhierarchy_transformed_parameters.stan"
```

- [ ] **Step 3: Substitute the scalar SD at line 48 of tr/transformed_parameters.stan**

Find the line (was line 48 before the insertion):

```stan
      tr_sd_level_intercept[lv] * tr_raw_level_intercept[r_lo:r_hi];
```

Replace with:

```stan
      tr_sd_intercept_pergroup[lv][1:(r_hi - r_lo + 1)] .* tr_raw_level_intercept[r_lo:r_hi];
```

- [ ] **Step 4: Substitute the scalar SD in `stan/psa/pioneer.stan` at line 156**

Find:

```stan
      tr_sd_level_intercept[lv] * tr_raw_level_intercept[lv_start_param:lv_end_param];
```

Replace with:

```stan
      tr_sd_intercept_pergroup[lv][1:(lv_end_param - lv_start_param + 1)] .* tr_raw_level_intercept[lv_start_param:lv_end_param];
```

- [ ] **Step 5: Verify Stan compiles**

Run: `~/.cmdstan/cmdstan-2.38.0/bin/stanc --include-paths=stan --include-paths=stan/psa stan/psa/pioneer.stan`
Expected: Exit 0. The code compiles; the new block is inert when all modes NONE.

- [ ] **Step 6: Commit**

```bash
git add stan/modules/tr/_sd_subhierarchy_transformed_parameters.stan \
        stan/modules/tr/transformed_parameters.stan stan/psa/pioneer.stan
git commit -m "stan(tr): build tr_sd_intercept_pergroup and substitute in location assembly"
```

---

## Task 8: Stan priors — substitute scalar in RE_CP branch and add sub-hierarchy priors

**Files:**
- Create: `stan/modules/tr/_sd_subhierarchy_priors.stan`
- Modify: `stan/modules/tr/priors.stan` (substitute lines 47, 50 — the RE_CP branch — and `#include` the new block)

- [ ] **Step 1: Create the sub-hierarchy priors block**

Create `stan/modules/tr/_sd_subhierarchy_priors.stan`:

```stan
// tr/_sd_subhierarchy_priors.stan
// Priors for SD sub-hierarchy parameters (issue #110).
// All statements are vacuously true when the mode matrix is all-NONE (size-0 vectors).

// Population log-SD per location level
for (lv in 1:n_levels) {
  if (has_sd_subhierarchy_tr_intercept[lv]) {
    int idx = pop_idx_tr_intercept[lv];
    tr_log_sd_level_intercept_pop[idx]
      ~ normal(tr_log_sd_level_intercept_pop_mean[lv], tr_log_sd_level_intercept_pop_sd[lv]);
  }
}

// Free hyperscales (half-normal via <lower=0>)
for (lv in 1:n_levels) {
  for (sub_lv in 1:(lv - 1)) {
    int m = enable_sd_level_intercept_mode_tr[lv, sub_lv];
    if (m == LEVEL_MODE_RE || m == LEVEL_MODE_RE_CP) {
      int idx = hyperscale_idx_tr_intercept[lv, sub_lv];
      tr_sd_hyperscale_level_intercept_raw[idx]
        ~ normal(0, tr_sd_hyperscale_level_intercept_sd[lv, sub_lv]);
    }
  }
}

// NCP raw deviations (FE and RE share raw bucket)
if (n_raw_groups_tr_log_sd_intercept > 0) {
  tr_raw_log_sd_level_intercept ~ std_normal();
}

// CP deviations — sampled at hyperscale directly
for (lv in 1:n_levels) {
  if (!has_sd_subhierarchy_tr_intercept[lv]) continue;
  for (sub_lv in 1:(lv - 1)) {
    if (enable_sd_level_intercept_mode_tr[lv, sub_lv] == LEVEL_MODE_RE_CP) {
      int c_lo = cp_pos_tr_log_sd_intercept[lv, sub_lv];
      int c_hi = cp_pos_tr_log_sd_intercept[lv, sub_lv + 1] - 1;
      real hs = tr_sd_hyperscale_level_intercept_raw[hyperscale_idx_tr_intercept[lv, sub_lv]];
      tr_cp_log_sd_level_intercept[c_lo:c_hi] ~ normal(0, hs);
    }
  }
}
```

- [ ] **Step 2: Substitute scalar SD in RE_CP branch of tr/priors.stan**

At `stan/modules/tr/priors.stan` line 47 (inside the RE_CP branch, student_t case):

Find:
```stan
          tr_cp_level_intercept[c_lo:c_hi]
            ~ student_t(tr_nu_level[lv], 0, tr_sd_level_intercept[lv]);
```

Replace with:
```stan
          tr_cp_level_intercept[c_lo:c_hi]
            ~ student_t(tr_nu_level[lv], 0, tr_sd_intercept_pergroup[lv][1:(c_hi - c_lo + 1)]);
```

At line 50 (normal case):
Find:
```stan
          tr_cp_level_intercept[c_lo:c_hi]
            ~ normal(0, tr_sd_level_intercept[lv]);
```

Replace with:
```stan
          tr_cp_level_intercept[c_lo:c_hi]
            ~ normal(0, tr_sd_intercept_pergroup[lv][1:(c_hi - c_lo + 1)]);
```

- [ ] **Step 3: Include the new priors block at the end of tr/priors.stan**

Append to `stan/modules/tr/priors.stan`:

```stan

// SD sub-hierarchy priors (inert when mode matrix is all-NONE)
#include "modules/tr/_sd_subhierarchy_priors.stan"
```

- [ ] **Step 4: Verify Stan compiles**

Run: `~/.cmdstan/cmdstan-2.38.0/bin/stanc --include-paths=stan --include-paths=stan/psa stan/psa/pioneer.stan`
Expected: Exit 0.

- [ ] **Step 5: Commit**

```bash
git add stan/modules/tr/priors.stan stan/modules/tr/_sd_subhierarchy_priors.stan
git commit -m "stan(tr): priors for SD sub-hierarchy; substitute per-group lookup in RE_CP branch"
```

---

## Task 9: R-side data wiring — priors, initializers, prepare_* helpers

**Files:**
- Modify: `r/priors.R` (add three new scalar defaults and the `tr_sd_intercept_modes` passthrough; search for the pioneer prior block around line 170)
- Modify: `r/pioneer/initializers.R` (add init routing for the four new parameter families)
- Modify: `r/pioneer/prepare_analysis_data.R` (accept `tr_sd_intercept_modes` in `prepare_psa_standalone_stan_data` and in the full-model prep; validate; serialize a 2D integer matrix into the Stan data list)

- [ ] **Step 1: Add prior defaults to r/priors.R**

Find the block that defines pioneer prior defaults (search for existing `tr_sd_level_intercept_sd` or similar). Add alongside:

```r
  # SD sub-hierarchy (issue #110): defaults produce a lognormal prior on sd_L
  # whose marginal approximately matches the current half-normal when hyperscales are small.
  # Only used when a non-empty mode tibble is supplied; otherwise not exercised.
  tr_log_sd_level_intercept_pop_mean = rep(log(0.5), n_levels),
  tr_log_sd_level_intercept_pop_sd   = rep(0.5, n_levels),
  tr_sd_hyperscale_level_intercept_sd = matrix(0.3, n_levels, n_levels),
  tr_fe_sd_hyperscale_level_intercept = matrix(0.3, n_levels, n_levels),
```

(Adjust `log(0.5)` if the current half-normal scale default differs — grep for `tr_sd_level_intercept_sd` in `r/priors.R` to match.)

- [ ] **Step 2: Add initializer routing in r/pioneer/initializers.R**

In the initializer function (`create_pioneer_initializer`), after the existing `tr_sd_level_intercept_raw` init logic, add:

```r
  # SD sub-hierarchy init (issue #110). Zero-sized when mode matrix is all-NONE.
  n_subhier_active_tr <- sum(rowSums(stan_data$enable_sd_level_intercept_mode_tr) > 0)
  n_hyperscales_tr <- sum(
    stan_data$enable_sd_level_intercept_mode_tr == 2L |
    stan_data$enable_sd_level_intercept_mode_tr == 4L
  )
  init$tr_log_sd_level_intercept_pop <- rep(log(0.5), n_subhier_active_tr)
  init$tr_sd_hyperscale_level_intercept_raw <- rep(0.1, n_hyperscales_tr)
  # Raw/cp buckets: their sizes are computed Stan-side; initializer just needs zero vectors
  # of matching length. Read them from stan_data if your initializer infrastructure already
  # passes these counts; otherwise compute with split_sd_cp_ncp_pos-equivalent in R:
  init$tr_raw_log_sd_level_intercept <- rep(0, stan_data$n_raw_groups_tr_log_sd_intercept %||% 0L)
  init$tr_cp_log_sd_level_intercept  <- rep(0, stan_data$n_cp_groups_tr_log_sd_intercept %||% 0L)
```

(The R function should compute `n_raw_groups_tr_log_sd_intercept` / `n_cp_groups_tr_log_sd_intercept` and include them in `stan_data` — see the next step.)

- [ ] **Step 3: Serialize the mode matrix in prepare_*_stan_data helpers**

In `r/pioneer/prepare_analysis_data.R`, find `prepare_psa_standalone_stan_data(...)` at line 948. Modify its signature to accept an `sd_modes` argument:

```r
prepare_psa_standalone_stan_data <- function(base_stan_data, fit_data, tr_sd_intercept_modes = NULL) {
```

At the top of the function body, validate and serialize:

```r
  # Default to empty (all-NONE) if not provided
  if (is.null(tr_sd_intercept_modes)) {
    tr_sd_intercept_modes <- tibble::tibble(
      location_level = character(),
      sub_level = character(),
      mode = character()
    )
  }
  level_stack <- colnames(base_stan_data$patient_level_groups)
  validate_sd_modes(tr_sd_intercept_modes, level_stack)

  # Build n_levels x n_levels integer mode matrix
  mode_codes <- c(none = 0L, fe = 1L, re = 2L, re_gp = 3L, re_cp = 4L)
  n_levels <- length(level_stack)
  enable_sd_level_intercept_mode_tr <- matrix(0L, n_levels, n_levels)
  for (i in seq_len(nrow(tr_sd_intercept_modes))) {
    row <- tr_sd_intercept_modes[i, ]
    L_idx   <- which(level_stack == row$location_level)
    sub_idx <- which(level_stack == row$sub_level)
    enable_sd_level_intercept_mode_tr[L_idx, sub_idx] <- mode_codes[[row$mode]]
  }
```

Add both the matrix and derived counts to the returned stan_data list:

```r
  stan_data$enable_sd_level_intercept_mode_tr <- enable_sd_level_intercept_mode_tr

  # Precompute raw/cp bucket sizes so the initializer has them (mirrors Stan's split_sd_cp_ncp_pos)
  n_groups_per_level <- attr(base_stan_data$patient_level_groups, "n_forecast_groups_per_level") %||%
    apply(base_stan_data$patient_level_groups, 2, \(x) length(unique(x)))
  stan_data$n_raw_groups_tr_log_sd_intercept <- sum(
    vapply(seq_len(n_levels), \(L) {
      sum(n_groups_per_level[seq_len(L - 1)] *
          (enable_sd_level_intercept_mode_tr[L, seq_len(L - 1)] %in% c(1L, 2L)))
    }, integer(1))
  )
  stan_data$n_cp_groups_tr_log_sd_intercept <- sum(
    vapply(seq_len(n_levels), \(L) {
      sum(n_groups_per_level[seq_len(L - 1)] *
          (enable_sd_level_intercept_mode_tr[L, seq_len(L - 1)] == 4L))
    }, integer(1))
  )
```

Apply the same serialization to any other `prepare_*_stan_data(...)` helpers that build pioneer stan_data. Grep for `enable_level_intercept_tr` in `r/pioneer/` to find them; add the mode-matrix serialization alongside.

- [ ] **Step 4: Confirm R side loads cleanly**

Run: `Rscript -e 'source("r/pioneer/prepare_analysis_data.R"); source("r/pioneer/initializers.R"); source("r/priors.R"); message("R sources load")'`
Expected: no errors.

- [ ] **Step 5: Commit**

```bash
git add r/priors.R r/pioneer/initializers.R r/pioneer/prepare_analysis_data.R
git commit -m "r(pioneer): wire SD sub-hierarchy through priors, initializer, stan_data prep"
```

---

## Task 10: Bit-exact regression test

**Files:**
- Create: `tests/testthat/test-sd-hierarchy-bit-exact.R`

This is the primary Phase 1 gate. The test runs a small pioneer `psa_standalone` fit with an empty mode tibble, then compares the Stan CSV draws against a frozen reference from `main` (the parent commit).

- [ ] **Step 1: Capture the reference CSVs on `main`**

Run (manually, once, before implementing the test):

```bash
# On a separate throwaway worktree anchored at main:
git worktree add /tmp/multilevel-sd-reference main
cd /tmp/multilevel-sd-reference
# Run a tiny deterministic fit (iter_warmup=50, iter_sampling=50, seed=42, 1 chain).
# Example target name:  pioneer_fit_posterior_combined_markov_fe_psa_re_ms
# Save the output CSV to tests/testthat/fixtures/sd-hierarchy/bitexact-reference.csv
mkdir -p /mnt/code/.worktrees/karim/multilevel-sd/tests/testthat/fixtures/sd-hierarchy
cp <output-dir>/<fit>-1.csv \
  /mnt/code/.worktrees/karim/multilevel-sd/tests/testthat/fixtures/sd-hierarchy/bitexact-reference.csv
cd - && git worktree remove /tmp/multilevel-sd-reference
```

(If this is impractical in a single session, the test can compare against any stored reference — document the reference's provenance in a sibling README.)

- [ ] **Step 2: Write the failing test**

```r
# tests/testthat/test-sd-hierarchy-bit-exact.R
# Phase 1 gate: with empty mode tibble, CSV output must match the main-branch reference byte-for-byte.
library(testthat)
library(cmdstanr)

test_that("empty mode tibble produces bit-identical draws to main reference", {
  skip_if(
    !file.exists("fixtures/sd-hierarchy/bitexact-reference.csv"),
    "no reference CSV"
  )

  # Recompile the model on this branch
  model <- cmdstan_model(
    "../../stan/psa/pioneer.stan",
    include_paths = c("../../stan", "../../stan/psa"),
    force_recompile = TRUE
  )

  # Build stan_data with EMPTY mode tibble (Phase 1 bit-exact case).
  source("../../r/pioneer/prepare_analysis_data.R")
  source("../../r/pioneer/initializers.R")
  # ... construct a minimal base_stan_data matching the reference run's data.
  # The exact construction depends on the fixture data; document it in
  # tests/testthat/fixtures/sd-hierarchy/README.md.

  stan_data <- prepare_psa_standalone_stan_data(base_stan_data, fit_data = TRUE)
  stopifnot(all(stan_data$enable_sd_level_intercept_mode_tr == 0))

  fit <- model$sample(
    data = stan_data,
    seed = 42,
    iter_warmup = 50,
    iter_sampling = 50,
    chains = 1,
    parallel_chains = 1,
    init = create_pioneer_initializer(stan_data),
    refresh = 0
  )

  # Byte-level diff against the reference
  produced <- fit$output_files()[1]
  reference <- "fixtures/sd-hierarchy/bitexact-reference.csv"
  expect_identical(
    tools::md5sum(produced),
    tools::md5sum(reference),
    info = "CSV differs — Phase 1 bit-exactness gate violated"
  )
})
```

- [ ] **Step 3: Run — expect PASS (or skip if reference not captured yet)**

Run: `Rscript -e 'testthat::test_file("tests/testthat/test-sd-hierarchy-bit-exact.R")'`
Expected: PASS (if reference captured) or SKIP with a clear message. If it FAILS, something in tasks 5–9 broke bit-exactness — diagnose via `diff` on the two CSVs and look for non-zero parameter additions or RNG-consumption shifts.

- [ ] **Step 4: Commit**

```bash
git add tests/testthat/test-sd-hierarchy-bit-exact.R \
        tests/testthat/fixtures/sd-hierarchy/
git commit -m "test: bit-exact regression gate for SD sub-hierarchy Phase 1"
```

---

## Task 11: Thread the mode tibble through `pioneer_targets.R`

**Files:**
- Modify: `targets/pioneer_targets.R` (tribble at lines 374-430+; and `prepare_psa_standalone_stan_data` call at line 1623)

- [ ] **Step 1: Add the `~tr_sd_intercept_modes` column to the outer tribble**

Locate line 374 in `targets/pioneer_targets.R`. The tribble header looks like:

```r
    tribble(
      ~model, ~hist, ~part_b_only, ~flatiron_only, ~flatiron_filtered, ~combined_unfiltered,
      ~model_formula, ~arm_level, ~patient_level, ~propensity, ~visit_gated_01,
      ~enable_02_tv_cov, ~n_tv_covar, ~latent_01, ~ms_time_scale_12, ~ms_arm_level,
      ~metric_files, ~standalone_metric_files, ~max_treedepth, ~fit_psa,
```

Add a new final column `~tr_sd_intercept_modes`:

```r
    tribble(
      ~model, ~hist, ~part_b_only, ~flatiron_only, ~flatiron_filtered, ~combined_unfiltered,
      ~model_formula, ~arm_level, ~patient_level, ~propensity, ~visit_gated_01,
      ~enable_02_tv_cov, ~n_tv_covar, ~latent_01, ~ms_time_scale_12, ~ms_arm_level,
      ~metric_files, ~standalone_metric_files, ~max_treedepth, ~fit_psa,
      ~tr_sd_intercept_modes,
```

For **every existing row** in the tribble, append an empty tibble to the row's trailing columns (representing all-NONE, bit-exact):

```r
      tibble::tibble(location_level = character(), sub_level = character(), mode = character()),
```

- [ ] **Step 2: Pass it through the PSA standalone tar_target at line ~1622**

At `targets/pioneer_targets.R` line 1622, modify:

```r
      tar_target(
        psa_standalone_stan_data,
        prepare_psa_standalone_stan_data(
          base_pioneer_stan_data,
          fit_data,
          tr_sd_intercept_modes = tr_sd_intercept_modes
        )
      ),
```

- [ ] **Step 3: Verify the pipeline parses (no actual build)**

Run: `Rscript -e 'source("targets/pioneer_targets.R"); message("Targets definition parses")'`
Expected: no errors. (If `source` has side effects that launch jobs, replace with `parse()` or use `targets::tar_validate()` after setting `TAR_PROJECT=pioneer`.)

- [ ] **Step 4: Commit**

```bash
git add targets/pioneer_targets.R
git commit -m "targets(pioneer): thread tr_sd_intercept_modes through psa_standalone tar_target"
```

---

## Task 12: Synthetic active-path test

**Files:**
- Create: `tests/testthat/test-sd-subhierarchy-active.R`

Proves the active code path runs end-to-end on a minimal synthetic input. This is a **structural / smoke test** — it does NOT assert on convergence quality or compare against the inactive path quantitatively. Its only job: "when the mode matrix is non-zero, parameters allocate, priors fire, the model samples without crashing, and `tr_sd_intercept_pergroup` varies across groups at the configured level."

- [ ] **Step 1: Write the test**

```r
# tests/testthat/test-sd-subhierarchy-active.R
# Smoke test: the SD sub-hierarchy active code path executes and produces a
# per-group SD vector that actually varies across groups.
# NOT a convergence test. NOT a diagnostic test. Structural only.

library(testthat)
library(cmdstanr)

test_that("active sub-hierarchy allocates parameters and samples without crashing", {
  # Build a minimal Stan harness that #includes hierarchy.stanfunctions and exposes
  # a tiny version of the tr assembly. To avoid the full pioneer fixture, run
  # a dry 5-iteration sample directly on stan/psa/pioneer.stan using the
  # smallest valid pioneer data list and the active mode matrix.
  #
  # If a suitable tiny fixture is not yet committed, SKIP with a TODO.
  skip_if(
    !file.exists("fixtures/sd-hierarchy/minimal-active-stan-data.rds"),
    "minimal active-path fixture not yet captured"
  )

  stan_data <- readRDS("fixtures/sd-hierarchy/minimal-active-stan-data.rds")

  # Confirm the fixture has an active mode matrix (patient/arm = RE_CP).
  expect_true(stan_data$enable_sd_level_intercept_mode_tr[3, 2] == 4L)

  model <- cmdstan_model(
    "../../stan/psa/pioneer.stan",
    include_paths = c("../../stan", "../../stan/psa")
  )

  fit <- model$sample(
    data = stan_data,
    seed = 1,
    iter_warmup = 5,
    iter_sampling = 5,
    chains = 1,
    parallel_chains = 1,
    refresh = 0
  )

  # Structural assertions:
  draws <- posterior::as_draws_df(fit$draws(
    variables = c("tr_log_sd_level_intercept_pop",
                  "tr_cp_log_sd_level_intercept",
                  "tr_sd_intercept_pergroup")
  ))

  # 1. Population log-SD is present (size > 0 because one level is active).
  expect_gt(length(grep("^tr_log_sd_level_intercept_pop", names(draws))), 0)

  # 2. CP bucket is present and has length = n_groups at sub_level (arm = 3 in fixture).
  cp_cols <- grep("^tr_cp_log_sd_level_intercept", names(draws), value = TRUE)
  expect_gt(length(cp_cols), 0)

  # 3. Per-group SD for the active level varies across groups within a draw.
  pg_cols <- grep("^tr_sd_intercept_pergroup\\[3,", names(draws), value = TRUE)
  expect_gt(length(pg_cols), 1L)
  row1 <- as.numeric(draws[1, pg_cols])
  expect_gt(stats::sd(row1), 0,
            info = "per-group SD should vary across groups when sub-hierarchy is active")
})
```

- [ ] **Step 2: Capture the minimal active-path fixture (one-time)**

If the fixture file does not yet exist, create it by running a minimal pioneer
`base_stan_data` through `prepare_psa_standalone_stan_data()` with an active mode
tibble:

```r
# Helper script (run once to save the fixture):
source("r/pioneer/prepare_analysis_data.R")
# base_stan_data is whatever minimal pioneer fixture already lives in the repo
# (or can be built from a small synthetic patient table — see tests for pioneer
# prep helpers for the pattern).
stan_data <- prepare_psa_standalone_stan_data(
  base_stan_data, fit_data = TRUE,
  tr_sd_intercept_modes = tibble::tibble(
    location_level = "patient", sub_level = "arm", mode = "re_cp"
  )
)
dir.create("tests/testthat/fixtures/sd-hierarchy", showWarnings = FALSE, recursive = TRUE)
saveRDS(stan_data, "tests/testthat/fixtures/sd-hierarchy/minimal-active-stan-data.rds")
```

If no minimal pioneer fixture exists in the repo, the test can `skip` for now;
the fixture can be captured in a follow-up PR. The plan's primary gate is the
bit-exact test in Task 10, not this smoke test.

- [ ] **Step 3: Run the test**

Run: `Rscript -e 'testthat::test_file("tests/testthat/test-sd-subhierarchy-active.R")'`
Expected: PASS (if fixture captured) or SKIP with the fixture-missing message.

- [ ] **Step 4: Commit**

```bash
git add tests/testthat/test-sd-subhierarchy-active.R
# If fixture was captured:
git add tests/testthat/fixtures/sd-hierarchy/
git commit -m "test: synthetic active-path smoke test for SD sub-hierarchy"
```

---

## Task 13: Close out the branch

**Files:** (none modified)

- [ ] **Step 1: Verify all tests pass**

Run: `Rscript -e 'testthat::test_dir("tests/testthat")'`
Expected: all tests pass (or skip with documented reason).

- [ ] **Step 2: Verify Stan compiles cleanly**

Run: `~/.cmdstan/cmdstan-2.38.0/bin/stanc --include-paths=stan --include-paths=stan/psa stan/psa/pioneer.stan`
Expected: Exit 0.

- [ ] **Step 3: Invoke the finishing skill**

Use the `superpowers:finishing-a-development-branch` skill to decide whether to
merge to `main`, open a PR, or leave for review.

---

## Self-review checklist

- [x] Every spec §4–§10 structural requirement maps to a task: helper (T1–T2), validator (T3–T4), Stan declarations (T5–T6), transformed_parameters (T7), priors (T8), R-side wiring (T9), bit-exact gate (T10), targets plumbing (T11), synthetic active-path smoke test (T12), branch close-out (T13).
- [x] Out of scope per user direction: Phase 2 activation variant, Domino job launch, convergence diagnostics, #109 baseline comparison, posterior analysis. These are explicitly excluded in the Goal section and belong to a separate diagnostic follow-up plan.
- [x] No placeholders except the reference-CSV capture step (Task 10) and the minimal active-path fixture (Task 12) — both are one-time operational setups with explicit instructions and gracefully degrade via `skip_if` when the fixture is absent.
- [x] Type consistency: `enable_sd_level_intercept_mode_tr[n_levels, n_levels]` is declared in Task 5, sized in Task 6, read in Task 7 and Task 8, serialized from R in Task 9. `tr_sd_intercept_pergroup[lv]` is declared in Task 7 and consumed in Tasks 7 (transformed_parameters) and 8 (priors). `pop_idx_tr_intercept`, `hyperscale_idx_tr_intercept`, `cp_pos_tr_log_sd_intercept`, `raw_pos_tr_log_sd_intercept`, `subgroup_idx_tr_intercept`, `n_raw_groups_tr_log_sd_intercept`, `n_cp_groups_tr_log_sd_intercept`, `n_subhier_active_tr_intercept`, `n_sd_hyperscales_tr_intercept`, `has_sd_subhierarchy_tr_intercept[lv]` are all defined in Task 6 and consumed in Tasks 7 and 8.
