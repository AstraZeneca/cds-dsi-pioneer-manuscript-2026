# Modular Initializer Refactor

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Extract shared multistate GP initialization into a reusable helper function so SCLC and pioneer initializers can't drift out of sync.

**Architecture:** Mirror the pattern already used in `r/initializers_fixed.R`, where `ms_init_values_fixed()` is a standalone helper called by both the joint model and standalone multistate initializers via `c(tumor_init, .ms_init(environment()))`. Create the equivalent `ms_init_values()` for the non-fixed (prior-sampling) initializers.

**Tech Stack:** R, rlang (`lst`), invgamma, purrr (`compact`)

---

## Current State

Three initializer files exist:

| File | Multistate init | Pattern |
|---|---|---|
| `r/initializers.R` (SCLC) | Inline in `lst()`, lines 572–705 | Monolithic |
| `r/pioneer/initializers.R` | Inline in `lst()`, lines 183–311 | Monolithic copy |
| `r/initializers_fixed.R` | `ms_init_values_fixed()` helper | **Modular** ← target pattern |

The duplication between `r/initializers.R` and `r/pioneer/initializers.R` is what caused the missing-inits bug. The `_fixed.R` file already solved this with a shared helper.

## Target State

```
r/initializers_ms.R              ← NEW: ms_init_values() shared helper
r/initializers.R                 ← MODIFIED: calls ms_init_values()
r/pioneer/initializers.R      ← MODIFIED: calls ms_init_values()
r/initializers_fixed.R           ← UNCHANGED (already modular)
```

Both SCLC and pioneer compose their init lists as:
```r
c(biomarker_init, .ms_init(environment()))
```

---

### Task 1: Create the shared multistate init helper

**Files:**
- Create: `r/initializers_ms.R`

**Step 1: Extract `ms_init_values()` from the SCLC initializer**

Create `r/initializers_ms.R` containing a single function `ms_init_values(env)` that:
- Takes an environment (from `with(stan_data, { ... })`) as input — same signature as `ms_init_values_fixed(env)`
- Contains the derived flags computation (`need_12_s_gp`, `need_12_t_gp`, enabled group counts)
- Contains all multistate GP parameter initialization (4 GP sets, covariate coefficients, random slopes)
- Returns a named list via `tibble::lst()` — NULL entries removed by caller via `compact()`

The body comes from `r/initializers.R` lines 428–705 (the derived flags block + the multistate section of the `lst()` return). The key difference from `ms_init_values_fixed()` is that this version samples from priors (`rnorm`, `rinvgamma`) rather than using fixed values.

```r
# Shared multistate initialization helper (prior-sampling version)
#
# Called by both the SCLC initializer (create_tumor_ssls_initializer) and
# the pioneer initializer (create_pioneer_initializer).
# Must be called inside a with(stan_data, { ... }) block so that all
# multistate fields (flags, hyperparameters, dimensions) are in scope.
#
# Compare with ms_init_values_fixed() in initializers_fixed.R which uses
# deterministic starting values instead of prior draws.

ms_init_values <- function(env) {
  with(env, {
    # Derived flags for 1→2 GPs
    need_12_s_gp <- enable_ms_12 && (ms_time_scale_12 == 1 || ms_time_scale_12 == 2)
    need_12_t_gp <- enable_ms_12 && (ms_time_scale_12 == 0 || ms_time_scale_12 == 2)

    # Multistate enabled group counts
    n_enabled_groups_ms_baseline_01 <- if (enable_ms_01) sum(n_groups_per_level[enable_ms_level_baseline_hazard == 1]) else 0L
    n_enabled_groups_ms_baseline_02 <- if (enable_ms_02) sum(n_groups_per_level[enable_ms_level_baseline_hazard == 1]) else 0L
    n_enabled_groups_ms_baseline_12_s <- if (need_12_s_gp) sum(n_groups_per_level[enable_ms_level_baseline_hazard == 1]) else 0L
    n_enabled_groups_ms_baseline_12_t <- if (need_12_t_gp) sum(n_groups_per_level[enable_ms_level_baseline_hazard == 1]) else 0L
    n_enabled_groups_ms_slope <- sum(n_groups_per_level[enable_ms_level_cov == 1])

    tibble::lst(
      # [All GP params for 0→1, 0→2, 1→2 sojourn, 1→2 clock-forward]
      # [Time-varying + time-invariant covariate coefficients]
      # [Multi-level random slope SDs and raw slopes]
      # ... exact code from r/initializers.R lines 580-705
    )
  })
}
```

**Step 2: Verify syntax**

Run: `Rscript -e 'parse("r/initializers_ms.R"); cat("OK\n")'`
Expected: `OK`

**Step 3: Commit**

```bash
git add r/initializers_ms.R
git commit -m "Extract ms_init_values() shared multistate init helper"
```

---

### Task 2: Refactor the SCLC initializer to use the shared helper

**Files:**
- Modify: `r/initializers.R` — `create_tumor_ssls_initializer()` (lines 405–709)

**Step 1: Add closure capture at top of `create_tumor_ssls_initializer()`**

After line 405 (`function(chain_id) {` is line 406), add:
```r
  .ms_init <- ms_init_values
```
This captures the helper in the closure, same pattern as `_fixed.R` line 133.

**Step 2: Remove inline multistate code**

Remove the following from `create_tumor_ssls_initializer()`:
- Lines 428–451: derived flags and enabled group counts (now inside `ms_init_values()`)
- Lines 572–705: the entire multistate `lst()` block

**Step 3: Split the `lst()` return into biomarker init + multistate composition**

Change the monolithic `lst(...)` into:
```r
      biomarker_init <- lst(
        # ... everything from tr_loc_pop through measure_sd_sld (lines 522-570)
      )

      c(biomarker_init, .ms_init(environment()))
```

This mirrors `_fixed.R` line 248: `c(tumor_init, .ms_init(environment()))`.

**Step 4: Verify syntax**

Run: `Rscript -e 'parse("r/initializers.R"); cat("OK\n")'`
Expected: `OK`

**Step 5: Commit**

```bash
git add r/initializers.R
git commit -m "Refactor SCLC initializer to use shared ms_init_values()"
```

---

### Task 3: Refactor the pioneer initializer to use the shared helper

**Files:**
- Modify: `r/pioneer/initializers.R` — `create_pioneer_initializer()` (lines 13–313)

**Step 1: Add closure capture at top of `create_pioneer_initializer()`**

After line 13, add:
```r
  .ms_init <- ms_init_values
```

**Step 2: Remove inline multistate code**

Remove from `create_pioneer_initializer()`:
- Lines 30–47: derived flags and enabled group counts
- Lines 183–311: the entire multistate block inside `lst()`

**Step 3: Split the `lst()` return into biomarker init + multistate composition**

```r
      biomarker_init <- lst(
        # ... everything from tr_loc_pop through measure_sd_psa (lines 110-181)
      )

      c(biomarker_init, .ms_init(environment()))
```

**Step 4: Verify syntax**

Run: `Rscript -e 'parse("r/pioneer/initializers.R"); cat("OK\n")'`
Expected: `OK`

**Step 5: Commit**

```bash
git add r/pioneer/initializers.R
git commit -m "Refactor pioneer initializer to use shared ms_init_values()"
```

---

### Task 4: Wire `initializers_ms.R` into the project sources

**Files:**
- Modify: `r/initializers.R` — add `source()` at top, or verify it's loaded via `init_project()`
- Modify: `.Rprofile` — `init_project()` function (if `initializers_ms.R` needs explicit sourcing)
- Modify: `targets/pioneer_targets.R` — verify `source()` chain includes new file

**Step 1: Check how initializers are currently loaded**

The targets pipeline sources `r/initializers.R` via `init_project()` in `.Rprofile`. The new `r/initializers_ms.R` must be sourced BEFORE both `r/initializers.R` and `r/pioneer/initializers.R` since they reference `ms_init_values`.

Add to `init_project()` in `.Rprofile`, before the existing `source(here("r", "initializers.R"))` line:
```r
source(here("r", "initializers_ms.R"))
```

**Step 2: Verify pioneer path**

Check that `targets/pioneer_targets.R` sources `r/pioneer/initializers.R` — if it goes through `init_project()`, the shared helper is already available. If it has its own `source()` calls, add `source(here("r", "initializers_ms.R"))` before the pioneer initializer source.

**Step 3: Verify syntax of all modified files**

Run:
```bash
Rscript -e 'parse("r/initializers_ms.R"); parse("r/initializers.R"); parse("r/pioneer/initializers.R"); cat("All OK\n")'
```
Expected: `All OK`

**Step 4: Commit**

```bash
git add .Rprofile
git commit -m "Source initializers_ms.R in init_project() before model-specific initializers"
```

---

### Task 5: Verify end-to-end (dry run)

**Step 1: Test that SCLC initializer produces valid output**

```r
Rscript -e '
  source(".Rprofile")
  init_project()
  # Load a cached stan_data from the store
  sd <- targets::tar_read(
    base_pioneer_stan_data_combined,
    store = "/mnt/data/analysis-results/karim_naguib/pioneer/main/_targets"
  )
  sd$fit_psa_data <- 1L
  sd$fit_multistate_data <- 1L
  init_fn <- create_pioneer_initializer(sd)
  init <- init_fn(1)
  # Check multistate params exist
  ms_params <- grep("^log_lambda_gp", names(init), value = TRUE)
  cat("Multistate params:", length(ms_params), "\n")
  stopifnot(length(ms_params) > 0)
  cat("SUCCESS: pioneer initializer produces multistate params\n")
'
```
Expected: `SUCCESS: pioneer initializer produces multistate params`

**Step 2: Commit final state**

```bash
git add -A
git commit -m "Verify modular initializer refactor end-to-end"
```

---

## Risk Notes

- The `environment()` trick used by `ms_init_values(env)` requires all multistate fields to be accessible in the calling `with(stan_data, {...})` scope. Both SCLC and pioneer stan_data include these fields when multistate is enabled, so this is safe.
- The `_fixed.R` file is left unchanged — it already has its own `ms_init_values_fixed()` with the same pattern. A future cleanup could unify the fixed/non-fixed helpers, but that's out of scope.
