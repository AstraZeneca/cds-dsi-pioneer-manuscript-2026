# Flag Transformed Data Testing Design

**Date**: 2026-04-05  
**Status**: Pending implementation

## Problem

The model has many feature flags across modules (tr, frac, init, multistate,
full-model). These flags gate parameter sizing, GP routing, transition
activation, and hierarchy structure through chains of nested logical
expressions in each module's `transformed_data.stan`. Bugs here are silent:
wrong things get enabled or disabled without any obvious error.

## Goal

Test that all flag-derived computations in `transformed_data.stan` produce
correct results across all valid flag permutations, and that invalid
combinations correctly trigger errors. Tests must use the actual production
Stan logic, not parallel R re-implementations.

---

## Approach

### Step 1 — Extract flag logic into tuple-returning Stan functions

Move all flag-derived computations out of `transformed_data.stan` inline code
into named tuple-returning functions. Each function takes a `strict` parameter:
- `strict=1` (production): calls `fatal_error()` on invalid inputs
- `strict=0` (test harness): returns an `is_valid=0` sentinel instead

`transformed_data.stan` becomes a thin caller that destructures tuples.
Test harnesses call the same functions in a loop over all combos.

### Step 2 — Batch test harnesses

One Stan test harness per module family, taking all flag combinations as
padded data arrays. A single `fixed_param` run with 1 iteration processes all
combos. Results are output arrays indexed by combo, checked in R.

Invalid combos are included in the same batch — the `is_valid` output
distinguishes them from valid ones, no separate runs needed.

### Step 3 — R test files with full expand.grid enumeration

R generates all flag permutations via `expand.grid()`, filters to valid/invalid
sets using simple predicates, runs the harness once, then asserts invariants.
For the hierarchy modules, the existing R oracle functions
(`r_compute_n_enabled_groups`, `r_create_enabled_pos`, `r_get_global_group_idx`
in `helper-pos.R`) are reused directly — no new oracle code needed.

---

## Modules in Scope

### A. Hierarchy modules: tr, frac, init

All three are structurally identical — same pattern, different flag names. A
single generic function handles all three.

**New function** in `stan/hierarchy.stanfunctions`:

```stan
tuple(
  int,              // n_enabled_groups_intercept
  array[] int,      // enabled_level_pos_intercept [n_levels+1]
  array[,] int,     // patient_intercept_flat_idx [n_patients, n_levels]
  int,              // n_enabled_groups_slope
  array[] int,      // enabled_level_pos_slope [n_levels+1]
  array[,] int      // patient_slope_flat_idx [n_patients, n_levels]
) compute_level_module_flags(
  int n_patients,
  int n_levels,
  int n_forecast_patients,
  array[] int n_forecast_groups_per_level,
  array[,] int patient_level_groups,
  array[] int enable_intercept,
  array[] int enable_slope
)
```

No `strict` parameter needed — all flag combos for hierarchy modules are valid.

**Modified** `stan/modules/tr/transformed_data.stan`,
`stan/modules/frac/transformed_data.stan`,
`stan/modules/init/transformed_data.stan`: each calls
`compute_level_module_flags` with its own flag names and destructures the tuple.

**New** `tests/testthat/stan/test_hierarchy_flags_all.stan`: takes padded
`array[n_combos, max_levels] int enable_intercept_combos` etc., loops over
combos, calls `compute_level_module_flags`, outputs result arrays.

**New** `tests/testthat/test-stan-hierarchy-flags.R`:
- `expand.grid` over `enable_intercept[1..n_levels]` and `enable_slope[1..n_levels]` with n_levels ∈ {2, 3}
- Oracle: `r_compute_n_enabled_groups`, `r_create_enabled_pos` (already in `helper-pos.R`)
- Assertions: group counts match oracle, `patient_flat_idx[i,lv] == 1` when level disabled, correct global index when enabled

---

### B. Multistate module

#### B1. Time-scale routing

**New function** in `stan/multistate.stanfunctions`:

```stan
tuple(int, int, int, int) compute_ms_time_scale_flags(
  int enable_ms_12,
  int ms_time_scale_12,
  int strict
)
// returns: (is_valid, need_12_s_gp, need_12_t_gp, ms_12_t_has_intercept)
```

Invalid input: `ms_time_scale_12 < 0 || ms_time_scale_12 > 2`.

**Invariants tested**:
- `need_12_s_gp == 0` when `enable_ms_12 == 0`
- `need_12_t_gp == 0` when `enable_ms_12 == 0`
- `need_12_s_gp == 1` iff `enable_ms_12 == 1 && ms_time_scale_12 ∈ {1,2}`
- `need_12_t_gp == 1` iff `enable_ms_12 == 1 && ms_time_scale_12 ∈ {0,2}`
- `ms_12_t_has_intercept == 1` iff pure Markov (`need_12_t_gp && !need_12_s_gp`)
- When extended mode (scale=2): both GPs on, `ms_12_t_has_intercept == 0`

#### B2. Level baseline hazard flags

**New function** in `stan/multistate.stanfunctions`:

```stan
tuple(int, array[] int, int, array[] int, array[] int) compute_ms_level_baseline_flags(
  int n_levels,
  array[] int n_forecast_groups_per_level,
  array[] int enable_ms_level_baseline_hazard,
  int strict
)
// returns: (is_valid, ms_level_baseline_is_gp[n_levels],
//           n_gp_groups_ms_baseline,
//           gp_level_pos_ms_baseline[n_levels+1],
//           enabled_level_pos_ms_baseline[n_levels+1])
```

Invalid input: any `enable_ms_level_baseline_hazard[lv] < 0 || > 3`.

**Invariants tested**:
- `ms_level_baseline_is_gp[lv] == 1` iff `enable_ms_level_baseline_hazard[lv] == 3`
- `any_re_level == 1` iff `max(enable_ms_level_baseline_hazard) >= 2`
- `n_gp_groups == 0` when all levels have `enable_ms_level_baseline_hazard ∈ {0,1,2}`
- Position arrays consistent with enabled group counts (via `r_create_enabled_pos`)

#### B3. Transition group counts

**New function** in `stan/multistate.stanfunctions`:

```stan
tuple(int, int, int, int, int, int) compute_ms_transition_group_counts(
  int enable_ms_01,
  int enable_ms_02,
  int need_12_s_gp,
  int need_12_t_gp,
  int enable_ms_03,
  int enable_ms_32,
  int n_enabled_groups_ms_baseline,
  int n_gp_groups_ms_baseline,
  int strict
)
// returns: (is_valid,
//           n_enabled_groups_ms_baseline_01,
//           n_enabled_groups_ms_baseline_02,
//           n_enabled_groups_ms_baseline_12_s,
//           n_enabled_groups_ms_baseline_12_t,
//           n_enabled_groups_ms_baseline_03,
//           n_enabled_groups_ms_baseline_32)
// (and mirrored n_gp_groups_* variants)
```

Invalid input: `enable_ms_32 == 1 && enable_ms_03 == 0`.

**Invariants tested**:
- `n_enabled_groups_*_XX == 0` when the corresponding `enable_ms_XX == 0`
- All enabled transitions share the same `n_enabled_groups` value (same level config)
- `n_gp_groups <= n_enabled_groups` for every transition
- `is_valid == 0` for the `enable_ms_32=1, enable_ms_03=0` combination

#### B4. Visit-gated PSA covariate size

**New function** in `stan/multistate.stanfunctions`:

```stan
int compute_ms_obs_psa_covar_size(
  int enable_ms_visit_gated_01,
  int enable_ms_visit_gated_latent_01,
  int n_visits
)
// returns: n_visits when visit_gated=1 && latent=0, else 0
```

**Invariants tested** (4 combinations of the two binary flags):
- Size is `n_visits` only when `visit_gated=1 && latent=0`
- Size is 0 in all other cases

---

### C. Full-model cross-module flags

**New function** in a new `stan/full_model.stanfunctions`:

```stan
tuple(int, int) compute_full_model_grid_flags(
  int enable_pop_process_noise_tr,
  int enable_patient_process_noise_tr,
  int enable_states_full_grid,
  int enable_ms_pop_time_varying_cov,
  int n_time_varying_covar,
  int enable_ms_01,
  int enable_ms_visit_gated_01,
  int enable_ms_02_time_varying_cov
)
// returns: (enable_any_process_noise_tr, need_states_full_grid)
```

**Invariants tested**:
- `enable_any_process_noise_tr == 1` iff either process noise flag is on
- `need_states_full_grid == 1` when process noise enabled
- `need_states_full_grid == 1` when `enable_states_full_grid == 1`
- `need_states_full_grid == 1` when ungated time-varying 0→1 covariate is active
  (i.e., `enable_ms_pop_time_varying_cov && n_time_varying_covar > 0 && enable_ms_01 && !enable_ms_visit_gated_01`)
- `need_states_full_grid == 0` when none of the above

---

## File Summary

### New files

| File | Purpose |
|------|---------|
| `stan/hierarchy.stanfunctions` | Generic `compute_level_module_flags()` |
| `stan/full_model.stanfunctions` | `compute_full_model_grid_flags()` |
| `tests/testthat/stan/test_hierarchy_flags_all.stan` | Batch harness: tr/frac/init |
| `tests/testthat/stan/test_multistate_flags_all.stan` | Batch harness: multistate |
| `tests/testthat/stan/test_full_model_flags_all.stan` | Batch harness: full-model |
| `tests/testthat/test-stan-hierarchy-flags.R` | R tests for hierarchy modules |
| `tests/testthat/test-stan-multistate-flags.R` | R tests for multistate |
| `tests/testthat/test-stan-full-model-flags.R` | R tests for full-model flags |

### Modified files

| File | Change |
|------|--------|
| `stan/modules/tr/transformed_data.stan` | Call `compute_level_module_flags()` |
| `stan/modules/frac/transformed_data.stan` | Same |
| `stan/modules/init/transformed_data.stan` | Same |
| `stan/modules/multistate/transformed_data.stan` | Call B1–B4 functions |
| `stan/_full_model_transformed_data.stan` | Call `compute_full_model_grid_flags()` |
| `stan/multistate.stanfunctions` | Add B1–B4 functions (already `#include`d everywhere) |

Note: `stan/hierarchy.stanfunctions` and `stan/full_model.stanfunctions` must be
added to `#include` lists in all model files that currently include
`modules/tr/transformed_data.stan` etc.

---

## What is NOT in scope

- Data validation guards (`ms_final_state > ms_max_state`) — these depend on
  patient-level data, not flags; tested separately if needed
- `ms_ic_gap_01` computation — depends on patient data, not flags
- GP knot grid computation — deterministic from `max_all_t` and `ms_gp_grid_step`, not gated by flags
- `patient_ms_baseline_flat_idx` / `patient_ms_slope_flat_idx` correctness — the
  underlying `get_global_group_idx` function is already tested in `test-stan-pos.R`

---

## Conventions

- All new `.stanfunctions` files follow the existing style in `stan/multistate.stanfunctions`
- Tuple destructuring uses Stan 2.31+ syntax: `(a, b, c) = f(...)`
- `strict=1` in all production `transformed_data.stan` calls
- `strict=0` in all test harness calls
- Padded data matrices follow the pattern established in `test-stan-pos.R`
  (`n_combos`, `max_*` dimensions, zero padding)
