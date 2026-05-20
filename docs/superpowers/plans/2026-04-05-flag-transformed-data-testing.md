# Flag Transformed Data Testing Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Extract flag-derived computations from `transformed_data.stan` files into tuple-returning Stan functions with a `strict` parameter, then test all flag permutations in a single batched Stan run per module.

**Architecture:** Each module's flag logic moves into named functions in `.stanfunctions` files. `strict=1` preserves production `fatal_error` behaviour; `strict=0` returns an `is_valid` sentinel instead, enabling invalid combos to be tested in the same batched run as valid ones. Test harnesses compile once, receive all flag combinations as padded data arrays, and return all derived quantities in a single `fixed_param` run.

**Tech Stack:** Stan 2.38 (tuple returns, dynamic local arrays), CmdStanR, testthat, posterior, tidyverse. R oracle functions `r_compute_n_enabled_groups`, `r_create_enabled_pos`, `r_get_global_group_idx` already exist in `tests/testthat/helper-pos.R`.

---

## File Map

### New files
| File | Responsibility |
|------|---------------|
| `stan/hierarchy.stanfunctions` | Generic `compute_level_module_flags()` for tr/frac/init |
| `stan/full_model.stanfunctions` | `compute_full_model_grid_flags()` |
| `tests/testthat/stan/test_hierarchy_flags_all.stan` | Batch harness: hierarchy module flags |
| `tests/testthat/stan/test_multistate_flags_all.stan` | Batch harness: all four multistate flag functions |
| `tests/testthat/stan/test_full_model_flags_all.stan` | Batch harness: full-model grid flags |
| `tests/testthat/test-stan-hierarchy-flags.R` | R tests for tr/frac/init flags |
| `tests/testthat/test-stan-multistate-flags.R` | R tests for multistate flags |
| `tests/testthat/test-stan-full-model-flags.R` | R tests for full-model flags |

### Modified files
| File | Change |
|------|--------|
| `stan/multistate.stanfunctions` | Add B1–B4 flag functions |
| `stan/modules/tr/transformed_data.stan` | Call `compute_level_module_flags` |
| `stan/modules/frac/transformed_data.stan` | Same |
| `stan/modules/init/transformed_data.stan` | Same |
| `stan/modules/multistate/transformed_data.stan` | Call B1–B4 functions |
| `stan/_full_model_transformed_data.stan` | Call `compute_full_model_grid_flags` |
| `stan/tumor/sf-ssm-log-space.stan` | Add `#include "hierarchy.stanfunctions"` and `#include "full_model.stanfunctions"` to `functions {}` |
| `stan/psa/pioneer.stan` | Same |

---

## Task 1: `compute_level_module_flags` — write failing test then implement

**Files:**
- Create: `stan/hierarchy.stanfunctions`
- Create: `tests/testthat/stan/test_hierarchy_flags_all.stan`
- Create: `tests/testthat/test-stan-hierarchy-flags.R`

- [ ] **Step 1: Create the Stan test harness**

Create `tests/testthat/stan/test_hierarchy_flags_all.stan`:

```stan
functions {
  #include "util.stanfunctions"
  #include "pos.stanfunctions"
  #include "hierarchy.stanfunctions"
}
data {
  int<lower=1> n_combos;
  int<lower=1> n_patients;
  int<lower=1> n_levels;
  int<lower=1> n_forecast_patients;
  array[n_levels] int n_forecast_groups_per_level;
  array[n_patients, n_levels] int patient_level_groups;

  // Per-combo flags (n_combos × n_levels)
  array[n_combos, n_levels] int enable_intercept;
  array[n_combos, n_levels] int enable_slope;
}
generated quantities {
  array[n_combos] int out_n_int;
  array[n_combos] int out_n_slp;
  array[n_combos, n_levels + 1] int out_pos_int;
  array[n_combos, n_levels + 1] int out_pos_slp;
  array[n_combos, n_patients, n_levels] int out_flat_int;
  array[n_combos, n_patients, n_levels] int out_flat_slp;

  for (c in 1:n_combos) {
    int n_ei;
    array[n_levels + 1] int pos_ei;
    array[n_patients, n_levels] int flat_ei;
    int n_es;
    array[n_levels + 1] int pos_es;
    array[n_patients, n_levels] int flat_es;

    (n_ei, pos_ei, flat_ei, n_es, pos_es, flat_es) = compute_level_module_flags(
      n_patients, n_levels, n_forecast_patients,
      n_forecast_groups_per_level, patient_level_groups,
      enable_intercept[c], enable_slope[c]
    );

    out_n_int[c]   = n_ei;
    out_n_slp[c]   = n_es;
    out_pos_int[c] = pos_ei;
    out_pos_slp[c] = pos_es;
    out_flat_int[c] = flat_ei;
    out_flat_slp[c] = flat_es;
  }
}
```

- [ ] **Step 2: Create R test file**

Create `tests/testthat/test-stan-hierarchy-flags.R`:

```r
library(testthat)
library(here)
library(posterior)
library(tidyverse)

# Structural data for n_levels=2: 2 trials × 3 patients = 6 total
make_hierarchy_struct <- function(n_levels) {
  if (n_levels == 2L) {
    list(
      n_patients                  = 6L,
      n_levels                    = 2L,
      n_forecast_patients         = 6L,
      n_forecast_groups_per_level = c(2L, 6L),
      patient_level_groups        = matrix(
        c(1L,1L, 1L,2L, 1L,3L, 2L,4L, 2L,5L, 2L,6L),
        nrow = 6L, ncol = 2L, byrow = TRUE
      )
    )
  } else {
    list(
      n_patients                  = 8L,
      n_levels                    = 3L,
      n_forecast_patients         = 8L,
      n_forecast_groups_per_level = c(2L, 4L, 8L),
      patient_level_groups        = matrix(
        c(1L,1L,1L, 1L,1L,2L, 1L,2L,3L, 1L,2L,4L,
          2L,3L,5L, 2L,3L,6L, 2L,4L,7L, 2L,4L,8L),
        nrow = 8L, ncol = 3L, byrow = TRUE
      )
    )
  }
}

# Build all 2^(2*n_levels) flag combos for a given n_levels
build_hierarchy_combos <- function(n_levels) {
  int_cols <- setNames(rep(list(0:1), n_levels), paste0("int_lv", seq_len(n_levels)))
  slp_cols <- setNames(rep(list(0:1), n_levels), paste0("slp_lv", seq_len(n_levels)))
  expand.grid(c(int_cols, slp_cols)) |> as.data.frame()
}

run_hierarchy_harness <- function(combos, struct) {
  nl <- struct$n_levels
  test_stan_function(
    here("tests", "testthat", "stan", "test_hierarchy_flags_all.stan"),
    c(
      struct,
      list(
        n_combos         = nrow(combos),
        enable_intercept = as.matrix(combos[, paste0("int_lv", seq_len(nl))]),
        enable_slope     = as.matrix(combos[, paste0("slp_lv", seq_len(nl))])
      )
    )
  )
}

for (nl in c(2L, 3L)) {
  struct <- make_hierarchy_struct(nl)
  combos <- build_hierarchy_combos(nl)
  nfg    <- struct$n_forecast_groups_per_level
  plg    <- struct$patient_level_groups
  np     <- struct$n_patients
  nfp    <- struct$n_forecast_patients

  # Compile + run once, reuse for all assertions
  fit <- run_hierarchy_harness(combos, struct)
  d   <- posterior::as_draws_df(fit$draws())

  test_that(sprintf("compute_level_module_flags: n_enabled matches oracle (n_levels=%d)", nl), {
    for (c in seq_len(nrow(combos))) {
      ei <- as.integer(combos[c, paste0("int_lv", seq_len(nl))])
      es <- as.integer(combos[c, paste0("slp_lv", seq_len(nl))])
      expect_equal(get_stan_val(d, "out_n_int", c),
                   r_compute_n_enabled_groups(nfg, ei),
                   label = sprintf("n_int combo %d", c))
      expect_equal(get_stan_val(d, "out_n_slp", c),
                   r_compute_n_enabled_groups(nfg, es),
                   label = sprintf("n_slp combo %d", c))
    }
  })

  test_that(sprintf("compute_level_module_flags: pos arrays match oracle (n_levels=%d)", nl), {
    for (c in seq_len(nrow(combos))) {
      ei       <- as.integer(combos[c, paste0("int_lv", seq_len(nl))])
      exp_pos  <- r_create_enabled_pos(nfg, ei)
      for (j in seq_along(exp_pos)) {
        expect_equal(get_stan_val(d, "out_pos_int", c, j), exp_pos[j],
                     label = sprintf("pos_int[%d,%d]", c, j))
      }
    }
  })

  test_that(sprintf("compute_level_module_flags: flat_idx==1 when level disabled (n_levels=%d)", nl), {
    for (c in seq_len(nrow(combos))) {
      ei <- as.integer(combos[c, paste0("int_lv", seq_len(nl))])
      for (lv in seq_len(nl)) {
        if (ei[lv] == 0L) {
          for (i in seq_len(np)) {
            expect_equal(get_stan_val(d, "out_flat_int", c, i, lv), 1L,
                         label = sprintf("flat_int[%d,%d,%d] disabled", c, i, lv))
          }
        }
      }
    }
  })

  test_that(sprintf("compute_level_module_flags: flat_idx correct when enabled (n_levels=%d)", nl), {
    for (c in seq_len(nrow(combos))) {
      ei      <- as.integer(combos[c, paste0("int_lv", seq_len(nl))])
      pos_int <- r_create_enabled_pos(nfg, ei)
      for (lv in seq_len(nl)) {
        if (ei[lv] != 0L) {
          for (i in seq_len(np)) {
            gid      <- plg[i, lv]
            expected <- if (lv == nl && gid > nfp) 1L
                        else r_get_global_group_idx(pos_int, lv, gid)
            expect_equal(get_stan_val(d, "out_flat_int", c, i, lv), expected,
                         label = sprintf("flat_int[%d,%d,%d]", c, i, lv))
          }
        }
      }
    }
  })
}
```

- [ ] **Step 3: Run test to confirm it fails (harness won't compile)**

```bash
Rscript -e 'testthat::test_file("tests/testthat/test-stan-hierarchy-flags.R")'
```
Expected: compilation error — `compute_level_module_flags` undefined.

- [ ] **Step 4: Implement `compute_level_module_flags`**

Create `stan/hierarchy.stanfunctions`:

```stan
/**
 * Generic hierarchy module derived quantities for tr/frac/init modules.
 * Computes enabled group counts, position arrays, and pre-computed flat
 * patient indices for both intercept and slope flag arrays.
 *
 * @param n_patients          Number of patients
 * @param n_levels            Number of hierarchy levels
 * @param n_forecast_patients Number of forecast (non-background) patients
 * @param n_forecast_groups_per_level Groups per level for parameter sizing
 * @param patient_level_groups Patient-to-group mapping [n_patients, n_levels]
 * @param enable_intercept    Per-level intercept enable flags [n_levels]
 * @param enable_slope        Per-level slope enable flags [n_levels]
 * @return Tuple (n_enabled_intercept, pos_intercept[n_levels+1],
 *               flat_intercept[n_patients,n_levels],
 *               n_enabled_slope, pos_slope[n_levels+1],
 *               flat_slope[n_patients,n_levels])
 */
tuple(int, array[] int, array[,] int, int, array[] int, array[,] int)
compute_level_module_flags(
  int n_patients,
  int n_levels,
  int n_forecast_patients,
  array[] int n_forecast_groups_per_level,
  array[,] int patient_level_groups,
  array[] int enable_intercept,
  array[] int enable_slope
) {
  int n_int = compute_n_enabled_groups(n_forecast_groups_per_level, enable_intercept);
  array[n_levels + 1] int pos_int = create_enabled_pos(n_forecast_groups_per_level, enable_intercept);
  array[n_patients, n_levels] int flat_int;
  for (i in 1:n_patients) {
    for (lv in 1:n_levels) {
      if (enable_intercept[lv]) {
        flat_int[i, lv] =
          (lv == n_levels && patient_level_groups[i, lv] > n_forecast_patients) ? 1
          : get_global_group_idx(pos_int, lv, patient_level_groups[i, lv]);
      } else {
        flat_int[i, lv] = 1;
      }
    }
  }

  int n_slp = compute_n_enabled_groups(n_forecast_groups_per_level, enable_slope);
  array[n_levels + 1] int pos_slp = create_enabled_pos(n_forecast_groups_per_level, enable_slope);
  array[n_patients, n_levels] int flat_slp;
  for (i in 1:n_patients) {
    for (lv in 1:n_levels) {
      if (enable_slope[lv]) {
        flat_slp[i, lv] =
          (lv == n_levels && patient_level_groups[i, lv] > n_forecast_patients) ? 1
          : get_global_group_idx(pos_slp, lv, patient_level_groups[i, lv]);
      } else {
        flat_slp[i, lv] = 1;
      }
    }
  }

  return (n_int, pos_int, flat_int, n_slp, pos_slp, flat_slp);
}
```

- [ ] **Step 5: Run test to confirm it passes**

```bash
Rscript -e 'testthat::test_file("tests/testthat/test-stan-hierarchy-flags.R")'
```
Expected: all tests PASS.

- [ ] **Step 6: Commit**

```bash
git add stan/hierarchy.stanfunctions \
        tests/testthat/stan/test_hierarchy_flags_all.stan \
        tests/testthat/test-stan-hierarchy-flags.R
git commit -m "test(hierarchy-flags): add compute_level_module_flags with batch flag tests"
```

---

## Task 2: Refactor tr/frac/init transformed_data to call the new function

**Files:**
- Modify: `stan/modules/tr/transformed_data.stan`
- Modify: `stan/modules/frac/transformed_data.stan`
- Modify: `stan/modules/init/transformed_data.stan`
- Modify: `stan/tumor/sf-ssm-log-space.stan`
- Modify: `stan/psa/pioneer.stan`

- [ ] **Step 1: Add `#include "hierarchy.stanfunctions"` to both model files**

In `stan/tumor/sf-ssm-log-space.stan`, add after `#include "pos.stanfunctions"`:
```stan
  #include "hierarchy.stanfunctions"
```

In `stan/psa/pioneer.stan`, add after `#include "pos.stanfunctions"`:
```stan
  #include "hierarchy.stanfunctions"
```

- [ ] **Step 2: Replace `stan/modules/tr/transformed_data.stan`**

```stan
// tr/transformed_data.stan
int n_enabled_groups_tr_intercept;
array[n_levels + 1] int enabled_level_pos_tr_intercept;
array[n_patients, n_levels] int patient_tr_intercept_flat_idx;
int n_enabled_groups_tr_slope;
array[n_levels + 1] int enabled_level_pos_tr_slope;
array[n_patients, n_levels] int patient_tr_slope_flat_idx;

(n_enabled_groups_tr_intercept, enabled_level_pos_tr_intercept, patient_tr_intercept_flat_idx,
 n_enabled_groups_tr_slope, enabled_level_pos_tr_slope, patient_tr_slope_flat_idx) =
  compute_level_module_flags(
    n_patients, n_levels, n_forecast_patients,
    n_forecast_groups_per_level, patient_level_groups,
    enable_level_intercept_tr, enable_level_cov_tr
  );
```

- [ ] **Step 3: Replace `stan/modules/frac/transformed_data.stan`**

```stan
// frac/transformed_data.stan
int n_enabled_groups_frac_intercept;
array[n_levels + 1] int enabled_level_pos_frac_intercept;
array[n_patients, n_levels] int patient_frac_intercept_flat_idx;
int n_enabled_groups_frac_slope;
array[n_levels + 1] int enabled_level_pos_frac_slope;
array[n_patients, n_levels] int patient_frac_slope_flat_idx;

(n_enabled_groups_frac_intercept, enabled_level_pos_frac_intercept, patient_frac_intercept_flat_idx,
 n_enabled_groups_frac_slope, enabled_level_pos_frac_slope, patient_frac_slope_flat_idx) =
  compute_level_module_flags(
    n_patients, n_levels, n_forecast_patients,
    n_forecast_groups_per_level, patient_level_groups,
    enable_level_intercept_frac, enable_level_cov_frac
  );
```

- [ ] **Step 4: Replace `stan/modules/init/transformed_data.stan`**

```stan
// init/transformed_data.stan
int n_enabled_groups_init_intercept;
array[n_levels + 1] int enabled_level_pos_init_intercept;
array[n_patients, n_levels] int patient_init_intercept_flat_idx;
int n_enabled_groups_init_slope;
array[n_levels + 1] int enabled_level_pos_init_slope;
array[n_patients, n_levels] int patient_init_slope_flat_idx;

(n_enabled_groups_init_intercept, enabled_level_pos_init_intercept, patient_init_intercept_flat_idx,
 n_enabled_groups_init_slope, enabled_level_pos_init_slope, patient_init_slope_flat_idx) =
  compute_level_module_flags(
    n_patients, n_levels, n_forecast_patients,
    n_forecast_groups_per_level, patient_level_groups,
    enable_level_intercept_init, enable_level_cov_init
  );
```

- [ ] **Step 5: Verify both models compile**

```bash
~/.cmdstan/cmdstan-2.38.0/bin/stanc \
  --include-paths=stan --include-paths=stan/tumor \
  stan/tumor/sf-ssm-log-space.stan

~/.cmdstan/cmdstan-2.38.0/bin/stanc \
  --include-paths=stan --include-paths=stan/psa \
  stan/psa/pioneer.stan
```
Expected: no errors.

- [ ] **Step 6: Run full test suite to confirm nothing regressed**

```bash
Rscript -e 'testthat::test_dir("tests/testthat")'
```
Expected: all tests pass.

- [ ] **Step 7: Commit**

```bash
git add stan/modules/tr/transformed_data.stan \
        stan/modules/frac/transformed_data.stan \
        stan/modules/init/transformed_data.stan \
        stan/tumor/sf-ssm-log-space.stan \
        stan/psa/pioneer.stan
git commit -m "refactor(hierarchy): tr/frac/init transformed_data use compute_level_module_flags"
```

---

## Task 3: Multistate B1–B4 flag functions — write failing tests then implement

This task adds all four multistate flag functions to `stan/multistate.stanfunctions` and a single batch test harness that covers all of them.

**Files:**
- Modify: `stan/multistate.stanfunctions`
- Create: `tests/testthat/stan/test_multistate_flags_all.stan`
- Create: `tests/testthat/test-stan-multistate-flags.R`

- [ ] **Step 1: Create the multistate batch test harness**

Create `tests/testthat/stan/test_multistate_flags_all.stan`:

```stan
functions {
  #include "util.stanfunctions"
  #include "pos.stanfunctions"
  #include "multistate.stanfunctions"
}
data {
  // === B1: time-scale routing ===
  int<lower=1> n_b1;
  array[n_b1] int b1_enable_ms_12;
  array[n_b1] int b1_ms_time_scale_12;

  // === B2: level baseline hazard ===
  int<lower=1> n_b2;
  int<lower=1> n_levels_b2;
  array[n_levels_b2] int n_forecast_groups_b2;
  array[n_b2, n_levels_b2] int b2_enable_level_baseline_hazard;

  // === B3: transition group counts ===
  int<lower=1> n_b3;
  array[n_b3] int b3_enable_01;
  array[n_b3] int b3_enable_02;
  array[n_b3] int b3_need_12_s;
  array[n_b3] int b3_need_12_t;
  array[n_b3] int b3_enable_03;
  array[n_b3] int b3_enable_32;
  int b3_n_enabled_baseline;
  int b3_n_gp_baseline;

  // === B4: visit-gated PSA covariate size ===
  int<lower=1> n_b4;
  array[n_b4] int b4_enable_visit_gated;
  array[n_b4] int b4_enable_latent;
  int b4_n_visits;
}
generated quantities {
  // B1 outputs
  array[n_b1] int b1_is_valid;
  array[n_b1] int b1_need_12_s_gp;
  array[n_b1] int b1_need_12_t_gp;
  array[n_b1] int b1_ms_12_t_has_intercept;

  // B2 outputs
  array[n_b2] int b2_is_valid;
  array[n_b2] int b2_any_re_level;
  array[n_b2, n_levels_b2] int b2_is_gp;
  array[n_b2] int b2_n_gp_groups;
  array[n_b2, n_levels_b2 + 1] int b2_gp_pos;
  array[n_b2, n_levels_b2 + 1] int b2_enabled_pos;

  // B3 outputs
  array[n_b3] int b3_is_valid;
  array[n_b3] int b3_n01;
  array[n_b3] int b3_n02;
  array[n_b3] int b3_n12s;
  array[n_b3] int b3_n12t;
  array[n_b3] int b3_n03;
  array[n_b3] int b3_n32;
  array[n_b3] int b3_gp01;
  array[n_b3] int b3_gp02;
  array[n_b3] int b3_gp12s;
  array[n_b3] int b3_gp12t;
  array[n_b3] int b3_gp03;
  array[n_b3] int b3_gp32;

  // B4 outputs
  array[n_b4] int b4_size;

  // B1
  for (c in 1:n_b1) {
    int is_v; int ns; int nt; int ti;
    (is_v, ns, nt, ti) = compute_ms_time_scale_flags(
      b1_enable_ms_12[c], b1_ms_time_scale_12[c], 0
    );
    b1_is_valid[c]              = is_v;
    b1_need_12_s_gp[c]          = ns;
    b1_need_12_t_gp[c]          = nt;
    b1_ms_12_t_has_intercept[c] = ti;
  }

  // B2
  for (c in 1:n_b2) {
    int is_v; int any_re;
    array[n_levels_b2] int is_gp;
    int n_gp;
    array[n_levels_b2 + 1] int gp_pos;
    array[n_levels_b2 + 1] int en_pos;
    (is_v, any_re, is_gp, n_gp, gp_pos, en_pos) = compute_ms_level_baseline_flags(
      n_levels_b2, n_forecast_groups_b2,
      b2_enable_level_baseline_hazard[c], 0
    );
    b2_is_valid[c]    = is_v;
    b2_any_re_level[c] = any_re;
    b2_is_gp[c]       = is_gp;
    b2_n_gp_groups[c] = n_gp;
    b2_gp_pos[c]      = gp_pos;
    b2_enabled_pos[c] = en_pos;
  }

  // B3
  for (c in 1:n_b3) {
    int is_v;
    int n01; int n02; int n12s; int n12t; int n03; int n32;
    int gp01; int gp02; int gp12s; int gp12t; int gp03; int gp32;
    (is_v, n01, n02, n12s, n12t, n03, n32,
     gp01, gp02, gp12s, gp12t, gp03, gp32) =
      compute_ms_transition_group_counts(
        b3_enable_01[c], b3_enable_02[c],
        b3_need_12_s[c], b3_need_12_t[c],
        b3_enable_03[c], b3_enable_32[c],
        b3_n_enabled_baseline, b3_n_gp_baseline, 0
      );
    b3_is_valid[c] = is_v;
    b3_n01[c]  = n01;  b3_n02[c]  = n02;
    b3_n12s[c] = n12s; b3_n12t[c] = n12t;
    b3_n03[c]  = n03;  b3_n32[c]  = n32;
    b3_gp01[c]  = gp01;  b3_gp02[c]  = gp02;
    b3_gp12s[c] = gp12s; b3_gp12t[c] = gp12t;
    b3_gp03[c]  = gp03;  b3_gp32[c]  = gp32;
  }

  // B4
  for (c in 1:n_b4) {
    b4_size[c] = compute_ms_obs_psa_covar_size(
      b4_enable_visit_gated[c], b4_enable_latent[c], b4_n_visits
    );
  }
}
```

- [ ] **Step 2: Create R test file**

Create `tests/testthat/test-stan-multistate-flags.R`:

```r
library(testthat)
library(here)
library(posterior)
library(tidyverse)

# Oracle functions
r_ms_time_scale_flags <- function(enable_12, time_scale) {
  if (time_scale < 0L || time_scale > 2L)
    return(list(is_valid=0L, need_s=0L, need_t=0L, t_int=0L))
  need_s <- as.integer(enable_12 && time_scale %in% c(1L, 2L))
  need_t <- as.integer(enable_12 && time_scale %in% c(0L, 2L))
  list(is_valid=1L, need_s=need_s, need_t=need_t,
       t_int=as.integer(need_t && !need_s))
}

r_ms_level_baseline_flags <- function(n_forecast_groups, enable_hazard) {
  if (any(enable_hazard < 0L) || any(enable_hazard > 3L))
    return(list(is_valid=0L))
  is_gp   <- as.integer(enable_hazard == 3L)
  any_re  <- as.integer(max(enable_hazard) >= 2L)
  n_gp    <- r_compute_n_enabled_groups(n_forecast_groups, is_gp)
  gp_pos  <- r_create_enabled_pos(n_forecast_groups, is_gp)
  en_pos  <- r_create_enabled_pos(n_forecast_groups, enable_hazard)
  list(is_valid=1L, any_re=any_re, is_gp=is_gp, n_gp=n_gp,
       gp_pos=gp_pos, en_pos=en_pos)
}

r_ms_transition_counts <- function(en01, en02, need_s, need_t, en03, en32,
                                    n_enabled, n_gp) {
  if (en32 && !en03) return(list(is_valid=0L))
  f <- function(gate) list(n=if(gate) n_enabled else 0L, gp=if(gate) n_gp else 0L)
  r01 <- f(en01);  r02 <- f(en02); r12s <- f(need_s); r12t <- f(need_t)
  r03 <- f(en03);  r32 <- f(en32)
  list(is_valid=1L,
       n01=r01$n,   n02=r02$n,   n12s=r12s$n, n12t=r12t$n, n03=r03$n, n32=r32$n,
       gp01=r01$gp, gp02=r02$gp, gp12s=r12s$gp, gp12t=r12t$gp, gp03=r03$gp, gp32=r32$gp)
}

# Build all combos and compile harness once
b1_combos <- expand.grid(enable_12=0:1, time_scale=c(-1L, 0L:2L, 3L))

n_levels_b2    <- 2L
n_groups_b2    <- c(2L, 4L)  # 2 trials, 4 arms
b2_combos      <- expand.grid(lv1=0:3, lv2=c(-1L, 0:3, 4L)) |>
  as.data.frame()

b3_combos <- expand.grid(
  en01=0:1, en02=0:1, need_s=0:1, need_t=0:1,
  en03=0:1, en32=0:1
) |> as.data.frame()

b4_combos   <- expand.grid(visit_gated=0:1, latent=0:1)
b4_n_visits <- 50L

stan_data <- list(
  n_b1                          = nrow(b1_combos),
  b1_enable_ms_12               = b1_combos$enable_12,
  b1_ms_time_scale_12           = b1_combos$time_scale,

  n_b2                          = nrow(b2_combos),
  n_levels_b2                   = n_levels_b2,
  n_forecast_groups_b2          = n_groups_b2,
  b2_enable_level_baseline_hazard = as.matrix(b2_combos),

  n_b3                          = nrow(b3_combos),
  b3_enable_01                  = b3_combos$en01,
  b3_enable_02                  = b3_combos$en02,
  b3_need_12_s                  = b3_combos$need_s,
  b3_need_12_t                  = b3_combos$need_t,
  b3_enable_03                  = b3_combos$en03,
  b3_enable_32                  = b3_combos$en32,
  b3_n_enabled_baseline         = 6L,  # 2 groups × 3 (arbitrary, tests gating)
  b3_n_gp_baseline              = 2L,  # subset of above that are GP mode

  n_b4                          = nrow(b4_combos),
  b4_enable_visit_gated         = b4_combos$visit_gated,
  b4_enable_latent              = b4_combos$latent,
  b4_n_visits                   = b4_n_visits
)

fit <- test_stan_function(
  here("tests", "testthat", "stan", "test_multistate_flags_all.stan"),
  stan_data
)
d <- posterior::as_draws_df(fit$draws())

# --- B1 tests ---
test_that("compute_ms_time_scale_flags: all combos correct", {
  for (c in seq_len(nrow(b1_combos))) {
    exp <- r_ms_time_scale_flags(b1_combos$enable_12[c], b1_combos$time_scale[c])
    expect_equal(get_stan_val(d, "b1_is_valid", c),              exp$is_valid, label=sprintf("b1_is_valid[%d]", c))
    if (exp$is_valid) {
      expect_equal(get_stan_val(d, "b1_need_12_s_gp", c),          exp$need_s,   label=sprintf("need_s[%d]", c))
      expect_equal(get_stan_val(d, "b1_need_12_t_gp", c),          exp$need_t,   label=sprintf("need_t[%d]", c))
      expect_equal(get_stan_val(d, "b1_ms_12_t_has_intercept", c), exp$t_int,    label=sprintf("t_int[%d]", c))
    }
  }
})

test_that("compute_ms_time_scale_flags: invalid time_scale returns is_valid=0", {
  invalid <- which(b1_combos$time_scale < 0L | b1_combos$time_scale > 2L)
  for (c in invalid)
    expect_equal(get_stan_val(d, "b1_is_valid", c), 0L,
                 label = sprintf("b1_is_valid invalid combo %d", c))
})

# --- B2 tests ---
test_that("compute_ms_level_baseline_flags: all combos correct", {
  for (c in seq_len(nrow(b2_combos))) {
    hz  <- as.integer(b2_combos[c, ])
    exp <- r_ms_level_baseline_flags(n_groups_b2, hz)
    expect_equal(get_stan_val(d, "b2_is_valid", c), exp$is_valid,
                 label = sprintf("b2_is_valid[%d]", c))
    if (exp$is_valid) {
      expect_equal(get_stan_val(d, "b2_any_re_level", c), exp$any_re,
                   label = sprintf("any_re[%d]", c))
      expect_equal(get_stan_val(d, "b2_n_gp_groups", c), exp$n_gp,
                   label = sprintf("n_gp[%d]", c))
      for (lv in seq_len(n_levels_b2)) {
        expect_equal(get_stan_val(d, "b2_is_gp", c, lv), exp$is_gp[lv],
                     label = sprintf("is_gp[%d,%d]", c, lv))
      }
      for (j in seq_along(exp$en_pos)) {
        expect_equal(get_stan_val(d, "b2_enabled_pos", c, j), exp$en_pos[j],
                     label = sprintf("en_pos[%d,%d]", c, j))
      }
    }
  }
})

# --- B3 tests ---
test_that("compute_ms_transition_group_counts: valid combos correct", {
  for (c in seq_len(nrow(b3_combos))) {
    r <- b3_combos[c, ]
    exp <- r_ms_transition_counts(r$en01, r$en02, r$need_s, r$need_t,
                                   r$en03, r$en32, 6L, 2L)
    expect_equal(get_stan_val(d, "b3_is_valid", c), exp$is_valid,
                 label = sprintf("b3_is_valid[%d]", c))
    if (exp$is_valid) {
      for (nm in c("n01","n02","n12s","n12t","n03","n32",
                   "gp01","gp02","gp12s","gp12t","gp03","gp32")) {
        expect_equal(get_stan_val(d, paste0("b3_", nm), c), exp[[nm]],
                     label = sprintf("%s[%d]", nm, c))
      }
    }
  }
})

test_that("compute_ms_transition_group_counts: enable_ms_32=1 && enable_ms_03=0 is invalid", {
  invalid <- which(b3_combos$en32 == 1L & b3_combos$en03 == 0L)
  for (c in invalid)
    expect_equal(get_stan_val(d, "b3_is_valid", c), 0L,
                 label = sprintf("b3_is_valid invalid[%d]", c))
})

# --- B4 tests ---
test_that("compute_ms_obs_psa_covar_size: all 4 combos correct", {
  for (c in seq_len(nrow(b4_combos))) {
    expected_size <- if (b4_combos$visit_gated[c] && !b4_combos$latent[c])
                       b4_n_visits else 0L
    expect_equal(get_stan_val(d, "b4_size", c), expected_size,
                 label = sprintf("b4_size[%d]", c))
  }
})
```

- [ ] **Step 3: Run test to confirm it fails**

```bash
Rscript -e 'testthat::test_file("tests/testthat/test-stan-multistate-flags.R")'
```
Expected: compilation error — the four functions are not yet defined.

- [ ] **Step 4: Add the four flag functions to `stan/multistate.stanfunctions`**

Append to the end of `stan/multistate.stanfunctions`:

```stan
// ============================================================================
// Flag-derived quantity functions — extracted from transformed_data.stan
// strict=1: call fatal_error on invalid input (production)
// strict=0: return is_valid=0 sentinel (test harness)
// ============================================================================

/**
 * Compute time-scale routing flags for the 1→2 transition.
 * @return (is_valid, need_12_s_gp, need_12_t_gp, ms_12_t_has_intercept)
 */
tuple(int, int, int, int) compute_ms_time_scale_flags(
  int enable_ms_12,
  int ms_time_scale_12,
  int strict
) {
  if (ms_time_scale_12 < 0 || ms_time_scale_12 > 2) {
    if (strict) fatal_error("ms_time_scale_12 must be 0-2; got ", ms_time_scale_12);
    return (0, 0, 0, 0);
  }
  int need_s   = enable_ms_12 && (ms_time_scale_12 == 1 || ms_time_scale_12 == 2);
  int need_t   = enable_ms_12 && (ms_time_scale_12 == 0 || ms_time_scale_12 == 2);
  int t_int    = need_t && !need_s;
  return (1, need_s, need_t, t_int);
}

/**
 * Compute level baseline hazard flags.
 * @return (is_valid, any_re_level, ms_level_baseline_is_gp[n_levels],
 *          n_gp_groups_ms_baseline,
 *          gp_level_pos_ms_baseline[n_levels+1],
 *          enabled_level_pos_ms_baseline[n_levels+1])
 */
tuple(int, int, array[] int, int, array[] int, array[] int)
compute_ms_level_baseline_flags(
  int n_levels,
  array[] int n_forecast_groups_per_level,
  array[] int enable_ms_level_baseline_hazard,
  int strict
) {
  for (lv in 1:n_levels) {
    if (enable_ms_level_baseline_hazard[lv] < 0 ||
        enable_ms_level_baseline_hazard[lv] > 3) {
      if (strict)
        fatal_error("enable_ms_level_baseline_hazard[", lv, "] must be 0-3");
      return (0, 0, rep_array(0, n_levels), 0,
              rep_array(0, n_levels + 1), rep_array(0, n_levels + 1));
    }
  }
  array[n_levels] int is_gp;
  for (lv in 1:n_levels)
    is_gp[lv] = (enable_ms_level_baseline_hazard[lv] == 3) ? 1 : 0;

  int any_re     = max(to_array_1d(enable_ms_level_baseline_hazard)) >= 2 ? 1 : 0;
  int n_gp       = compute_n_enabled_groups(n_forecast_groups_per_level, is_gp);
  array[n_levels + 1] int gp_pos =
    create_enabled_pos(n_forecast_groups_per_level, is_gp);
  array[n_levels + 1] int en_pos =
    create_enabled_pos(n_forecast_groups_per_level, enable_ms_level_baseline_hazard);

  return (1, any_re, is_gp, n_gp, gp_pos, en_pos);
}

/**
 * Compute per-transition enabled/GP group counts from baseline hazard summary.
 * @return (is_valid, n_01, n_02, n_12s, n_12t, n_03, n_32,
 *          gp_01, gp_02, gp_12s, gp_12t, gp_03, gp_32)
 */
tuple(int, int, int, int, int, int, int, int, int, int, int, int, int)
compute_ms_transition_group_counts(
  int enable_ms_01, int enable_ms_02,
  int need_12_s_gp, int need_12_t_gp,
  int enable_ms_03, int enable_ms_32,
  int n_enabled_groups_ms_baseline,
  int n_gp_groups_ms_baseline,
  int strict
) {
  if (enable_ms_32 && !enable_ms_03) {
    if (strict) fatal_error("enable_ms_32=1 requires enable_ms_03=1");
    return (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
  }
  int N = n_enabled_groups_ms_baseline;
  int G = n_gp_groups_ms_baseline;
  return (1,
    enable_ms_01  ? N : 0, enable_ms_02  ? N : 0,
    need_12_s_gp  ? N : 0, need_12_t_gp  ? N : 0,
    enable_ms_03  ? N : 0, enable_ms_32  ? N : 0,
    enable_ms_01  ? G : 0, enable_ms_02  ? G : 0,
    need_12_s_gp  ? G : 0, need_12_t_gp  ? G : 0,
    enable_ms_03  ? G : 0, enable_ms_32  ? G : 0);
}

/**
 * Compute size of flat observed-PSA covariate vector for visit-gated 0→1.
 * Non-zero only when visit-gated mode is on and latent PSA is NOT used.
 */
int compute_ms_obs_psa_covar_size(
  int enable_ms_visit_gated_01,
  int enable_ms_visit_gated_latent_01,
  int n_visits
) {
  return (enable_ms_visit_gated_01 && !enable_ms_visit_gated_latent_01)
    ? n_visits : 0;
}
```

- [ ] **Step 5: Run test to confirm it passes**

```bash
Rscript -e 'testthat::test_file("tests/testthat/test-stan-multistate-flags.R")'
```
Expected: all tests PASS.

- [ ] **Step 6: Commit**

```bash
git add stan/multistate.stanfunctions \
        tests/testthat/stan/test_multistate_flags_all.stan \
        tests/testthat/test-stan-multistate-flags.R
git commit -m "test(multistate-flags): add B1-B4 flag functions with batch permutation tests"
```

---

## Task 4: Refactor `multistate/transformed_data.stan` to call the new functions

**Files:**
- Modify: `stan/modules/multistate/transformed_data.stan`

- [ ] **Step 1: Replace the inline flag derivations with function calls**

In `stan/modules/multistate/transformed_data.stan`, replace lines 20–100 (from the `fatal_error` guard through the position array declarations) with:

```stan
// --- Validate transitions and derive time-scale routing ---
if (!enable_ms_03 && enable_ms_32) {
  fatal_error("enable_ms_32=1 requires enable_ms_03=1");
}

int is_valid_time_scale; int need_12_s_gp; int need_12_t_gp; int ms_12_t_has_intercept;
(is_valid_time_scale, need_12_s_gp, need_12_t_gp, ms_12_t_has_intercept) =
  compute_ms_time_scale_flags(enable_ms_12, ms_time_scale_12, 1);

// --- Level baseline hazard derived quantities ---
int any_re_level;
array[n_levels] int ms_level_baseline_is_gp;
int n_gp_groups_ms_baseline;
array[n_levels + 1] int gp_level_pos_ms_baseline;
array[n_levels + 1] int enabled_level_pos_ms_baseline;
(, any_re_level, ms_level_baseline_is_gp, n_gp_groups_ms_baseline,
 gp_level_pos_ms_baseline, enabled_level_pos_ms_baseline) =
  compute_ms_level_baseline_flags(
    n_levels, n_forecast_groups_per_level, enable_ms_level_baseline_hazard, 1
  );

// --- Per-transition group counts ---
int n_enabled_groups_ms_baseline = compute_n_enabled_groups(
  n_forecast_groups_per_level, enable_ms_level_baseline_hazard
);
int n_enabled_groups_ms_baseline_01; int n_enabled_groups_ms_baseline_02;
int n_enabled_groups_ms_baseline_12_s; int n_enabled_groups_ms_baseline_12_t;
int n_enabled_groups_ms_baseline_03; int n_enabled_groups_ms_baseline_32;
int n_gp_groups_ms_baseline_01; int n_gp_groups_ms_baseline_02;
int n_gp_groups_ms_baseline_12_s; int n_gp_groups_ms_baseline_12_t;
int n_gp_groups_ms_baseline_03; int n_gp_groups_ms_baseline_32;
(,
 n_enabled_groups_ms_baseline_01, n_enabled_groups_ms_baseline_02,
 n_enabled_groups_ms_baseline_12_s, n_enabled_groups_ms_baseline_12_t,
 n_enabled_groups_ms_baseline_03, n_enabled_groups_ms_baseline_32,
 n_gp_groups_ms_baseline_01, n_gp_groups_ms_baseline_02,
 n_gp_groups_ms_baseline_12_s, n_gp_groups_ms_baseline_12_t,
 n_gp_groups_ms_baseline_03, n_gp_groups_ms_baseline_32) =
  compute_ms_transition_group_counts(
    enable_ms_01, enable_ms_02,
    need_12_s_gp, need_12_t_gp,
    enable_ms_03, enable_ms_32,
    n_enabled_groups_ms_baseline, n_gp_groups_ms_baseline, 1
  );

// --- Covariate slope counts ---
int n_enabled_groups_ms_slope = compute_n_enabled_groups(
  n_forecast_groups_per_level, enable_ms_level_cov
);
array[n_levels + 1] int enabled_level_pos_ms_slope = create_enabled_pos(
  n_forecast_groups_per_level, enable_ms_level_cov
);
```

Also replace the `ms_obs_psa_covar_flat` declaration (line ~199):
```stan
// Before:
vector[enable_ms_visit_gated_01 && !enable_ms_visit_gated_latent_01 ? size(t_patient_visits) : 0] ms_obs_psa_covar_flat;

// After:
vector[compute_ms_obs_psa_covar_size(
  enable_ms_visit_gated_01, enable_ms_visit_gated_latent_01, size(t_patient_visits)
)] ms_obs_psa_covar_flat;
```

Note: Stan's tuple destructuring uses `,` to skip fields: `(, b, c) = f()` skips the first return value (`is_valid` sentinel — always 1 in strict mode, so can be discarded in production).

- [ ] **Step 2: Verify both models compile**

```bash
~/.cmdstan/cmdstan-2.38.0/bin/stanc \
  --include-paths=stan --include-paths=stan/tumor \
  stan/tumor/sf-ssm-log-space.stan

~/.cmdstan/cmdstan-2.38.0/bin/stanc \
  --include-paths=stan --include-paths=stan/psa \
  stan/psa/pioneer.stan

~/.cmdstan/cmdstan-2.38.0/bin/stanc \
  --include-paths=stan \
  stan/ms-standalone.stan
```
Expected: no errors from any model.

- [ ] **Step 3: Run full test suite**

```bash
Rscript -e 'testthat::test_dir("tests/testthat")'
```
Expected: all tests pass.

- [ ] **Step 4: Commit**

```bash
git add stan/modules/multistate/transformed_data.stan
git commit -m "refactor(multistate): transformed_data calls B1-B4 flag functions"
```

---

## Task 5: Full-model grid flags — write failing test then implement

**Files:**
- Create: `stan/full_model.stanfunctions`
- Create: `tests/testthat/stan/test_full_model_flags_all.stan`
- Create: `tests/testthat/test-stan-full-model-flags.R`
- Modify: `stan/tumor/sf-ssm-log-space.stan`
- Modify: `stan/psa/pioneer.stan`
- Modify: `stan/_full_model_transformed_data.stan`

- [ ] **Step 1: Create the Stan test harness**

Create `tests/testthat/stan/test_full_model_flags_all.stan`:

```stan
functions {
  #include "util.stanfunctions"
  #include "full_model.stanfunctions"
}
data {
  int<lower=1> n_combos;
  array[n_combos] int enable_pop_pn;        // enable_pop_process_noise_tr
  array[n_combos] int enable_patient_pn;    // enable_patient_process_noise_tr
  array[n_combos] int enable_states_grid;   // enable_states_full_grid
  array[n_combos] int enable_ms_tv_cov;     // enable_ms_pop_time_varying_cov
  array[n_combos] int n_tv_covar;           // n_time_varying_covar
  array[n_combos] int enable_ms_01;
  array[n_combos] int enable_visit_gated;   // enable_ms_visit_gated_01
  array[n_combos] int enable_02_tv_cov;     // enable_ms_02_time_varying_cov
}
generated quantities {
  array[n_combos] int out_any_process_noise;
  array[n_combos] int out_need_states_full_grid;

  for (c in 1:n_combos) {
    int apn; int nsg;
    (apn, nsg) = compute_full_model_grid_flags(
      enable_pop_pn[c], enable_patient_pn[c],
      enable_states_grid[c],
      enable_ms_tv_cov[c], n_tv_covar[c],
      enable_ms_01[c], enable_visit_gated[c],
      enable_02_tv_cov[c]
    );
    out_any_process_noise[c]      = apn;
    out_need_states_full_grid[c]  = nsg;
  }
}
```

- [ ] **Step 2: Create R test file**

Create `tests/testthat/test-stan-full-model-flags.R`:

```r
library(testthat)
library(here)
library(posterior)
library(tidyverse)

r_full_model_grid_flags <- function(pop_pn, pat_pn, states_grid,
                                     ms_tv_cov, n_tv, ms_01, visit_gated,
                                     tv_02_cov) {
  any_pn    <- as.integer(pop_pn || pat_pn)
  need_grid <- as.integer(
    any_pn || states_grid ||
    (ms_tv_cov && n_tv > 0 &&
       (ms_01 && !visit_gated || tv_02_cov))
  )
  list(any_pn = any_pn, need_grid = need_grid)
}

combos <- expand.grid(
  pop_pn      = 0:1,
  pat_pn      = 0:1,
  states_grid = 0:1,
  ms_tv_cov   = 0:1,
  n_tv        = c(0L, 1L, 3L),
  ms_01       = 0:1,
  visit_gated = 0:1,
  tv_02_cov   = 0:1
) |> as.data.frame()

fit <- test_stan_function(
  here("tests", "testthat", "stan", "test_full_model_flags_all.stan"),
  list(
    n_combos          = nrow(combos),
    enable_pop_pn     = combos$pop_pn,
    enable_patient_pn = combos$pat_pn,
    enable_states_grid= combos$states_grid,
    enable_ms_tv_cov  = combos$ms_tv_cov,
    n_tv_covar        = combos$n_tv,
    enable_ms_01      = combos$ms_01,
    enable_visit_gated= combos$visit_gated,
    enable_02_tv_cov  = combos$tv_02_cov
  )
)
d <- posterior::as_draws_df(fit$draws())

test_that("compute_full_model_grid_flags: all combos correct", {
  for (c in seq_len(nrow(combos))) {
    r   <- combos[c, ]
    exp <- r_full_model_grid_flags(r$pop_pn, r$pat_pn, r$states_grid,
                                    r$ms_tv_cov, r$n_tv, r$ms_01,
                                    r$visit_gated, r$tv_02_cov)
    expect_equal(get_stan_val(d, "out_any_process_noise",     c), exp$any_pn,
                 label = sprintf("any_pn[%d]", c))
    expect_equal(get_stan_val(d, "out_need_states_full_grid", c), exp$need_grid,
                 label = sprintf("need_grid[%d]", c))
  }
})

test_that("compute_full_model_grid_flags: need_grid=1 when any process noise", {
  pn_combos <- which(combos$pop_pn == 1L | combos$pat_pn == 1L)
  for (c in pn_combos)
    expect_equal(get_stan_val(d, "out_need_states_full_grid", c), 1L,
                 label = sprintf("need_grid pn combo %d", c))
})

test_that("compute_full_model_grid_flags: need_grid=0 when all gates off", {
  off_combos <- which(
    combos$pop_pn == 0L & combos$pat_pn == 0L &
    combos$states_grid == 0L &
    !(combos$ms_tv_cov == 1L & combos$n_tv > 0L &
      (combos$ms_01 == 1L & combos$visit_gated == 0L | combos$tv_02_cov == 1L))
  )
  for (c in off_combos)
    expect_equal(get_stan_val(d, "out_need_states_full_grid", c), 0L,
                 label = sprintf("need_grid off combo %d", c))
})
```

- [ ] **Step 3: Run test to confirm it fails**

```bash
Rscript -e 'testthat::test_file("tests/testthat/test-stan-full-model-flags.R")'
```
Expected: compilation error — `compute_full_model_grid_flags` undefined.

- [ ] **Step 4: Create `stan/full_model.stanfunctions`**

```stan
/**
 * Compute process-noise and full-grid routing flags for the full biomarker model.
 *
 * @return (enable_any_process_noise_tr, need_states_full_grid)
 */
tuple(int, int) compute_full_model_grid_flags(
  int enable_pop_process_noise_tr,
  int enable_patient_process_noise_tr,
  int enable_states_full_grid,
  int enable_ms_pop_time_varying_cov,
  int n_time_varying_covar,
  int enable_ms_01,
  int enable_ms_visit_gated_01,
  int enable_ms_02_time_varying_cov
) {
  int any_pn = enable_pop_process_noise_tr || enable_patient_process_noise_tr;
  int need_grid = any_pn || enable_states_full_grid ||
    (enable_ms_pop_time_varying_cov && n_time_varying_covar > 0 &&
     (enable_ms_01 && !enable_ms_visit_gated_01 || enable_ms_02_time_varying_cov));
  return (any_pn, need_grid);
}
```

- [ ] **Step 5: Add `#include "full_model.stanfunctions"` to both model files**

In `stan/tumor/sf-ssm-log-space.stan` `functions {}` block, add after `#include "hierarchy.stanfunctions"`:
```stan
  #include "full_model.stanfunctions"
```

In `stan/psa/pioneer.stan` `functions {}` block, add after `#include "hierarchy.stanfunctions"`:
```stan
  #include "full_model.stanfunctions"
```

- [ ] **Step 6: Run tests to confirm they pass**

```bash
Rscript -e 'testthat::test_file("tests/testthat/test-stan-full-model-flags.R")'
```
Expected: all tests PASS.

- [ ] **Step 7: Refactor `stan/_full_model_transformed_data.stan`**

Replace lines 44–57 (the inline `enable_any_process_noise_tr` and `need_states_full_grid` derivations):

```stan
// Before:
// int enable_any_process_noise_tr = enable_pop_process_noise_tr || enable_patient_process_noise_tr;
// int need_states_full_grid = enable_any_process_noise_tr || enable_states_full_grid || ...

// After:
int enable_any_process_noise_tr;
int need_states_full_grid;
(enable_any_process_noise_tr, need_states_full_grid) = compute_full_model_grid_flags(
  enable_pop_process_noise_tr, enable_patient_process_noise_tr,
  enable_states_full_grid,
  enable_ms_pop_time_varying_cov, n_time_varying_covar,
  enable_ms_01, enable_ms_visit_gated_01,
  enable_ms_02_time_varying_cov
);
print("need_states_full_grid = ", need_states_full_grid,
      " (process_noise=", enable_any_process_noise_tr,
      ", full_grid_flag=", enable_states_full_grid,
      ", ungated_01=", enable_ms_01 && !enable_ms_visit_gated_01, ")");
```

- [ ] **Step 8: Verify all models compile**

```bash
~/.cmdstan/cmdstan-2.38.0/bin/stanc \
  --include-paths=stan --include-paths=stan/tumor \
  stan/tumor/sf-ssm-log-space.stan

~/.cmdstan/cmdstan-2.38.0/bin/stanc \
  --include-paths=stan --include-paths=stan/psa \
  stan/psa/pioneer.stan
```
Expected: no errors.

- [ ] **Step 9: Run full test suite**

```bash
Rscript -e 'testthat::test_dir("tests/testthat")'
```
Expected: all tests pass.

- [ ] **Step 10: Commit**

```bash
git add stan/full_model.stanfunctions \
        stan/_full_model_transformed_data.stan \
        stan/tumor/sf-ssm-log-space.stan \
        stan/psa/pioneer.stan \
        tests/testthat/stan/test_full_model_flags_all.stan \
        tests/testthat/test-stan-full-model-flags.R
git commit -m "test(full-model-flags): add compute_full_model_grid_flags with batch permutation tests"
```
