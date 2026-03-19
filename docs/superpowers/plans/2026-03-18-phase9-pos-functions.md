# Plan: Phase 9 — pos.stanfunctions Batch Coverage

**Date:** 2026-03-18
**Branch:** karim/tests
**Source file:** `stan/pos.stanfunctions`
**Target assertions:** ~80 new assertions

## Background

`pos.stanfunctions` implements the position-array pattern — Stan's idiomatic replacement for
ragged arrays. Position arrays are used *everywhere* in the model (patient visits, tumors,
trials). These functions are the lowest-level building blocks; testing them thoroughly provides
the foundation for all higher-level tests.

The key invariant: for groups with sizes `[s_1, s_2, ..., s_n]`, `create_pos` returns
`[1, 1+s_1, 1+s_1+s_2, ..., 1+sum(s_i)]` (1-based). `get_pos(pos, i)` returns
`(pos[i], pos[i+1]-1)`.

## Functions to cover

| Function | Overloads |
|----------|-----------|
| `create_pos` | 3 (from sizes, from sub_pos, from range) |
| `get_pos` | 3 (single group, range, double-pos) |
| `get_pos_size` | 2 (single group, all groups) |
| `get_pos_total_size` | 1 |
| `get_int_sub_array` | 2 (single group, range) |
| `get_real_sub_array` | 2 |
| `get_sub_vert_matrix` | 1 |
| `get_min_pos` / `get_max_pos` | 2 each |
| `get_int` | 1 |
| `validate_pos` | 1 |
| `create_enabled_pos` | 1 |
| `compute_n_enabled_groups` | 1 |
| `get_global_group_idx` | 1 |

---

## Task 1 — R helpers `helper-pos.R`

**File:** `tests/testthat/helper-pos.R`

```r
# create_pos: cumulative sum, 1-based
r_create_pos <- function(n_x) {
  c(1L, 1L + cumsum(as.integer(n_x)))
}

# get_pos: (start, end) for group i (1-based)
r_get_pos <- function(pos, i) c(pos[i], pos[i + 1L] - 1L)

# get_pos_size
r_get_pos_size_single <- function(pos, i) pos[i + 1L] - pos[i]
r_get_pos_size_all    <- function(pos) diff(as.integer(pos))

# get_int_sub_array
r_get_int_sub_array <- function(full, pos, n) {
  s <- r_get_pos(pos, n)
  full[s[1]:s[2]]
}

# create_enabled_pos
r_create_enabled_pos <- function(n_groups, enabled) {
  n <- length(n_groups)
  pos <- integer(n + 1L)
  pos[1] <- 1L
  for (i in seq_len(n)) pos[i + 1L] <- pos[i] + if (enabled[i]) n_groups[i] else 0L
  pos
}

# compute_n_enabled_groups
r_compute_n_enabled_groups <- function(n_groups, enabled) {
  sum(n_groups[enabled == 1L])
}

# get_global_group_idx
r_get_global_group_idx <- function(level_pos, level, group_id) {
  level_pos[level] + group_id - 1L
}

# validate_pos: sorted check
r_validate_pos <- function(pos) {
  stopifnot(all(diff(pos) >= 0L))
  pos
}
```

**Test file:** `tests/testthat/test-helper-pos.R`
**Assertions (~20):**
- `r_create_pos(c(2, 0, 3))` → `c(1, 3, 3, 6)` (zero-size group handled correctly)
- `r_get_pos(c(1,3,3,6), 1)` → `c(1, 2)`; group 2 → `c(3, 2)` (empty: start > end)
- `r_get_pos_size_single(pos, 2)` → 0 for empty group
- `r_get_int_sub_array(1:5, c(1,3,5,6), 3)` → `5`
- `r_create_enabled_pos(c(3,5), c(1,0))` → `c(1,4,4)` (level 2 disabled)
- `r_compute_n_enabled_groups(c(3,5), c(1,0))` → 3
- `r_get_global_group_idx(c(1,4,9), 2, 3)` → 6

---

## Task 2 — Stan harness `test_pos_all.stan`

**File:** `tests/testthat/stan/test_pos_all.stan`

```stan
functions {
  #include util.stanfunctions
  #include pos.stanfunctions
}
data {
  int<lower=1> N_GROUPS;
  array[N_GROUPS] int group_sizes;
  array[N_GROUPS + 1] int pos_input;  // pre-built position array for get_pos tests
  int<lower=1> N_DATA;
  array[N_DATA] int full_int;
  array[N_DATA] real full_real;
  // Enabled pos
  array[N_GROUPS] int enabled;
  // level mapping
  int<lower=1> N_LEVELS;
  array[N_LEVELS + 1] int level_pos;
  int query_level;
  int query_group;
  // validate test (valid input)
  array[4] int valid_pos;
  // get_int test
  int get_int_group;
  int get_int_elem;
}
generated quantities {
  // create_pos
  array[N_GROUPS + 1] int created_pos = create_pos(group_sizes);
  // get_pos (single group)
  int gp_start; int gp_end;
  (gp_start, gp_end) = get_pos(pos_input, 1);
  // get_pos_size
  int gps_single = get_pos_size(pos_input, 1);
  array[N_GROUPS] int gps_all = get_pos_size(pos_input);
  // get_pos_total_size
  int total_size = get_pos_total_size(pos_input);
  // get_int_sub_array (int, group 2)
  array[get_pos_size(pos_input, 2)] int sub_int = get_int_sub_array(full_int, pos_input, 2);
  // get_real_sub_array
  array[get_pos_size(pos_input, 2)] real sub_real = get_real_sub_array(full_real, pos_input, 2);
  // get_min_pos / get_max_pos
  array[N_GROUPS] int min_vals = get_min_pos(full_int, pos_input);
  array[N_GROUPS] int max_vals = get_max_pos(full_int, pos_input);
  // create_enabled_pos
  array[N_GROUPS + 1] int enabled_pos = create_enabled_pos(group_sizes, enabled);
  int n_enabled = compute_n_enabled_groups(group_sizes, enabled);
  // get_global_group_idx
  int global_idx = get_global_group_idx(level_pos, query_level, query_group);
  // get_int (element within group)
  int elem = get_int(full_int, pos_input, get_int_group, get_int_elem);
  // validate_pos: returns pos unchanged for valid input
  array[4] int validated = validate_pos(valid_pos);
}
```

**Test file:** `tests/testthat/test-stan-pos.R`

```r
library(testthat)

test_that("create_pos, get_pos, get_pos_size match R oracles", {
  group_sizes <- c(2L, 0L, 3L, 1L)
  pos         <- r_create_pos(group_sizes)          # c(1,3,3,6,7)
  full_int    <- c(10L, 20L, 30L, 40L, 50L, 60L)   # 6 elements total
  full_real   <- as.numeric(full_int)

  stan_data <- list(
    N_GROUPS = length(group_sizes), group_sizes = group_sizes,
    pos_input = pos,
    N_DATA = length(full_int), full_int = full_int, full_real = full_real,
    enabled = c(1L, 0L, 1L, 1L),
    N_LEVELS = 2L, level_pos = c(1L, 4L, 9L),
    query_level = 2L, query_group = 3L,
    valid_pos = c(1L, 3L, 3L, 6L),
    get_int_group = 3L, get_int_elem = 2L
  )
  fit <- test_stan_function("test_pos_all", stan_data)
  d   <- posterior::as_draws_df(fit$draws())

  # create_pos
  for (i in seq_along(pos)) {
    expect_equal(get_stan_val(d, "created_pos", i), pos[i],
      label = paste0("created_pos[", i, "]"))
  }

  # get_pos for group 1: start=1, end=2
  expect_equal(get_stan_val(d, "gp_start"), 1L)
  expect_equal(get_stan_val(d, "gp_end"),   2L)

  # get_pos_size
  expect_equal(get_stan_val(d, "gps_single"), 2L)  # group 1 has size 2
  expect_equal(get_stan_val(d, "gps_all", 2),  0L)  # group 2 is empty
  expect_equal(get_stan_val(d, "gps_all", 3),  3L)  # group 3 has size 3

  # get_pos_total_size
  expect_equal(get_stan_val(d, "total_size"), 6L)

  # get_int_sub_array: group 2 is empty, so skip; group 3 = elements 3,4,5
  # (test with group 2 = first element of group 3 region)
  # sub_int for group 2 is empty (size 0) — Stan won't output it; skip
  # Verify sub_int via group min/max instead:
  expect_equal(get_stan_val(d, "min_vals", 1), 10L)  # group 1: min(10,20)=10
  expect_equal(get_stan_val(d, "max_vals", 1), 20L)  # group 1: max(10,20)=20
  expect_equal(get_stan_val(d, "min_vals", 3), 30L)  # group 3: min(30,40,50)=30
  expect_equal(get_stan_val(d, "max_vals", 3), 50L)  # group 3: max(30,40,50)=50

  # create_enabled_pos: groups 2 disabled
  expected_enabled <- r_create_enabled_pos(group_sizes, c(1L, 0L, 1L, 1L))
  for (i in seq_along(expected_enabled)) {
    expect_equal(get_stan_val(d, "enabled_pos", i), expected_enabled[i],
      label = paste0("enabled_pos[", i, "]"))
  }

  # compute_n_enabled_groups: groups 1+3+4 = 2+3+1 = 6
  expect_equal(get_stan_val(d, "n_enabled"), 6L)

  # get_global_group_idx: level=2 (pos[2]=4), group=3 → 4+3-1=6
  expect_equal(get_stan_val(d, "global_idx"), 6L)

  # get_int: group 3, element 2 → full_int[pos[3] + 2 - 1] = full_int[3+2-1] = full_int[4] = 40
  expect_equal(get_stan_val(d, "elem"), 40L)

  # validate_pos: returns input unchanged
  for (i in 1:4) {
    expect_equal(get_stan_val(d, "validated", i), c(1L, 3L, 3L, 6L)[i])
  }
})

test_that("create_pos with all-zero sizes: pos = c(1,1,...,1)", {
  group_sizes <- c(0L, 0L, 0L)
  pos_expected <- c(1L, 1L, 1L, 1L)
  stan_data <- list(
    N_GROUPS = 3L, group_sizes = group_sizes,
    pos_input = pos_expected,
    N_DATA = 1L, full_int = 999L, full_real = 999.0,
    enabled = c(1L, 1L, 1L),
    N_LEVELS = 2L, level_pos = c(1L, 2L, 4L),
    query_level = 1L, query_group = 1L,
    valid_pos = c(1L, 1L, 1L, 1L),
    get_int_group = 1L, get_int_elem = 1L
  )
  fit <- test_stan_function("test_pos_all", stan_data)
  d   <- posterior::as_draws_df(fit$draws())
  for (i in 1:4) {
    expect_equal(get_stan_val(d, "created_pos", i), pos_expected[i])
  }
  expect_equal(get_stan_val(d, "total_size"), 0L)
})
```

**Assertion count:** ~35 Stan assertions

---

## Summary

| Task | File | Assertions |
|------|------|-----------|
| 1 | `helper-pos.R` + `test-helper-pos.R` | ~20 R |
| 2 | `test_pos_all.stan` + `test-stan-pos.R` | ~35 Stan |
| **Total** | | **~55 new assertions** |
