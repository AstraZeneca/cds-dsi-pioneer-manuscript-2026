# Plan: Phase 8 — util.stanfunctions Batch Coverage

**Date:** 2026-03-18
**Branch:** karim/tests
**Source file:** `stan/util.stanfunctions`
**Target assertions:** ~90 new assertions

## Background

`util.stanfunctions` contains ~25 distinct functions (many with overloads) covering integer array
utilities, date conversion, uniqueness counting, level-mapping, missing-value detection, and
patient compaction. Most are pure and easily oracle-tested. The highest-value targets are the
ones used most widely in the pipeline: `which`, `id2idx`, `calendar_date_to_study_date`,
`num_unique`, `unique`, `find_first`, and `get_level2level_idx`.

## Functions to cover (grouped by harness)

### Harness A — Simple array utilities
| Function | Signature |
|----------|-----------|
| `count_positive` | `array[] int → int` |
| `which` / `which(mask, inverse)` | `array[] int → array[] int` |
| `rep_each` | `(array[] int, int) → array[] int` |
| `months_to_weeks` | `int → real` |
| `calendar_date_to_study_date` | 3 overloads |
| `study_date_to_calendar_date` | 2 overloads |
| `num_unique` | 2 primary overloads (single array + by-pos) |
| `unique` | 2 overloads |
| `id2idx` | 2 overloads |

### Harness B — Level mapping + compaction
| Function | Signature |
|----------|-----------|
| `get_level2level_idx` | 3 overloads |
| `find_first` | 2 overloads |
| `create_compact_patient_mapping` | `(mask, n_obs) → (patients, idx_map)` |
| `create_compact_group_pos` | `(compact_patients, patient_group, n_groups) → pos` |
| `diag_matrix` | `(real, int) → matrix` |
| `standardize_tumor_sizes` | `vector → (mean, sd, standardized)` |

---

## Task 1 — R helpers `helper-util.R`

**File:** `tests/testthat/helper-util.R`

```r
r_count_positive <- function(arr) sum(arr > 0L)

r_which <- function(mask, inverse = FALSE) {
  if (inverse) which(mask <= 0L) else which(mask > 0L)
}

r_rep_each <- function(to_repeat, repeats) {
  as.integer(rep(to_repeat, each = repeats))
}

r_months_to_weeks <- function(mon) mon * 365.25 / (7 * 12)

r_calendar_date_to_study_date_scalar <- function(first, cal) cal - first + 1L

r_calendar_date_to_study_date_vec <- function(first_vec, cal_vec) {
  as.integer(cal_vec - first_vec + 1L)
}

r_num_unique <- function(x, post_treatment = TRUE) {
  if (post_treatment) x <- x[x > 0L]
  length(unique(x))
}

r_unique_sorted <- function(x, post_treatment = TRUE) {
  if (post_treatment) x <- x[x > 0L]
  sort(unique(x))
}

r_id2idx <- function(id) as.integer(id - min(id) + 1L)

r_get_level2level_idx <- function(hi_level, low_level) {
  # For each element in low_level, find its position in hi_level (sorted match)
  sorted_low <- sort(low_level)
  idx <- integer(length(low_level))
  for (k in seq_along(sorted_low)) {
    pos <- match(sorted_low[k], hi_level)
    idx[k] <- pos
  }
  idx
}

r_find_first_arr <- function(all, what, n_succ) {
  n <- length(all)
  for (i in seq_len(n - n_succ + 1L)) {
    if (all(all[i:(i + n_succ - 1L)] %in% what)) return(i)
  }
  0L
}

r_create_compact_patient_mapping <- function(observed_mask) {
  observed_patients <- which(observed_mask == 1L)
  n_patients <- length(observed_mask)
  patient_to_compact <- integer(n_patients)
  patient_to_compact[observed_patients] <- seq_along(observed_patients)
  list(observed_patients = observed_patients, patient_to_compact_idx = patient_to_compact)
}
```

**Test file:** `tests/testthat/test-helper-util.R`
**Assertions (~20):**
- `r_count_positive`: mixed array; all-negative → 0; all-positive → n
- `r_which`: basic mask; inverse mask; empty result
- `r_rep_each`: `c(1,2)` × 3 = `c(1,1,1,2,2,2)`
- `r_months_to_weeks`: 12 months ≈ 52.18 weeks (365.25/7)
- `r_calendar_date_to_study_date_scalar`: day 100, first_day=90 → study_day=11
- `r_num_unique`: with duplicates; post_treatment flag removes negatives
- `r_id2idx`: `c(3,1,2)` → `c(3,1,2)`; offset case
- `r_find_first_arr`: found vs not found; n_succ=2

---

## Task 2 — Stan harness A `test_util_array_all.stan`

**File:** `tests/testthat/stan/test_util_array_all.stan`

```stan
functions {
  #include util.stanfunctions
  #include pos.stanfunctions
}
data {
  int<lower=1> N;
  array[N] int arr_int;
  array[N] int mask;
  array[N] int to_repeat;
  int repeats;
  int mon;
  int first_calendar;
  int calendar_date;
  array[N] int id_arr;
  // find_first
  int<lower=1> N_ALL;
  array[N_ALL] int search_arr;
  int search_what;
  int n_succ;
}
generated quantities {
  int cnt_pos = count_positive(arr_int);
  array[count_positive(mask)] int which_idx = which(mask);
  array[count_positive(mask)] int which_inv_idx = which(mask, 0);
  array[N * repeats] int rep_each_out = rep_each(to_repeat, repeats);
  real weeks = months_to_weeks(mon);
  int study_date = calendar_date_to_study_date(first_calendar, calendar_date);
  array[N] int study_dates = calendar_date_to_study_date(arr_int, calendar_date);
  int num_uniq = num_unique(arr_int);
  array[num_unique(arr_int)] int uniq_vals = unique(arr_int);
  array[N] int idx_arr = id2idx(id_arr);
  int ff = find_first(search_arr, search_what);
}
```

**Test file:** `tests/testthat/test-stan-util-array.R`

```r
library(testthat)

test_that("Stan util: count_positive, which, rep_each, months_to_weeks", {
  arr  <- c(3L, -1L, 0L, 5L, 2L)
  mask <- c(1L, 0L, 1L, 0L, 1L)

  stan_data <- list(
    N = 5L, arr_int = arr, mask = mask,
    to_repeat = c(10L, 20L), repeats = 2L,
    mon = 12L, first_calendar = 100L, calendar_date = 110L,
    id_arr = c(3L, 1L, 4L, 1L, 5L),
    N_ALL = 6L, search_arr = c(1L, 3L, 2L, 3L, 3L, 1L),
    search_what = 3L, n_succ = 2L
  )
  fit <- test_stan_function("test_util_array_all", stan_data)
  d   <- posterior::as_draws_df(fit$draws())

  # count_positive
  expect_equal(get_stan_val(d, "cnt_pos"), r_count_positive(arr))
  # which
  expect_equal(get_stan_val(d, "which_idx", 1), 1L)
  expect_equal(get_stan_val(d, "which_idx", 2), 3L)
  expect_equal(get_stan_val(d, "which_idx", 3), 5L)
  # rep_each
  expect_equal(get_stan_val(d, "rep_each_out", 1), 10L)
  expect_equal(get_stan_val(d, "rep_each_out", 2), 10L)
  expect_equal(get_stan_val(d, "rep_each_out", 3), 20L)
  expect_equal(get_stan_val(d, "rep_each_out", 4), 20L)
  # months_to_weeks
  expect_equal(get_stan_val(d, "weeks"), r_months_to_weeks(12L), tolerance = 1e-6)
  # calendar_date_to_study_date
  expect_equal(get_stan_val(d, "study_date"), 11L)
  # num_unique (arr has values 3,-1,0,5,2; all unique → 5)
  expect_equal(get_stan_val(d, "num_uniq"), 5L)
  # id2idx: c(3,1,4,1,5) → min=1 → c(3,1,4,1,5)
  expect_equal(get_stan_val(d, "idx_arr", 1), 3L)
  expect_equal(get_stan_val(d, "idx_arr", 2), 1L)
  # find_first: search for 3 with n_succ=2 → first consecutive pair is at index 4
  expect_equal(get_stan_val(d, "ff"), 4L)
})
```

**Assertion count (Harness A):** ~22 Stan assertions

---

## Task 3 — Stan harness B `test_util_level_map_all.stan`

**File:** `tests/testthat/stan/test_util_level_map_all.stan`

```stan
functions {
  #include util.stanfunctions
  #include pos.stanfunctions
}
data {
  int N_HI; int N_LOW;
  array[N_HI] int hi_level;
  array[N_LOW] int low_level;
  // compact patient mapping
  int N_PAT;
  int N_OBS;
  array[N_PAT] int observed_mask;
  // standardize_tumor_sizes
  int<lower=1> N_TUM;
  vector[N_TUM] tumor_sizes;
  // diag_matrix
  real diag_val;
  int diag_n;
}
generated quantities {
  array[N_LOW] int l2l_idx = get_level2level_idx(hi_level, low_level);
  array[N_OBS] int compact_pats;
  array[N_PAT] int compact_idx;
  (compact_pats, compact_idx) = create_compact_patient_mapping(observed_mask, N_OBS);
  real tum_mean; real tum_sd; vector[N_TUM] tum_std;
  (tum_mean, tum_sd, tum_std) = standardize_tumor_sizes(tumor_sizes);
  matrix[diag_n, diag_n] diag_out = diag_matrix(diag_val, diag_n);
}
```

**Test file assertions (~20):**
- `get_level2level_idx`: each `low_level` value maps to its position in `hi_level`
- `create_compact_patient_mapping`: compact_pats = original IDs where mask=1; compact_idx
  nonzero only for observed patients; compact_idx[compact_pats[k]] = k
- `standardize_tumor_sizes`: mean ≈ 0 for standardized; sd = 1
- `diag_matrix(v, n)`: diagonal = v; off-diagonal = 0

---

## Summary

| Task | File | Assertions |
|------|------|-----------|
| 1 | `helper-util.R` + `test-helper-util.R` | ~20 R |
| 2 | `test_util_array_all.stan` + `test-stan-util-array.R` | ~22 Stan |
| 3 | `test_util_level_map_all.stan` + `test-stan-util-level-map.R` | ~20 Stan |
| **Total** | | **~62 new assertions** |
