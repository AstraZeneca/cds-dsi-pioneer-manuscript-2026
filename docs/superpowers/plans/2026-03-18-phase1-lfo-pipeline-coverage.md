# Phase 1 LFO Pipeline Coverage — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add tests for 4 untested functions that form the LFO pipeline backbone: `get_oos_patients_idx`, `truncate_at_max_time`, `find_first_week`, and `find_first_forecast_week`.

**Architecture:** Each task follows TDD: write the R test first (it fails), then create/fix the Stan harness (it passes). R reference implementations (`helper-pfs.R`) serve as oracles for expected values, mirroring Stan logic exactly. Reuses the `test_stan_function` + `get_stan_val` pattern established in `helper-stan.R`.

**Tech Stack:** R + testthat, Stan (cmdstanr), posterior, existing helpers (`helper-stan.R`, `helper-lfo.R`).

---

## File Structure

| File | Action | Purpose |
|---|---|---|
| `tests/testthat/stan/test_get_oos_patients_idx_all.stan` | **Create** | Stan harness for `get_oos_patients_idx` |
| `tests/testthat/test-stan-get_oos_patients_idx.R` | **Create** | R test: 6 cases, asserts `oos_idx_out` |
| `tests/testthat/stan/test_pfs_functions_all.stan` | **Modify** | Remove broken inline section (does not compile) |
| `tests/testthat/helper-pfs.R` | **Create** | R references: `r_truncate_at_max_time`, `r_find_first`, `r_find_first_week`, `r_find_first_forecast_week` |
| `tests/testthat/test-stan-pfs-truncate_at_max_time.R` | **Create** | R test: 6 cases for `truncate_at_max_time` |
| `tests/testthat/test-helper-pfs.R` | **Create** | R-only unit tests for `helper-pfs.R` references |
| `tests/testthat/stan/test_pfs_find_first_week_all.stan` | **Create** | Stan harness for `find_first_week` + `find_first_forecast_week` |
| `tests/testthat/test-stan-pfs-find_first_week.R` | **Create** | R test: 7 `find_first_week` cases + 3 `find_first_forecast_week` cases |

---

## Background: The Four Functions

All four are in `stan/` (one in `lfo.stanfunctions`, three in `pfs.stanfunctions`):

```
lfo.stanfunctions:
  get_oos_patients_idx(sorted_last_visit_calendar_day, cutoff_calendar_day)
    → array[] int   (index of first OOS patient for each cutoff; 0 = none)

pfs.stanfunctions:
  truncate_at_max_time(pfs, right_censored, max_time)
    → tuple(array[] int, array[] int)   (truncated pfs, updated right_censored)

  find_first_week(arr, values, min_run_length, curr_visits, forecast_time, max_all_t)
    → tuple(int, int)   (week, right_censored)

  find_first_forecast_week(arr, values, min_run_length, forecast_time, max_all_t)
    → tuple(int, int)   delegates to find_first_week with empty curr_visits
```

**Key semantic notes:**
- `get_oos_patients_idx`: walks sorted patients; for each cutoff finds first patient whose last visit calendar day is **strictly greater** than the cutoff. Patients with last_visit == cutoff are NOT OOS.
- `find_first_week`: calls `find_first(arr, values, min_run_length)` → 1-based index → maps to week via `map_idx_to_week`. Index 0 means not found → week = max_all_t, right_censored = 1.
- `map_idx_to_week`: if `idx <= n_obs`, returns `curr_visits[idx]`; else `forecast_time[idx - n_obs]`. Forecast indexing starts at 1 (not overlapping the last obs visit).
- `truncate_at_max_time`: `pfs[i] > max_time[i]` → set `pfs[i] = max_time[i]`, `right_censored[i] = 1`. Boundary `pfs[i] == max_time[i]` is NOT truncated.

---

## Task 1: `get_oos_patients_idx` — Stan harness + R test

**Files:**
- Create: `tests/testthat/stan/test_get_oos_patients_idx_all.stan`
- Create: `tests/testthat/test-stan-get_oos_patients_idx.R`

---

- [ ] **Step 1: Write the failing R test**

Create `tests/testthat/test-stan-get_oos_patients_idx.R`:

```r
library(testthat)
library(cmdstanr)
library(posterior)
source(here::here("tests/testthat/helper-stan.R"))

test_that("get_oos_patients_idx: all cases produce correct OOS start indices", {
  cases <- list(
    list(
      case_name        = "single_patient_after_cutoff",
      n_patients       = 1L,
      n_cutoffs        = 1L,
      sorted_last_visit = c(20L),
      cutoffs          = c(15L),
      expected         = c(1L)   # patient 1: last=20 > 15 → OOS from idx 1
    ),
    list(
      case_name        = "single_patient_before_cutoff",
      n_patients       = 1L,
      n_cutoffs        = 1L,
      sorted_last_visit = c(10L),
      cutoffs          = c(15L),
      expected         = c(0L)   # patient 1: last=10 ≤ 15 → no OOS patients
    ),
    list(
      case_name        = "two_patients_split",
      n_patients       = 2L,
      n_cutoffs        = 1L,
      sorted_last_visit = c(10L, 20L),
      cutoffs          = c(15L),
      expected         = c(2L)   # patient 1 skipped (10 ≤ 15), patient 2 is first OOS
    ),
    list(
      case_name        = "two_patients_both_oos",
      n_patients       = 2L,
      n_cutoffs        = 1L,
      sorted_last_visit = c(20L, 30L),
      cutoffs          = c(15L),
      expected         = c(1L)   # patient 1: last=20 > 15 → OOS from idx 1
    ),
    list(
      case_name        = "exact_boundary_not_oos",
      n_patients       = 3L,
      n_cutoffs        = 1L,
      sorted_last_visit = c(10L, 15L, 20L),
      cutoffs          = c(15L),
      expected         = c(3L)   # patients 1,2 have last ≤ 15; patient 3 (last=20) is first OOS
    ),
    list(
      case_name        = "two_cutoffs_three_patients",
      n_patients       = 3L,
      n_cutoffs        = 2L,
      sorted_last_visit = c(10L, 20L, 30L),
      cutoffs          = c(15L, 25L),
      expected         = c(2L, 3L)  # cutoff=15 → idx 2 (last=20); cutoff=25 → idx 3 (last=30)
    )
  )

  N_CASES      <- length(cases)
  MAX_PATIENTS <- max(sapply(cases, function(x) x$n_patients))
  MAX_CUTOFFS  <- max(sapply(cases, function(x) x$n_cutoffs))

  # Rectangularize into padded 2D matrices
  sorted_last_visit_mat <- matrix(0L, N_CASES, MAX_PATIENTS)
  cutoffs_mat           <- matrix(0L, N_CASES, MAX_CUTOFFS)
  n_patients_vec        <- integer(N_CASES)
  n_cutoffs_vec         <- integer(N_CASES)

  for (i in seq_along(cases)) {
    cc <- cases[[i]]
    n_patients_vec[i] <- cc$n_patients
    n_cutoffs_vec[i]  <- cc$n_cutoffs
    sorted_last_visit_mat[i, seq_len(cc$n_patients)] <- cc$sorted_last_visit
    cutoffs_mat[i, seq_len(cc$n_cutoffs)]            <- cc$cutoffs
  }

  stan_data <- list(
    N_CASES    = N_CASES,
    MAX_PATIENTS = MAX_PATIENTS,
    MAX_CUTOFFS  = MAX_CUTOFFS,
    n_patients = n_patients_vec,
    n_cutoffs  = n_cutoffs_vec,
    sorted_last_visit_calendar_day = sorted_last_visit_mat,
    cutoff_calendar_day            = cutoffs_mat
  )

  fit      <- test_stan_function(
    "tests/testthat/stan/test_get_oos_patients_idx_all.stan",
    data = stan_data
  )
  draws_df <- posterior::as_draws_df(fit$draws())
  get_val  <- function(var, ...) get_stan_val(draws_df, var, ...)

  for (i in seq_along(cases)) {
    cc <- cases[[i]]
    for (c in seq_len(cc$n_cutoffs)) {
      expect_equal(
        get_val("oos_idx_out", i, c),
        cc$expected[c],
        label = paste0("case ", i, " (", cc$case_name, ") cutoff ", c)
      )
    }
  }
})
```

- [ ] **Step 2: Run test — verify it fails (Stan file not found)**

```bash
cd /mnt/code/.worktrees/karim/tests
Rscript -e 'testthat::test_file("tests/testthat/test-stan-get_oos_patients_idx.R")'
```

Expected: FAIL — `test_get_oos_patients_idx_all.stan` not found.

- [ ] **Step 3: Create the Stan harness**

Create `tests/testthat/stan/test_get_oos_patients_idx_all.stan`:

```stan
functions {
  #include "util.stanfunctions"
  #include "pos.stanfunctions"
  #include "lfo.stanfunctions"
}

data {
  int<lower=1> N_CASES;
  int<lower=1> MAX_PATIENTS;
  int<lower=1> MAX_CUTOFFS;
  array[N_CASES] int<lower=1> n_patients;
  array[N_CASES] int<lower=1> n_cutoffs;
  array[N_CASES, MAX_PATIENTS] int sorted_last_visit_calendar_day;
  array[N_CASES, MAX_CUTOFFS]  int cutoff_calendar_day;
}

generated quantities {
  array[N_CASES, MAX_CUTOFFS] int oos_idx_out = rep_array(0, N_CASES, MAX_CUTOFFS);

  for (case in 1:N_CASES) {
    int np = n_patients[case];
    int nc = n_cutoffs[case];
    array[nc] int result = get_oos_patients_idx(
      sorted_last_visit_calendar_day[case, 1:np],
      cutoff_calendar_day[case, 1:nc]
    );
    for (c in 1:nc) {
      oos_idx_out[case, c] = result[c];
    }
  }
}
```

- [ ] **Step 4: Run test — verify it passes**

```bash
Rscript -e 'testthat::test_file("tests/testthat/test-stan-get_oos_patients_idx.R")'
```

Expected: PASS — 8 assertions (sum of n_cutoffs across 6 cases).

- [ ] **Step 5: Commit**

```bash
git add tests/testthat/stan/test_get_oos_patients_idx_all.stan \
        tests/testthat/test-stan-get_oos_patients_idx.R
git commit -m "test: add get_oos_patients_idx tests — 8 assertions"
```

---

## Task 2: `truncate_at_max_time` — fix Stan harness + R test

**Files:**
- Modify: `tests/testthat/stan/test_pfs_functions_all.stan`
- Create: `tests/testthat/helper-pfs.R`
- Create: `tests/testthat/test-stan-pfs-truncate_at_max_time.R`

**Context:** `test_pfs_functions_all.stan` already exists with a correct data-driven section, but has a broken inline section (deprecated `int arr[N]` syntax, 5-argument call) that prevents it from compiling. The fix is a full rewrite of just the Stan file with the inline section removed.

---

- [ ] **Step 1: Create `helper-pfs.R` with the R reference**

Create `tests/testthat/helper-pfs.R`:

```r
#' R reference implementation of Stan truncate_at_max_time.
#'
#' Mirrors: pfs.stanfunctions::truncate_at_max_time(pfs, right_censored, max_time)
#'
#' Semantics: if pfs[i] > max_time[i], truncate to max_time[i] and censor.
#' Boundary (pfs[i] == max_time[i]) is NOT truncated.
#'
#' @param pfs           Integer vector of PFS times.
#' @param right_censored Integer vector (0/1) of censoring flags.
#' @param max_time      Integer vector of per-patient cutoff times.
#' @return Named list with elements `pfs` and `right_censored`.
r_truncate_at_max_time <- function(pfs, right_censored, max_time) {
  n     <- length(pfs)
  o_pfs <- integer(n)
  o_rc  <- integer(n)
  for (i in seq_len(n)) {
    if (pfs[i] > max_time[i]) {
      o_pfs[i] <- max_time[i]
      o_rc[i]  <- 1L
    } else {
      o_pfs[i] <- pfs[i]
      o_rc[i]  <- right_censored[i]
    }
  }
  list(pfs = o_pfs, right_censored = o_rc)
}

#' R port of Stan find_first(all, what, n_succ).
#'
#' Returns the 1-based start index of the first run of n_succ consecutive
#' elements of `all` that are each present in `what`. Returns 0 if not found.
r_find_first <- function(all, what, n_succ) {
  n <- length(all)
  if (n < n_succ) return(0L)
  if (length(what) == 0L) stop("r_find_first: 'what' cannot be empty")
  for (i in seq_len(n - n_succ + 1L)) {
    matches <- 0L
    for (j in seq_len(n_succ)) {
      if (all[i + j - 1L] %in% what) {
        matches <- matches + 1L
      } else {
        break
      }
    }
    if (matches == n_succ) return(i)
  }
  return(0L)
}

#' R port of Stan map_idx_to_week.
#'
#' Maps a 1-based combined (obs + forecast) index to an actual week number.
#' idx <= n_obs → curr_visits[idx]; else → forecast_time[idx - n_obs].
#' idx == 0 → max_all_t (censoring sentinel).
r_map_idx_to_week <- function(idx, curr_visits, forecast_time, max_all_t) {
  n_obs <- length(curr_visits)
  if (idx == 0L) return(max_all_t)
  if (idx <= n_obs) return(curr_visits[idx])
  forecast_idx <- idx - n_obs
  if (forecast_idx < 1L || forecast_idx > length(forecast_time)) return(max_all_t)
  forecast_time[forecast_idx]
}

#' R reference implementation of Stan find_first_week.
#'
#' Mirrors: pfs.stanfunctions::find_first_week(arr, values, min_run_length,
#'                                              curr_visits, forecast_time, max_all_t)
#'
#' @param arr            Integer vector (RECIST codes or similar).
#' @param values         Integer vector of target values to search for.
#' @param min_run_length Minimum consecutive matches required.
#' @param curr_visits    Observed treatment visit weeks (may be integer(0)).
#' @param forecast_time  Forecast visit weeks.
#' @param max_all_t      Censoring sentinel returned when not found.
#' @return Named list with elements `week` (int) and `right_censored` (0/1 int).
r_find_first_week <- function(arr, values, min_run_length,
                              curr_visits, forecast_time, max_all_t) {
  idx <- r_find_first(arr, values, min_run_length)
  if (idx == 0L) {
    list(week = max_all_t, right_censored = 1L)
  } else {
    list(
      week          = r_map_idx_to_week(idx, curr_visits, forecast_time, max_all_t),
      right_censored = 0L
    )
  }
}

#' R reference implementation of Stan find_first_forecast_week.
#'
#' Delegates to r_find_first_week with empty curr_visits.
#'
#' @param arr            Integer vector.
#' @param values         Integer vector of target values.
#' @param min_run_length Minimum consecutive matches.
#' @param forecast_time  Forecast visit weeks.
#' @param max_all_t      Censoring sentinel.
#' @return Named list with `week` and `right_censored`.
r_find_first_forecast_week <- function(arr, values, min_run_length,
                                       forecast_time, max_all_t) {
  r_find_first_week(arr, values, min_run_length, integer(0), forecast_time, max_all_t)
}
```

- [ ] **Step 2: Write the failing R test for `truncate_at_max_time`**

Create `tests/testthat/test-stan-pfs-truncate_at_max_time.R`:

```r
library(testthat)
library(cmdstanr)
library(posterior)
source(here::here("tests/testthat/helper-stan.R"))
source(here::here("tests/testthat/helper-pfs.R"))

test_that("truncate_at_max_time: all cases produce correct truncated pfs and right_censored", {
  cases <- list(
    list(
      case_name     = "no_truncation_needed",
      pfs           = c(5L, 10L, 15L),
      right_censored = c(0L, 0L, 0L),
      max_time      = c(20L, 20L, 20L)
    ),
    list(
      case_name     = "all_truncated",
      pfs           = c(25L, 30L, 40L),
      right_censored = c(0L, 0L, 0L),
      max_time      = c(20L, 20L, 20L)
    ),
    list(
      case_name     = "mixed_truncation_and_censoring",
      pfs           = c(5L, 15L, 25L, 35L),
      right_censored = c(0L, 1L, 0L, 1L),
      max_time      = c(10L, 20L, 20L, 30L)
    ),
    list(
      case_name     = "exact_boundary_not_truncated",
      pfs           = c(10L, 20L, 30L),
      right_censored = c(0L, 0L, 0L),
      max_time      = c(10L, 20L, 30L)
    ),
    list(
      case_name     = "already_censored_no_truncation",
      pfs           = c(8L, 12L),
      right_censored = c(1L, 1L),
      max_time      = c(20L, 20L)
    ),
    list(
      case_name     = "already_censored_truncation_needed",
      pfs           = c(25L, 30L),
      right_censored = c(1L, 1L),
      max_time      = c(20L, 20L)
    )
  )

  N_CASES <- length(cases)
  N       <- max(sapply(cases, function(x) length(x$pfs)))

  # Pad all cases to N patients using dummy values that won't be truncated
  pfs_mat <- matrix(1L, N_CASES, N)
  rc_mat  <- matrix(0L, N_CASES, N)
  mt_mat  <- matrix(100L, N_CASES, N)

  for (i in seq_along(cases)) {
    cc <- cases[[i]]
    n  <- length(cc$pfs)
    pfs_mat[i, seq_len(n)] <- cc$pfs
    rc_mat[i,  seq_len(n)] <- cc$right_censored
    mt_mat[i,  seq_len(n)] <- cc$max_time
  }

  stan_data <- list(
    N_CASES        = N_CASES,
    N              = N,
    pfs            = pfs_mat,
    right_censored = rc_mat,
    max_time       = mt_mat
  )

  fit      <- test_stan_function(
    "tests/testthat/stan/test_pfs_functions_all.stan",
    data = stan_data
  )
  draws_df <- posterior::as_draws_df(fit$draws())
  get_val  <- function(var, ...) get_stan_val(draws_df, var, ...)

  for (i in seq_along(cases)) {
    cc       <- cases[[i]]
    expected <- r_truncate_at_max_time(cc$pfs, cc$right_censored, cc$max_time)
    n        <- length(cc$pfs)
    for (j in seq_len(n)) {
      expect_equal(
        get_val("out_pfs", i, j),
        expected$pfs[j],
        label = paste0("case ", i, " (", cc$case_name, ") out_pfs[", j, "]")
      )
      expect_equal(
        get_val("out_right_censored", i, j),
        expected$right_censored[j],
        label = paste0("case ", i, " (", cc$case_name, ") out_right_censored[", j, "]")
      )
    }
  }
})
```

- [ ] **Step 3: Run — verify it fails (Stan harness does not compile)**

```bash
Rscript -e 'testthat::test_file("tests/testthat/test-stan-pfs-truncate_at_max_time.R")'
```

Expected: FAIL — `test_pfs_functions_all.stan` syntax error (deprecated `int arr[N]` syntax in inline section).

- [ ] **Step 4: Fix the Stan harness (overwrite with clean version)**

Replace `tests/testthat/stan/test_pfs_functions_all.stan` with the following (removes the broken inline section, also removes `expected_pfs`/`expected_right_censored`/`n_match` since expected values are computed in R):

```stan
functions {
  #include "pfs.stanfunctions"
  #include "util.stanfunctions"
  #include "pos.stanfunctions"
}

data {
  int<lower=1> N_CASES;
  int<lower=1> N;
  array[N_CASES, N] int pfs;
  array[N_CASES, N] int right_censored;
  array[N_CASES, N] int max_time;
}

generated quantities {
  array[N_CASES, N] int out_pfs;
  array[N_CASES, N] int out_right_censored;

  for (case in 1:N_CASES) {
    tuple(array[N] int, array[N] int) res =
      truncate_at_max_time(pfs[case], right_censored[case], max_time[case]);
    out_pfs[case]            = res.1;
    out_right_censored[case] = res.2;
  }
}
```

- [ ] **Step 5: Run — verify it passes**

```bash
Rscript -e 'testthat::test_file("tests/testthat/test-stan-pfs-truncate_at_max_time.R")'
```

Expected: PASS — 36 assertions (6 cases × avg 3 patients × 2 outputs).

- [ ] **Step 6: Commit**

```bash
git add tests/testthat/stan/test_pfs_functions_all.stan \
        tests/testthat/helper-pfs.R \
        tests/testthat/test-stan-pfs-truncate_at_max_time.R
git commit -m "test: add truncate_at_max_time tests + fix broken Stan harness"
```

---

## Task 3: `find_first_week` + `find_first_forecast_week`

**Files:**
- Create: `tests/testthat/test-helper-pfs.R`
- Create: `tests/testthat/stan/test_pfs_find_first_week_all.stan`
- Create: `tests/testthat/test-stan-pfs-find_first_week.R`

**Context:** `helper-pfs.R` already exists from Task 2. `find_first_week` calls `find_first` (from `util.stanfunctions`) to get a 1-based index, then maps it to a week via `map_idx_to_week`. `find_first_forecast_week` is a thin wrapper that passes `zeros_int_array(0)` as `curr_visits`.

---

- [ ] **Step 1: Write R-only unit tests for `helper-pfs.R`**

Create `tests/testthat/test-helper-pfs.R`:

```r
library(testthat)
source(here::here("tests/testthat/helper-pfs.R"))

test_that("r_find_first: basic find_first semantics", {
  # found at position 3
  expect_equal(r_find_first(c(3L, 3L, 4L), c(4L), 1L), 3L)
  # found at position 1
  expect_equal(r_find_first(c(4L, 3L, 3L), c(4L), 1L), 1L)
  # not found
  expect_equal(r_find_first(c(3L, 3L, 3L), c(4L), 1L), 0L)
  # run of 2 required
  expect_equal(r_find_first(c(3L, 4L, 4L, 3L), c(4L), 2L), 2L)
  # run of 2, not found (only one 4)
  expect_equal(r_find_first(c(3L, 4L, 3L, 4L), c(4L), 2L), 0L)
  # multiple values in 'what'
  expect_equal(r_find_first(c(3L, 1L, 3L), c(1L, 2L), 1L), 2L)
})

test_that("r_map_idx_to_week: maps obs vs forecast correctly", {
  curr   <- c(4L, 8L, 12L)
  fore   <- c(16L, 20L)
  max_t  <- 52L
  # idx in obs range
  expect_equal(r_map_idx_to_week(1L, curr, fore, max_t), 4L)
  expect_equal(r_map_idx_to_week(3L, curr, fore, max_t), 12L)
  # idx in forecast range
  expect_equal(r_map_idx_to_week(4L, curr, fore, max_t), 16L)
  expect_equal(r_map_idx_to_week(5L, curr, fore, max_t), 20L)
  # idx == 0 → censoring sentinel
  expect_equal(r_map_idx_to_week(0L, curr, fore, max_t), 52L)
  # out of forecast range → censoring sentinel
  expect_equal(r_map_idx_to_week(6L, curr, fore, max_t), 52L)
})

test_that("r_find_first_week: week and right_censored", {
  arr  <- c(3L, 3L, 4L)
  curr <- c(4L, 8L, 12L)
  fore <- c(16L, 20L)
  mt   <- 52L

  # found in obs → week = curr_visits[3] = 12, rc = 0
  res <- r_find_first_week(arr, c(4L), 1L, curr, fore, mt)
  expect_equal(res$week, 12L)
  expect_equal(res$right_censored, 0L)

  # not found → week = max_all_t, rc = 1
  res <- r_find_first_week(c(3L, 3L, 3L), c(4L), 1L, curr, fore, mt)
  expect_equal(res$week, 52L)
  expect_equal(res$right_censored, 1L)

  # found in forecast (idx=4 > n_obs=3) → week = fore[1] = 16
  arr4 <- c(3L, 3L, 3L, 4L)
  res  <- r_find_first_week(arr4, c(4L), 1L, curr, fore, mt)
  expect_equal(res$week, 16L)
  expect_equal(res$right_censored, 0L)

  # empty curr_visits → idx maps directly to forecast
  res <- r_find_first_week(c(3L, 4L), c(4L), 1L, integer(0), c(8L, 16L), mt)
  expect_equal(res$week, 16L)   # idx=2, forecast_idx=2
  expect_equal(res$right_censored, 0L)
})

test_that("r_find_first_forecast_week: delegates with empty curr_visits", {
  # found at forecast position 1
  res <- r_find_first_forecast_week(c(4L), c(4L), 1L, c(8L), 52L)
  expect_equal(res$week, 8L)
  expect_equal(res$right_censored, 0L)

  # not found
  res <- r_find_first_forecast_week(c(3L, 3L), c(4L), 1L, c(8L, 16L), 52L)
  expect_equal(res$week, 52L)
  expect_equal(res$right_censored, 1L)
})
```

- [ ] **Step 2: Run R-only tests — verify they pass**

```bash
Rscript -e 'testthat::test_file("tests/testthat/test-helper-pfs.R")'
```

Expected: PASS — 20 assertions. (No Stan required.)

- [ ] **Step 3: Create the Stan harness**

Create `tests/testthat/stan/test_pfs_find_first_week_all.stan`:

```stan
functions {
  #include "pfs.stanfunctions"
  #include "util.stanfunctions"
  #include "pos.stanfunctions"
}

// Harness tests two functions:
//   find_first_week(arr, values, min_run_length, curr_visits, forecast_time, max_all_t)
//   find_first_forecast_week(arr, values, min_run_length, forecast_time, max_all_t)
//
// N_CASES / N_FC_CASES use padded 2D arrays; case-specific lengths are passed via n_* arrays.
// Slicing to 1:0 produces a valid empty array when n_obs[case] = 0.

data {
  // ---- find_first_week cases ----
  int<lower=1> N_CASES;
  int<lower=1> MAX_ARR;
  int<lower=1> MAX_VALS;
  int<lower=1> MAX_OBS;
  int<lower=1> MAX_FORE;
  array[N_CASES] int n_arr;
  array[N_CASES] int n_vals;
  array[N_CASES] int n_obs;
  array[N_CASES] int n_fore;
  array[N_CASES] int min_run_length;
  array[N_CASES] int max_all_t;
  array[N_CASES, MAX_ARR]  int arr;
  array[N_CASES, MAX_VALS] int vals;
  array[N_CASES, MAX_OBS]  int curr_visits;
  array[N_CASES, MAX_FORE] int forecast_time;

  // ---- find_first_forecast_week cases ----
  int<lower=1> N_FC_CASES;
  int<lower=1> MAX_ARR_FC;
  int<lower=1> MAX_VALS_FC;
  int<lower=1> MAX_FORE_FC;
  array[N_FC_CASES] int n_arr_fc;
  array[N_FC_CASES] int n_vals_fc;
  array[N_FC_CASES] int n_fore_fc;
  array[N_FC_CASES] int min_run_fc;
  array[N_FC_CASES] int max_all_t_fc;
  array[N_FC_CASES, MAX_ARR_FC]  int arr_fc;
  array[N_FC_CASES, MAX_VALS_FC] int vals_fc;
  array[N_FC_CASES, MAX_FORE_FC] int forecast_time_fc;
}

generated quantities {
  // find_first_week outputs
  array[N_CASES] int out_week;
  array[N_CASES] int out_rc;

  for (case in 1:N_CASES) {
    (out_week[case], out_rc[case]) = find_first_week(
      arr[case, 1:n_arr[case]],
      vals[case, 1:n_vals[case]],
      min_run_length[case],
      curr_visits[case, 1:n_obs[case]],
      forecast_time[case, 1:n_fore[case]],
      max_all_t[case]
    );
  }

  // find_first_forecast_week outputs
  array[N_FC_CASES] int out_week_fc;
  array[N_FC_CASES] int out_rc_fc;

  for (case in 1:N_FC_CASES) {
    (out_week_fc[case], out_rc_fc[case]) = find_first_forecast_week(
      arr_fc[case, 1:n_arr_fc[case]],
      vals_fc[case, 1:n_vals_fc[case]],
      min_run_fc[case],
      forecast_time_fc[case, 1:n_fore_fc[case]],
      max_all_t_fc[case]
    );
  }
}
```

- [ ] **Step 4: Write the R-Stan integration test**

Create `tests/testthat/test-stan-pfs-find_first_week.R`:

```r
library(testthat)
library(cmdstanr)
library(posterior)
source(here::here("tests/testthat/helper-stan.R"))
source(here::here("tests/testthat/helper-pfs.R"))

test_that("find_first_week and find_first_forecast_week: Stan matches R reference", {

  # ---- find_first_week cases ----
  # Each case: arr, values, min_run_length, curr_visits (may be empty), forecast_time, max_all_t
  ffw_cases <- list(
    list(
      name           = "found_in_obs_visits",
      arr            = c(3L, 3L, 4L),
      values         = c(4L),
      min_run_length = 1L,
      curr_visits    = c(4L, 8L, 12L),
      forecast_time  = c(16L, 20L),
      max_all_t      = 52L
      # find_first=3, idx=3 ≤ n_obs=3 → week=12, rc=0
    ),
    list(
      name           = "found_in_forecast",
      arr            = c(3L, 3L, 3L, 4L),
      values         = c(4L),
      min_run_length = 1L,
      curr_visits    = c(4L, 8L, 12L),
      forecast_time  = c(16L, 20L),
      max_all_t      = 52L
      # find_first=4, idx=4 > n_obs=3 → forecast_idx=1 → week=16, rc=0
    ),
    list(
      name           = "not_found",
      arr            = c(3L, 3L, 3L),
      values         = c(4L),
      min_run_length = 1L,
      curr_visits    = c(4L, 8L, 12L),
      forecast_time  = c(16L, 20L),
      max_all_t      = 52L
      # find_first=0 → week=52, rc=1
    ),
    list(
      name           = "min_run_length_2",
      arr            = c(3L, 4L, 4L, 3L),
      values         = c(4L),
      min_run_length = 2L,
      curr_visits    = c(4L, 8L, 12L, 16L),
      forecast_time  = c(20L),
      max_all_t      = 52L
      # find_first([3,4,4,3],[4],2)=2 → week=curr_visits[2]=8, rc=0
    ),
    list(
      name           = "found_at_first_position",
      arr            = c(4L, 3L, 3L),
      values         = c(4L),
      min_run_length = 1L,
      curr_visits    = c(4L, 8L, 12L),
      forecast_time  = c(16L),
      max_all_t      = 52L
      # find_first=1 → week=4, rc=0
    ),
    list(
      name           = "multiple_values_in_what",
      arr            = c(3L, 1L, 3L),
      values         = c(1L, 2L),
      min_run_length = 1L,
      curr_visits    = c(4L, 8L, 12L),
      forecast_time  = c(16L),
      max_all_t      = 52L
      # arr[2]=1 ∈ {1,2} → find_first=2 → week=8, rc=0
    ),
    list(
      name           = "empty_curr_visits_match_in_forecast",
      arr            = c(3L, 4L),
      values         = c(4L),
      min_run_length = 1L,
      curr_visits    = integer(0),   # n_obs = 0
      forecast_time  = c(8L, 16L),
      max_all_t      = 52L
      # find_first=2, n_obs=0 → forecast_idx=2 → week=16, rc=0
    )
  )

  # ---- find_first_forecast_week cases ----
  fffc_cases <- list(
    list(
      name           = "found_at_position_2",
      arr            = c(3L, 4L),
      values         = c(4L),
      min_run_length = 1L,
      forecast_time  = c(8L, 16L),
      max_all_t      = 52L
      # find_first=2, forecast_idx=2 → week=16, rc=0
    ),
    list(
      name           = "found_at_first_position",
      arr            = c(4L),
      values         = c(4L),
      min_run_length = 1L,
      forecast_time  = c(8L),
      max_all_t      = 52L
      # find_first=1 → week=8, rc=0
    ),
    list(
      name           = "not_found",
      arr            = c(3L, 3L),
      values         = c(4L),
      min_run_length = 1L,
      forecast_time  = c(8L, 16L),
      max_all_t      = 52L
      # find_first=0 → week=52, rc=1
    )
  )

  # ---- Compute R expected values ----
  ffw_expected <- lapply(ffw_cases, function(cc) {
    r_find_first_week(
      cc$arr, cc$values, cc$min_run_length,
      cc$curr_visits, cc$forecast_time, cc$max_all_t
    )
  })
  fffc_expected <- lapply(fffc_cases, function(cc) {
    r_find_first_forecast_week(
      cc$arr, cc$values, cc$min_run_length,
      cc$forecast_time, cc$max_all_t
    )
  })

  # ---- Build Stan data ----
  N_CASES    <- length(ffw_cases)
  N_FC_CASES <- length(fffc_cases)

  MAX_ARR  <- max(sapply(ffw_cases,  function(x) length(x$arr)))
  MAX_VALS <- max(sapply(ffw_cases,  function(x) length(x$values)))
  MAX_OBS  <- max(1L, max(sapply(ffw_cases,  function(x) length(x$curr_visits))))
  MAX_FORE <- max(sapply(ffw_cases,  function(x) length(x$forecast_time)))

  MAX_ARR_FC  <- max(sapply(fffc_cases, function(x) length(x$arr)))
  MAX_VALS_FC <- max(sapply(fffc_cases, function(x) length(x$values)))
  MAX_FORE_FC <- max(sapply(fffc_cases, function(x) length(x$forecast_time)))

  # Padding helper
  pad <- function(x, n, fill = 0L) {
    if (length(x) >= n) return(x[seq_len(n)])
    c(x, rep(fill, n - length(x)))
  }

  arr_mat  <- matrix(0L, N_CASES, MAX_ARR)
  vals_mat <- matrix(0L, N_CASES, MAX_VALS)
  obs_mat  <- matrix(0L, N_CASES, MAX_OBS)
  fore_mat <- matrix(0L, N_CASES, MAX_FORE)
  n_arr    <- integer(N_CASES)
  n_vals   <- integer(N_CASES)
  n_obs    <- integer(N_CASES)
  n_fore   <- integer(N_CASES)
  min_run  <- integer(N_CASES)
  max_t    <- integer(N_CASES)

  for (i in seq_along(ffw_cases)) {
    cc         <- ffw_cases[[i]]
    n_arr[i]   <- length(cc$arr)
    n_vals[i]  <- length(cc$values)
    n_obs[i]   <- length(cc$curr_visits)
    n_fore[i]  <- length(cc$forecast_time)
    min_run[i] <- cc$min_run_length
    max_t[i]   <- cc$max_all_t
    arr_mat[i,  seq_len(n_arr[i])]  <- cc$arr
    vals_mat[i, seq_len(n_vals[i])] <- cc$values
    if (n_obs[i] > 0) obs_mat[i, seq_len(n_obs[i])] <- cc$curr_visits
    fore_mat[i, seq_len(n_fore[i])] <- cc$forecast_time
  }

  arr_fc_mat  <- matrix(0L, N_FC_CASES, MAX_ARR_FC)
  vals_fc_mat <- matrix(0L, N_FC_CASES, MAX_VALS_FC)
  fore_fc_mat <- matrix(0L, N_FC_CASES, MAX_FORE_FC)
  n_arr_fc    <- integer(N_FC_CASES)
  n_vals_fc   <- integer(N_FC_CASES)
  n_fore_fc   <- integer(N_FC_CASES)
  min_run_fc  <- integer(N_FC_CASES)
  max_t_fc    <- integer(N_FC_CASES)

  for (i in seq_along(fffc_cases)) {
    cc              <- fffc_cases[[i]]
    n_arr_fc[i]     <- length(cc$arr)
    n_vals_fc[i]    <- length(cc$values)
    n_fore_fc[i]    <- length(cc$forecast_time)
    min_run_fc[i]   <- cc$min_run_length
    max_t_fc[i]     <- cc$max_all_t
    arr_fc_mat[i,   seq_len(n_arr_fc[i])]  <- cc$arr
    vals_fc_mat[i,  seq_len(n_vals_fc[i])] <- cc$values
    fore_fc_mat[i,  seq_len(n_fore_fc[i])] <- cc$forecast_time
  }

  stan_data <- list(
    N_CASES    = N_CASES,
    MAX_ARR    = MAX_ARR,
    MAX_VALS   = MAX_VALS,
    MAX_OBS    = MAX_OBS,
    MAX_FORE   = MAX_FORE,
    n_arr      = n_arr,
    n_vals     = n_vals,
    n_obs      = n_obs,
    n_fore     = n_fore,
    min_run_length = min_run,
    max_all_t  = max_t,
    arr        = arr_mat,
    vals       = vals_mat,
    curr_visits     = obs_mat,
    forecast_time   = fore_mat,
    N_FC_CASES      = N_FC_CASES,
    MAX_ARR_FC      = MAX_ARR_FC,
    MAX_VALS_FC     = MAX_VALS_FC,
    MAX_FORE_FC     = MAX_FORE_FC,
    n_arr_fc        = n_arr_fc,
    n_vals_fc       = n_vals_fc,
    n_fore_fc       = n_fore_fc,
    min_run_fc      = min_run_fc,
    max_all_t_fc    = max_t_fc,
    arr_fc          = arr_fc_mat,
    vals_fc         = vals_fc_mat,
    forecast_time_fc = fore_fc_mat
  )

  fit      <- test_stan_function(
    "tests/testthat/stan/test_pfs_find_first_week_all.stan",
    data = stan_data
  )
  draws_df <- posterior::as_draws_df(fit$draws())
  get_val  <- function(var, ...) get_stan_val(draws_df, var, ...)

  # ---- Assert find_first_week ----
  for (i in seq_along(ffw_cases)) {
    cc  <- ffw_cases[[i]]
    exp <- ffw_expected[[i]]
    expect_equal(
      get_val("out_week", i),
      exp$week,
      label = paste0("find_first_week case ", i, " (", cc$name, ") week")
    )
    expect_equal(
      get_val("out_rc", i),
      exp$right_censored,
      label = paste0("find_first_week case ", i, " (", cc$name, ") right_censored")
    )
  }

  # ---- Assert find_first_forecast_week ----
  for (i in seq_along(fffc_cases)) {
    cc  <- fffc_cases[[i]]
    exp <- fffc_expected[[i]]
    expect_equal(
      get_val("out_week_fc", i),
      exp$week,
      label = paste0("find_first_forecast_week case ", i, " (", cc$name, ") week")
    )
    expect_equal(
      get_val("out_rc_fc", i),
      exp$right_censored,
      label = paste0("find_first_forecast_week case ", i, " (", cc$name, ") right_censored")
    )
  }
})
```

- [ ] **Step 5: Run — verify it passes**

```bash
Rscript -e 'testthat::test_file("tests/testthat/test-stan-pfs-find_first_week.R")'
```

Expected: PASS — 20 assertions (14 for `find_first_week` + 6 for `find_first_forecast_week`).

- [ ] **Step 6: Run full test suite — confirm no regressions**

```bash
Rscript -e 'testthat::test_dir("tests/testthat")'
```

Expected: all previously passing tests still pass; new total ≥ 762 assertions (698 + ~64 new).

- [ ] **Step 7: Commit**

```bash
git add tests/testthat/test-helper-pfs.R \
        tests/testthat/stan/test_pfs_find_first_week_all.stan \
        tests/testthat/test-stan-pfs-find_first_week.R
git commit -m "test: add find_first_week + find_first_forecast_week tests — 20 assertions"
```

---

## Summary

| Task | New files | Assertions added |
|---|---|---|
| 1. `get_oos_patients_idx` | 2 | ~8 |
| 2. `truncate_at_max_time` | 2 + 1 fix | ~36 |
| 3. `find_first_week` + forecast variant | 3 + helper | ~20 + 20 (R-only) |
| **Total** | **7 files** | **~84 new** |

All Phase 1 LFO pipeline functions covered: `get_oos_patients_idx`, `truncate_at_max_time`, `find_first_week`, `find_first_forecast_week`.
