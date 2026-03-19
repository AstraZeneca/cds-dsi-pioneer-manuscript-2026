# Test Design Improvements — `get_testing_visit_week_bounds` Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Eliminate three classes of test-design bugs: (1) copy-paste draws-extraction boilerplate that breaks on multi-dim Stan arrays, (2) hand-computed expected values that drift from the Stan implementation, (3) untested `last_testing_visit_week`/`testing_end_idx` output fields.

**Architecture:** A shared `get_stan_val()` helper goes into the existing `helper-stan.R`. A new `helper-lfo.R` provides a pure-R mirror of `get_testing_visit_week_bounds()` that is the single source of truth for expected values. Both Stan test files are refactored to source `helper-lfo.R` and compute all expected values programmatically rather than hardcoding matrices.

**Tech Stack:** R, testthat, cmdstanr, posterior

---

## File Map

| File | Action | Responsibility |
|------|--------|---------------|
| `tests/testthat/helper-stan.R` | Modify | Add `get_stan_val(draws_df, var, ...)` |
| `tests/testthat/helper-lfo.R` | **Create** | R mirror of `get_testing_visit_week_bounds()` |
| `tests/testthat/test-helper-lfo.R` | **Create** | Unit tests for the R mirror |
| `tests/testthat/test-stan-get_testing_visit_week_bounds.R` | Modify | Use helpers; compute all expected values; validate all 4 output arrays |
| `tests/testthat/test-multi-patient-get_testing_visit_week_bounds.R` | Modify | Use helpers; compute expected values |

---

## Task 1 — Add `get_stan_val()` to `helper-stan.R`

**Files:**
- Modify: `tests/testthat/helper-stan.R`

`draws_array` is always 3D `[iterations, chains, variables]`. Multi-dimensional Stan arrays flatten into variable names like `x[1,2,3]`. A shared helper removes all ad-hoc extraction code.

- [ ] **Step 1: Append the helper to `helper-stan.R` after `create_mock_km_data`**

```r
#' Extract a scalar value from a posterior draws data frame by named Stan variable.
#'
#' draws_array is always 3D [iterations, chains, variables]. Multi-dimensional
#' Stan arrays are stored as named scalars, e.g. x[1,2,3]. Use this helper
#' instead of numeric indexing which breaks when variable names sort
#' alphabetically (e.g. x[10,1] sorts before x[2,1]).
#'
#' @param draws_df A draws_df from posterior::as_draws_df(fit$draws())
#' @param var Stan variable name (string)
#' @param ... Integer indices (one per dimension); omit for scalar variables
#' @return Numeric scalar (first iteration value)
get_stan_val <- function(draws_df, var, ...) {
  indices <- c(...)
  vname   <- if (length(indices) == 0) var else sprintf("%s[%s]", var, paste(indices, collapse = ","))
  as.numeric(draws_df[[vname]][1])
}
```

- [ ] **Step 2: Run tests to confirm no regression**

```bash
Rscript -e 'testthat::test_dir("tests/testthat")'
```
Expected: `[ FAIL 0 | PASS 334 ]`

- [ ] **Step 3: Commit**

```bash
git add tests/testthat/helper-stan.R
git commit -m "feat: add get_stan_val() helper to helper-stan.R for robust named draws access"
```

---

## Task 2 — Create `helper-lfo.R`: pure-R mirror of `get_testing_visit_week_bounds()`

**Files:**
- Create: `tests/testthat/helper-lfo.R`
- Create: `tests/testthat/test-helper-lfo.R`

This R function mirrors the Stan logic exactly, serving as the single authoritative source for expected values in all `get_testing_visit_week_bounds` tests.

Key semantics faithfully replicated from `stan/lfo.stanfunctions` lines 241–343:
- `calendar_date_to_study_date = cutoff - patient_entry + 1` (note +1)
- `first_testing_visit_week`: first visit where `day > study_day` (strict greater-than)
- Patients at sort positions `< oos_patient_idx[n]` are skipped; their arrays stay at 0
- `last_testing_visit_week[n, m, i]`: last visit where `day <= upper_study_day` (inclusive), computed for `m > n` only using the reversed-walk loop
- Initialisation: `first = 0`, `last = min(t_patient_visits_week)`, `start = 0`, `end = 0`
- **Pre-condition**: both `t_patient_visits_week` and `t_patient_visits_day` must be strictly ascending within each patient's block. The function validates this and stops if violated.
- **Pre-condition**: all patients must have at least one visit (`patient_visit_pos[i+1] > patient_visit_pos[i]`).

- [ ] **Step 1: Create `helper-lfo.R`**

```r
# tests/testthat/helper-lfo.R
#
# Pure R mirror of get_testing_visit_week_bounds() from stan/lfo.stanfunctions.
# Used to compute reference expected values for Stan tests.
#
# Pre-conditions (same as Stan asserts):
#   - t_patient_visits_week and t_patient_visits_day are strictly ascending
#     within each patient's visit block
#   - Every patient has at least 1 visit (patient_visit_pos[i+1] > patient_visit_pos[i])
#   - patient_visit_pos has length n_patients + 1


#' R mirror of Stan's get_testing_visit_week_bounds().
#'
#' @param oos_patient_idx     integer vector [n_cutoffs]
#' @param last_visit_calendar_day_sort_idx integer vector [n_patients]
#' @param cutoff_calendar_day integer vector [n_cutoffs]
#' @param patient_calendar_day integer vector [n_patients]
#' @param t_patient_visits_week integer vector [n_visits]
#' @param t_patient_visits_day  integer vector [n_visits]
#' @param patient_visit_pos     integer vector [n_patients + 1]
#'
#' @return named list:
#'   first_testing_visit_week  matrix[n_cutoffs, n_patients]
#'   last_testing_visit_week   array [n_cutoffs, n_cutoffs, n_patients]
#'   testing_start_idx         matrix[n_cutoffs, n_patients]
#'   testing_end_idx           array [n_cutoffs, n_cutoffs, n_patients]
r_get_testing_visit_week_bounds <- function(
  oos_patient_idx,
  last_visit_calendar_day_sort_idx,
  cutoff_calendar_day,
  patient_calendar_day,
  t_patient_visits_week,
  t_patient_visits_day,
  patient_visit_pos
) {
  n_patients <- length(patient_calendar_day)
  n_futures  <- length(oos_patient_idx)
  min_all_t  <- min(t_patient_visits_week)

  # Validate pre-conditions (mirrors Stan's assert_strict_ascending calls)
  for (i in seq_len(n_patients)) {
    vs <- patient_visit_pos[i]
    ve <- patient_visit_pos[i + 1L] - 1L
    stopifnot(ve >= vs)  # at least one visit
    if (ve > vs) {
      weeks <- t_patient_visits_week[vs:ve]
      days  <- t_patient_visits_day[vs:ve]
      stopifnot(all(diff(weeks) > 0), all(diff(days) > 0))
    }
  }

  first_testing_visit_week <- matrix(0L, nrow = n_futures, ncol = n_patients)
  last_testing_visit_week  <- array(min_all_t, dim = c(n_futures, n_futures, n_patients))
  testing_start_idx        <- matrix(0L, nrow = n_futures, ncol = n_patients)
  testing_end_idx          <- array(0L,  dim = c(n_futures, n_futures, n_patients))

  for (n in seq_len(n_futures)) {
    n_curr_patients <- n_patients - oos_patient_idx[n] + 1L
    curr_patients   <- last_visit_calendar_day_sort_idx[oos_patient_idx[n]:n_patients]

    # Lower cutoff: patient-specific study day for cutoff n (+1 is intentional)
    lower_sd <- cutoff_calendar_day[n] - patient_calendar_day[curr_patients] + 1L

    for (i_idx in seq_len(n_curr_patients)) {
      i           <- curr_patients[i_idx]
      visit_start <- patient_visit_pos[i]
      visit_end   <- patient_visit_pos[i + 1L] - 1L
      n_visits    <- visit_end - visit_start + 1L
      p_days      <- t_patient_visits_day[visit_start:visit_end]
      p_weeks     <- t_patient_visits_week[visit_start:visit_end]

      # First visit STRICTLY AFTER lower cutoff (mirrors Stan while loop)
      t_idx <- 1L
      while (t_idx <= n_visits && p_days[t_idx] <= max(0L, lower_sd[i_idx])) {
        t_idx <- t_idx + 1L
      }
      if (t_idx <= n_visits) {
        first_testing_visit_week[n, i] <- p_weeks[t_idx]
        testing_start_idx[n, i]        <- visit_start + t_idx - 1L
      }
    }

    # Upper cutoff: last visit on or before cutoff m — only computed for m > n
    # NOTE: guard required because R's (n+1):n gives c(n+1, n), not an empty range
    if (n < n_futures) {
      for (m in (n + 1L):n_futures) {
        upper_sd <- cutoff_calendar_day[m] - patient_calendar_day[curr_patients] + 1L

        for (i_idx in seq_len(n_curr_patients)) {
          i           <- curr_patients[i_idx]
          visit_start <- patient_visit_pos[i]
          visit_end   <- patient_visit_pos[i + 1L] - 1L
          n_visits    <- visit_end - visit_start + 1L
          rev_days    <- rev(t_patient_visits_day[visit_start:visit_end])
          rev_weeks   <- rev(t_patient_visits_week[visit_start:visit_end])

          # Walk reversed visits until day <= upper cutoff (mirrors Stan while loop)
          t_idx <- 1L
          while (t_idx <= n_visits && rev_days[t_idx] > max(0L, upper_sd[i_idx])) {
            t_idx <- t_idx + 1L
          }
          if (t_idx <= n_visits) {
            last_testing_visit_week[n, m, i] <- rev_weeks[t_idx]
            testing_end_idx[n, m, i]         <- visit_end - t_idx + 1L
          }
        }
      }
    }
  }

  list(
    first_testing_visit_week = first_testing_visit_week,
    last_testing_visit_week  = last_testing_visit_week,
    testing_start_idx        = testing_start_idx,
    testing_end_idx          = testing_end_idx
  )
}
```

- [ ] **Step 2: Create `test-helper-lfo.R`**

```r
library(testthat)
source(here::here("tests/testthat/helper-lfo.R"))

test_that("r_get_testing_visit_week_bounds matches manually verified results", {
  # ── Case A: single patient, single cutoff ─────────────────────────────────
  # study_day = 110 - 100 + 1 = 11
  # visits days 7 (wk1), 21 (wk3) — first > 11: day 21 → wk3, idx 2
  resA <- r_get_testing_visit_week_bounds(
    oos_patient_idx                  = c(1L),
    last_visit_calendar_day_sort_idx = c(1L),
    cutoff_calendar_day              = c(110L),
    patient_calendar_day             = c(100L),
    t_patient_visits_week            = c(1L, 3L),
    t_patient_visits_day             = c(7L, 21L),
    patient_visit_pos                = c(1L, 3L)
  )
  expect_equal(resA$first_testing_visit_week, matrix(3L, 1, 1))
  expect_equal(resA$testing_start_idx,        matrix(2L, 1, 1))
  # Single cutoff: last_testing_visit_week stays at min_all_t = 1; testing_end_idx = 0
  expect_equal(resA$last_testing_visit_week, array(1L, dim = c(1, 1, 1)))
  expect_equal(resA$testing_end_idx,         array(0L, dim = c(1, 1, 1)))

  # ── Case B: single patient, two cutoffs ───────────────────────────────────
  # visits days 7,14,21,28,35 (weeks 1,2,3,4,5); entry=100
  # cutoff 1 (cal=110): sd=11, first > 11 → day14 wk2, idx2
  # cutoff 2 (cal=120): sd=21, first > 21 → day28 wk4, idx4
  # last[n=1, m=2]: last day <= 21 → day21 wk3, end_idx=3
  # diagonal [1,1] and [2,2] stay at 0 for testing_end_idx
  resB <- r_get_testing_visit_week_bounds(
    oos_patient_idx                  = c(1L, 1L),
    last_visit_calendar_day_sort_idx = c(1L),
    cutoff_calendar_day              = c(110L, 120L),
    patient_calendar_day             = c(100L),
    t_patient_visits_week            = c(1L, 2L, 3L, 4L, 5L),
    t_patient_visits_day             = c(7L, 14L, 21L, 28L, 35L),
    patient_visit_pos                = c(1L, 6L)
  )
  expect_equal(resB$first_testing_visit_week, matrix(c(2L, 4L), nrow = 2, ncol = 1))
  expect_equal(resB$testing_start_idx,        matrix(c(2L, 4L), nrow = 2, ncol = 1))
  expect_equal(resB$last_testing_visit_week[1, 2, 1], 3L)
  expect_equal(resB$testing_end_idx[1, 2, 1],         3L)
  # Diagonal cells must stay at their initialised values
  expect_equal(resB$testing_end_idx[1, 1, 1], 0L)
  expect_equal(resB$testing_end_idx[2, 2, 1], 0L)

  # ── Case C: no visit after lower cutoff (retention path) ─────────────────
  # entry=100, cutoff=140 (sd=41), all visits at days 7,14,21,28,35 (<= 41)
  # first_testing_visit_week should stay at 0 (no visit after cutoff)
  resC <- r_get_testing_visit_week_bounds(
    oos_patient_idx                  = c(1L),
    last_visit_calendar_day_sort_idx = c(1L),
    cutoff_calendar_day              = c(140L),
    patient_calendar_day             = c(100L),
    t_patient_visits_week            = c(1L, 2L, 3L, 4L, 5L),
    t_patient_visits_day             = c(7L, 14L, 21L, 28L, 35L),
    patient_visit_pos                = c(1L, 6L)
  )
  expect_equal(resC$first_testing_visit_week[1, 1], 0L)
  expect_equal(resC$testing_start_idx[1, 1],        0L)

  # ── Case D: oos_patient_idx skips patient 1 for cutoff 2 ─────────────────
  # Two patients; cutoff 2 has oos_idx=2 so only patient 2 is tested
  resD <- r_get_testing_visit_week_bounds(
    oos_patient_idx                  = c(1L, 2L),
    last_visit_calendar_day_sort_idx = c(1L, 2L),
    cutoff_calendar_day              = c(110L, 120L),
    patient_calendar_day             = c(100L, 105L),
    t_patient_visits_week            = c(1L, 3L, 2L, 4L),
    t_patient_visits_day             = c(7L, 21L, 14L, 28L),
    patient_visit_pos                = c(1L, 3L, 5L)
  )
  # cutoff 2, patient 1 not tested → stays at 0
  expect_equal(resD$first_testing_visit_week[2, 1], 0L)
  expect_equal(resD$testing_start_idx[2, 1],        0L)
  # cutoff 2, patient 2: sd = 120-105+1 = 16, first > 16: day28 wk4, global idx 4
  expect_equal(resD$first_testing_visit_week[2, 2], 4L)
  expect_equal(resD$testing_start_idx[2, 2],        4L)
})
```

- [ ] **Step 3: Run the helper unit tests**

```bash
Rscript -e 'testthat::test_file("tests/testthat/test-helper-lfo.R")'
```
Expected: `[ FAIL 0 | PASS ... ]`

- [ ] **Step 4: Commit**

```bash
git add tests/testthat/helper-lfo.R tests/testthat/test-helper-lfo.R
git commit -m "feat: add r_get_testing_visit_week_bounds() R reference + unit tests"
```

---

## Task 3 — Refactor `test-stan-get_testing_visit_week_bounds.R`

**Files:**
- Modify: `tests/testthat/test-stan-get_testing_visit_week_bounds.R`

Three changes:
1. Source `helper-lfo.R` at the top
2. Remove all hardcoded `expected_*` fields from the 16 non-error test cases
3. Compute expected values programmatically; extend validation loop to all 4 arrays

- [ ] **Step 1: Add `source(helper-lfo.R)` at the top of the file**

After line 3 (`source(here::here("tests/testthat/helper-stan.R"))`), add:
```r
source(here::here("tests/testthat/helper-lfo.R"))
```

- [ ] **Step 2: Remove `expected_*` fields from all non-error test cases**

The 16 non-error cases (Cases 1–5, 8–18) each end with four fields:
```r
expected_first_testing_visit_week = <matrix or NULL>,
expected_testing_start_idx = <matrix or NULL>,
expected_last_testing_visit_week = NULL,
expected_testing_end_idx = NULL
```
Delete all four fields from every non-error case. Cases 6 and 7 (`should_error = TRUE`) never had these fields — leave them untouched.

After editing, the last field in every non-error case should be `patient_visit_pos = c(...)`.

- [ ] **Step 3: Compute expected values before the `tryCatch` block**

Insert the following block between line 594 (`cat("*** This test pushes...")`) and line 596 (`tryCatch(`):

```r
  # Compute expected values using the R reference implementation.
  # This is the authoritative source of truth — not hand-computed matrices.
  expected_by_case <- lapply(seq_len(N_cases), function(ci) {
    case <- test_cases[[ci]]
    if (case$should_error) return(NULL)
    r_get_testing_visit_week_bounds(
      oos_patient_idx                  = case$oos_patient_idx,
      last_visit_calendar_day_sort_idx = case$last_visit_calendar_day_sort_idx,
      cutoff_calendar_day              = case$cutoff_calendar_day,
      patient_calendar_day             = case$patient_calendar_day,
      t_patient_visits_week            = case$t_patient_visits,
      t_patient_visits_day             = case$t_patient_visits_day,
      patient_visit_pos                = case$patient_visit_pos
    )
  })
```

- [ ] **Step 4: Replace the draws extraction and validation loop**

Replace lines 606–659 (from `draws_df <- posterior::...` through the closing brace of the valid-case loop) with:

```r
      draws_df <- posterior::as_draws_df(fit$draws())
      get_val  <- function(var, ...) get_stan_val(draws_df, var, ...)

      cat("Successfully extracted test results\n")

      # Validate all valid cases — all 4 output arrays
      for (case_idx in valid_cases) {
        case     <- test_cases[[case_idx]]
        expected <- expected_by_case[[case_idx]]
        lbl      <- paste("Case", case_idx, case$case_name)

        expect_equal(get_val("case_status", case_idx), 0, label = paste(lbl, "no error"))

        for (n in 1:case$n_cutoffs) {
          for (i in 1:case$n_patients) {
            expect_equal(
              get_val("first_testing_visit_week", case_idx, n, i),
              expected$first_testing_visit_week[n, i],
              label = paste(lbl, "first_testing_visit_week n", n, "i", i)
            )
            expect_equal(
              get_val("testing_start_idx", case_idx, n, i),
              expected$testing_start_idx[n, i],
              label = paste(lbl, "testing_start_idx n", n, "i", i)
            )
            for (m in 1:case$n_cutoffs) {
              expect_equal(
                get_val("last_testing_visit_week", case_idx, n, m, i),
                expected$last_testing_visit_week[n, m, i],
                label = paste(lbl, "last_testing_visit_week n", n, "m", m, "i", i)
              )
              expect_equal(
                get_val("testing_end_idx", case_idx, n, m, i),
                expected$testing_end_idx[n, m, i],
                label = paste(lbl, "testing_end_idx n", n, "m", m, "i", i)
              )
            }
          }
        }
      }
```

- [ ] **Step 5: Run the full test suite**

```bash
Rscript -e 'testthat::test_dir("tests/testthat")'
```
Expected: `[ FAIL 0 | PASS > 334 ]` (higher pass count due to new assertions)

- [ ] **Step 6: Commit**

```bash
git add tests/testthat/test-stan-get_testing_visit_week_bounds.R
git commit -m "refactor: compute expected values from R reference; validate all 4 output arrays"
```

---

## Task 4 — Refactor `test-multi-patient-get_testing_visit_week_bounds.R`

**Files:**
- Modify: `tests/testthat/test-multi-patient-get_testing_visit_week_bounds.R`

Same three changes as Task 3 applied to the per-case-loop structure.

- [ ] **Step 1: Add sources at the top**

After line 4 (`source(here::here("tests/testthat/helper-stan.R"))`), add:
```r
source(here::here("tests/testthat/helper-lfo.R"))
```

- [ ] **Step 2: Remove hardcoded `expected_*` fields from all 3 cases**

Each case ends with:
```r
expected_first_testing_visit_week = matrix(...),
expected_testing_start_idx = matrix(...)
```
Delete both fields. The last field in each case should be `patient_visit_pos = c(...)`.

- [ ] **Step 3: Replace the extraction + validation block inside the `for (case_idx ...)` loop**

The current block (lines 155–228) runs a Stan model, extracts draws, then validates. Replace lines 160–228 (from `# Extract results...` through `cat("Case", case_idx, "completed successfully!\n")`) with:

```r
    # Compute expected values from R reference
    expected <- r_get_testing_visit_week_bounds(
      oos_patient_idx                  = case$oos_patient_idx,
      last_visit_calendar_day_sort_idx = case$last_visit_calendar_day_sort_idx,
      cutoff_calendar_day              = case$cutoff_calendar_day,
      patient_calendar_day             = case$patient_calendar_day,
      t_patient_visits_week            = case$t_patient_visits,
      t_patient_visits_day             = case$t_patient_visits_day,
      patient_visit_pos                = case$patient_visit_pos
    )

    draws_df <- posterior::as_draws_df(fit$draws())
    get_val  <- function(var, ...) get_stan_val(draws_df, var, ...)

    cat("Stan execution successful!\n")

    case_status <- get_val("case_status", 1)
    expect_equal(case_status, 0, label = paste("Case", case_idx, "should not error"))

    for (n in 1:case$n_cutoffs) {
      for (i in 1:case$n_patients) {
        expect_equal(
          get_val("first_testing_visit_week", 1L, n, i),
          expected$first_testing_visit_week[n, i],
          label = paste("Case", case_idx, case$case_name, "first_testing_visit_week n=", n, "i=", i)
        )
        expect_equal(
          get_val("testing_start_idx", 1L, n, i),
          expected$testing_start_idx[n, i],
          label = paste("Case", case_idx, case$case_name, "testing_start_idx n=", n, "i=", i)
        )
        for (m in 1:case$n_cutoffs) {
          expect_equal(
            get_val("last_testing_visit_week", 1L, n, m, i),
            expected$last_testing_visit_week[n, m, i],
            label = paste("Case", case_idx, case$case_name, "last_testing_visit_week n=", n, "m=", m, "i=", i)
          )
          expect_equal(
            get_val("testing_end_idx", 1L, n, m, i),
            expected$testing_end_idx[n, m, i],
            label = paste("Case", case_idx, case$case_name, "testing_end_idx n=", n, "m=", m, "i=", i)
          )
        }
      }
    }

    cat("Case", case_idx, "completed successfully!\n")
```

- [ ] **Step 4: Run the full test suite**

```bash
Rscript -e 'testthat::test_dir("tests/testthat")'
```
Expected: `[ FAIL 0 | PASS > 334 ]`

- [ ] **Step 5: Commit**

```bash
git add tests/testthat/test-multi-patient-get_testing_visit_week_bounds.R
git commit -m "refactor: use R reference for expected values in multi-patient test; validate all 4 output arrays"
```

---

## Verification

After all tasks complete:

```bash
Rscript -e 'testthat::test_dir("tests/testthat")'
```

Pass count should be noticeably higher than 334 (new assertions for `last_testing_visit_week` and `testing_end_idx` across 19 cases × n_cutoffs × n_cutoffs × n_patients combinations). Zero failures, zero warnings.
