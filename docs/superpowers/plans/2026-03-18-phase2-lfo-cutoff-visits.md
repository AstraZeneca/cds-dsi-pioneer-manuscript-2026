# Phase 2: LFO `cutoff_visits` + `fine_cutoff_visits` Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add Stan unit tests for `cutoff_visits` and `fine_cutoff_visits` from `stan/lfo.stanfunctions`, completing LFO pipeline test coverage.

**Architecture:** Same oracle-driven pattern as Phase 1: add R reference implementations to `helper-lfo.R`, write R-only unit tests, then write a Stan harness + R-Stan integration test. Both functions determine the last visit/measurement before a calendar cutoff for each patient — `cutoff_visits` works on patient visits (5-tuple output), `fine_cutoff_visits` on tumor measurements (3-tuple output).

**Tech Stack:** R + testthat, Stan (cmdstanr), posterior, `helper-stan.R`, `helper-lfo.R`.

---

## Background: The Two Functions

Both live in `stan/lfo.stanfunctions`.

```
cutoff_visits(cutoff_calendar_day, patient_calendar_day, t_patient_visits,
              t_patient_visits_day, patient_visit_pos)
  → tuple(last_visit_day[], last_visit_week[], last_visit_calendar_day[],
           last_visit_week_observed[], cutoff_last_visit_idx[])

fine_cutoff_visits(cutoff_calendar_day, patient_calendar_day,
                   t_measure, t_day_measure, patient_tumor_measure_pos)
  → tuple(last_visit_day[], last_visit_week[], last_visit_calendar_day[])
```

**Key semantics (both functions):**
- `study_day = cutoff_calendar_day - patient_calendar_day[i] + 1`
- Finds the last visit/measure with `day <= study_day` (inclusive boundary)
- If no visit before cutoff: `last_visit_day=0`, `last_visit_week=0`
- `last_visit_calendar_day` = calendar day of patient's **final** visit (ignores cutoff — always set)

**`cutoff_visits`-only semantics:**
- `last_visit_week_observed = last_visit_week > 0` (0 or 1 flag)
- `cutoff_last_visit_idx`: index in `t_patient_visits` (global flat array) of the last visit with week ≤ last_visit_week; 0 if no visit before cutoff. Walks original (unsorted) visit array — requires visits to be in ascending week order.

**Difference in loop style:**
- `cutoff_visits`: `t_idx=0; while(t_idx < n && day[sort[t_idx+1]] <= study_day) t_idx++`
- `fine_cutoff_visits`: `t_idx=1; while(t_idx <= n && day[sort[t_idx]] <= study_day) t_idx++`
Both yield the same semantics; the off-by-one accounting differs.

---

## File Structure

| File | Action | Purpose |
|---|---|---|
| `tests/testthat/helper-lfo.R` | Modify | Add `r_cutoff_visits`, `r_fine_cutoff_visits` |
| `tests/testthat/test-helper-lfo-cutoff-visits.R` | Create | R-only unit tests for both R references |
| `tests/testthat/stan/test_lfo_cutoff_visits_all.stan` | Create | Stan harness for both functions |
| `tests/testthat/test-stan-lfo-cutoff-visits.R` | Create | R-Stan integration test |

---

## Task 1: R references + R-only unit tests

**Files:**
- Modify: `tests/testthat/helper-lfo.R` (append to end of file)
- Create: `tests/testthat/test-helper-lfo-cutoff-visits.R`

---

- [ ] **Step 1: Append R references to `helper-lfo.R`**

Append these two functions after the existing `r_get_testing_visit_week_bounds` function:

```r
#' R mirror of Stan's cutoff_visits().
#'
#' Finds the last visit before or on the study-day equivalent of cutoff_calendar_day
#' for each patient. Visits may be in any order in the array — they are sorted by day.
#'
#' @param cutoff_calendar_day Scalar int: global calendar cutoff day.
#' @param patient_calendar_day Integer vector [n_patients]: study entry calendar days.
#' @param t_patient_visits Integer vector: visit week numbers (flat, all patients).
#' @param t_patient_visits_day Integer vector: visit day numbers (flat, all patients).
#' @param patient_visit_pos Integer vector [n_patients+1]: position array.
#' @return Named list with 5 integer vectors of length n_patients.
r_cutoff_visits <- function(cutoff_calendar_day, patient_calendar_day,
                            t_patient_visits, t_patient_visits_day,
                            patient_visit_pos) {
  n_patients             <- length(patient_calendar_day)
  last_visit_day         <- integer(n_patients)
  last_visit_week        <- integer(n_patients)
  last_visit_calendar_day <- integer(n_patients)
  last_visit_week_observed <- integer(n_patients)
  cutoff_last_visit_idx  <- integer(n_patients)

  for (i in seq_len(n_patients)) {
    vs  <- patient_visit_pos[i]
    ve  <- patient_visit_pos[i + 1L] - 1L
    n_v <- ve - vs + 1L
    if (n_v <= 0L) next

    study_day <- cutoff_calendar_day - patient_calendar_day[i] + 1L
    wks  <- t_patient_visits[vs:ve]
    days <- t_patient_visits_day[vs:ve]
    srt  <- order(days)  # mirrors sort_indices_asc(patient_t_visits_day)

    # t_idx=0; while(t_idx < n && day[sort[t_idx+1]] <= study_day) t_idx++
    t_idx <- 0L
    while (t_idx < n_v && days[srt[t_idx + 1L]] <= study_day) {
      t_idx <- t_idx + 1L
    }

    if (t_idx > 0L) {
      last_visit_day[i]  <- days[srt[t_idx]]
      last_visit_week[i] <- wks[srt[t_idx]]
    }

    last_visit_week_observed[i] <- as.integer(last_visit_week[i] > 0L)
    # Final visit calendar day (independent of cutoff — always the last visit)
    last_visit_calendar_day[i] <- patient_calendar_day[i] + days[srt[n_v]] - 1L

    if (last_visit_week[i] > 0L) {
      # Walk original (unsorted) array to find last index with week <= last_visit_week
      j <- vs
      while (j <= ve && t_patient_visits[j] <= last_visit_week[i]) {
        cutoff_last_visit_idx[i] <- j
        j <- j + 1L
      }
    }
  }

  list(
    last_visit_day            = last_visit_day,
    last_visit_week           = last_visit_week,
    last_visit_calendar_day   = last_visit_calendar_day,
    last_visit_week_observed  = last_visit_week_observed,
    cutoff_last_visit_idx     = cutoff_last_visit_idx
  )
}

#' R mirror of Stan's fine_cutoff_visits().
#'
#' Like r_cutoff_visits but for tumor measurements (not visits). Returns 3 arrays.
#'
#' @param cutoff_calendar_day Scalar int: global calendar cutoff day.
#' @param patient_calendar_day Integer vector [n_patients].
#' @param t_measure Integer vector: measurement week numbers (flat).
#' @param t_day_measure Integer vector: measurement day numbers (flat).
#' @param patient_tumor_measure_pos Integer vector [n_patients+1]: position array.
#' @return Named list with 3 integer vectors of length n_patients.
r_fine_cutoff_visits <- function(cutoff_calendar_day, patient_calendar_day,
                                  t_measure, t_day_measure,
                                  patient_tumor_measure_pos) {
  n_patients              <- length(patient_calendar_day)
  last_visit_day          <- integer(n_patients)
  last_visit_week         <- integer(n_patients)
  last_visit_calendar_day <- integer(n_patients)

  for (i in seq_len(n_patients)) {
    ms  <- patient_tumor_measure_pos[i]
    me  <- patient_tumor_measure_pos[i + 1L] - 1L
    n_m <- me - ms + 1L

    if (n_m <= 0L) {
      # last_visit_calendar_day[i] stays 0
      next
    }

    study_day <- cutoff_calendar_day - patient_calendar_day[i] + 1L
    wks  <- t_measure[ms:me]
    days <- t_day_measure[ms:me]
    srt  <- order(days)

    # t_idx=1; while(t_idx <= n && day[sort[t_idx]] <= study_day) t_idx++
    t_idx <- 1L
    while (t_idx <= n_m && days[srt[t_idx]] <= study_day) {
      t_idx <- t_idx + 1L
    }

    if (t_idx > 1L) {
      last_visit_day[i]  <- days[srt[t_idx - 1L]]
      last_visit_week[i] <- wks[srt[t_idx - 1L]]
    }

    last_visit_calendar_day[i] <- patient_calendar_day[i] + days[srt[n_m]] - 1L
  }

  list(
    last_visit_day          = last_visit_day,
    last_visit_week         = last_visit_week,
    last_visit_calendar_day = last_visit_calendar_day
  )
}
```

- [ ] **Step 2: Create R-only unit tests**

Create `tests/testthat/test-helper-lfo-cutoff-visits.R`:

```r
library(testthat)

test_that("r_cutoff_visits: correct last visit before cutoff", {
  # Case 1: one visit before cutoff
  res <- r_cutoff_visits(
    cutoff_calendar_day   = 120L,
    patient_calendar_day  = c(100L),
    t_patient_visits      = c(2L),
    t_patient_visits_day  = c(14L),
    patient_visit_pos     = c(1L, 2L)
  )
  expect_equal(res$last_visit_day,           c(14L))
  expect_equal(res$last_visit_week,          c(2L))
  expect_equal(res$last_visit_calendar_day,  c(113L))  # 100 + 14 - 1
  expect_equal(res$last_visit_week_observed, c(1L))
  expect_equal(res$cutoff_last_visit_idx,    c(1L))

  # Case 2: all visits after cutoff → sentinels
  res2 <- r_cutoff_visits(
    cutoff_calendar_day  = 110L,
    patient_calendar_day = c(100L),
    t_patient_visits     = c(2L),
    t_patient_visits_day = c(14L),
    patient_visit_pos    = c(1L, 2L)
  )
  expect_equal(res2$last_visit_day,           c(0L))
  expect_equal(res2$last_visit_week,          c(0L))
  expect_equal(res2$last_visit_calendar_day,  c(113L))  # still the final visit cal day
  expect_equal(res2$last_visit_week_observed, c(0L))
  expect_equal(res2$cutoff_last_visit_idx,    c(0L))

  # Case 3: exact boundary (day == study_day) is included
  res3 <- r_cutoff_visits(
    cutoff_calendar_day  = 113L,   # study_day = 113-100+1 = 14
    patient_calendar_day = c(100L),
    t_patient_visits     = c(2L),
    t_patient_visits_day = c(14L),
    patient_visit_pos    = c(1L, 2L)
  )
  expect_equal(res3$last_visit_week, c(2L))

  # Case 4: multiple visits, partial
  # visits wk=[1,2,4,6], day=[7,14,28,42]; cutoff study_day=31
  # last before cutoff: day28 wk4; final visit: day42 → cal=141
  res4 <- r_cutoff_visits(
    cutoff_calendar_day  = 130L,
    patient_calendar_day = c(100L),
    t_patient_visits     = c(1L, 2L, 4L, 6L),
    t_patient_visits_day = c(7L, 14L, 28L, 42L),
    patient_visit_pos    = c(1L, 5L)
  )
  expect_equal(res4$last_visit_week,         c(4L))
  expect_equal(res4$last_visit_day,          c(28L))
  expect_equal(res4$last_visit_calendar_day, c(141L))  # 100 + 42 - 1
  expect_equal(res4$cutoff_last_visit_idx,   c(3L))    # wk 1,2,4 all ≤ 4; wk 6 > 4

  # Case 5: two patients
  res5 <- r_cutoff_visits(
    cutoff_calendar_day  = 125L,
    patient_calendar_day = c(100L, 110L),
    t_patient_visits     = c(1L, 3L, 2L, 4L),
    t_patient_visits_day = c(7L, 21L, 14L, 28L),
    patient_visit_pos    = c(1L, 3L, 5L)
  )
  # P1: study_day=26, day7+day21 both ≤ 26 → last_wk=3; final=day21 → cal=120
  expect_equal(res5$last_visit_week[1],         3L)
  expect_equal(res5$last_visit_calendar_day[1], 120L)  # 100 + 21 - 1
  expect_equal(res5$cutoff_last_visit_idx[1],   2L)    # pos 1→wk1≤3, pos 2→wk3≤3
  # P2: study_day=16, day14≤16, day28>16 → last_wk=2; final=day28 → cal=137
  expect_equal(res5$last_visit_week[2],         2L)
  expect_equal(res5$last_visit_calendar_day[2], 137L)  # 110 + 28 - 1
  expect_equal(res5$cutoff_last_visit_idx[2],   3L)    # pos 3→wk2≤2; pos 4→wk4>2
})

test_that("r_fine_cutoff_visits: correct last measure before cutoff", {
  # Case 1: one measure before cutoff
  res <- r_fine_cutoff_visits(
    cutoff_calendar_day        = 120L,
    patient_calendar_day       = c(100L),
    t_measure                  = c(2L),
    t_day_measure              = c(14L),
    patient_tumor_measure_pos  = c(1L, 2L)
  )
  expect_equal(res$last_visit_day,          c(14L))
  expect_equal(res$last_visit_week,         c(2L))
  expect_equal(res$last_visit_calendar_day, c(113L))

  # Case 2: all measures after cutoff → sentinels
  res2 <- r_fine_cutoff_visits(
    cutoff_calendar_day       = 110L,
    patient_calendar_day      = c(100L),
    t_measure                 = c(2L),
    t_day_measure             = c(14L),
    patient_tumor_measure_pos = c(1L, 2L)
  )
  expect_equal(res2$last_visit_day,          c(0L))
  expect_equal(res2$last_visit_week,         c(0L))
  expect_equal(res2$last_visit_calendar_day, c(113L))  # final measure cal day

  # Case 3: exact boundary included
  res3 <- r_fine_cutoff_visits(
    cutoff_calendar_day       = 113L,
    patient_calendar_day      = c(100L),
    t_measure                 = c(2L),
    t_day_measure             = c(14L),
    patient_tumor_measure_pos = c(1L, 2L)
  )
  expect_equal(res3$last_visit_week, c(2L))

  # Case 4: multiple measures, partial
  res4 <- r_fine_cutoff_visits(
    cutoff_calendar_day       = 130L,
    patient_calendar_day      = c(100L),
    t_measure                 = c(1L, 2L, 4L, 6L),
    t_day_measure             = c(7L, 14L, 28L, 42L),
    patient_tumor_measure_pos = c(1L, 5L)
  )
  expect_equal(res4$last_visit_week,         c(4L))
  expect_equal(res4$last_visit_day,          c(28L))
  expect_equal(res4$last_visit_calendar_day, c(141L))

  # Case 5: two patients
  res5 <- r_fine_cutoff_visits(
    cutoff_calendar_day       = 125L,
    patient_calendar_day      = c(100L, 110L),
    t_measure                 = c(1L, 3L, 2L, 4L),
    t_day_measure             = c(7L, 21L, 14L, 28L),
    patient_tumor_measure_pos = c(1L, 3L, 5L)
  )
  expect_equal(res5$last_visit_week[1],         3L)
  expect_equal(res5$last_visit_calendar_day[1], 120L)
  expect_equal(res5$last_visit_week[2],         2L)
  expect_equal(res5$last_visit_calendar_day[2], 137L)
})
```

- [ ] **Step 3: Run R-only tests**

```bash
cd /mnt/code/.worktrees/karim/tests
Rscript -e 'testthat::test_file("tests/testthat/test-helper-lfo-cutoff-visits.R")'
```

Expected: PASS — 28 assertions (13 for `r_cutoff_visits` + 15 for `r_fine_cutoff_visits`).

- [ ] **Step 4: Commit**

```bash
git add tests/testthat/helper-lfo.R tests/testthat/test-helper-lfo-cutoff-visits.R
git commit -m "test: add r_cutoff_visits and r_fine_cutoff_visits R references + unit tests"
```

---

## Task 2: Stan harness + R-Stan integration test

**Files:**
- Create: `tests/testthat/stan/test_lfo_cutoff_visits_all.stan`
- Create: `tests/testthat/test-stan-lfo-cutoff-visits.R`

---

- [ ] **Step 1: Create the Stan harness**

Create `tests/testthat/stan/test_lfo_cutoff_visits_all.stan`:

```stan
functions {
  #include "util.stanfunctions"
  #include "pos.stanfunctions"
  #include "lfo.stanfunctions"
}

// Tests cutoff_visits (5-tuple) and fine_cutoff_visits (3-tuple) from lfo.stanfunctions.
// Each set of cases uses padded 2D arrays; per-case actual sizes are passed via n_* arrays.

data {
  // ---- cutoff_visits cases ----
  int<lower=1> N_CV;
  int<lower=1> MAX_P_CV;
  int<lower=1> MAX_V_CV;
  array[N_CV] int n_patients_cv;
  array[N_CV] int n_visits_cv;
  array[N_CV] int cutoff_cal_day_cv;
  array[N_CV, MAX_P_CV]   int patient_cal_day_cv;
  array[N_CV, MAX_V_CV]   int t_visits_cv;
  array[N_CV, MAX_V_CV]   int t_visits_day_cv;
  array[N_CV, MAX_P_CV+1] int visit_pos_cv;

  // ---- fine_cutoff_visits cases ----
  int<lower=1> N_FCV;
  int<lower=1> MAX_P_FCV;
  int<lower=1> MAX_M_FCV;
  array[N_FCV] int n_patients_fcv;
  array[N_FCV] int n_measures_fcv;
  array[N_FCV] int cutoff_cal_day_fcv;
  array[N_FCV, MAX_P_FCV]   int patient_cal_day_fcv;
  array[N_FCV, MAX_M_FCV]   int t_measure_fcv;
  array[N_FCV, MAX_M_FCV]   int t_day_measure_fcv;
  array[N_FCV, MAX_P_FCV+1] int measure_pos_fcv;
}

generated quantities {
  // cutoff_visits outputs (initialize with sentinel -999)
  array[N_CV, MAX_P_CV] int cv_last_visit_day       = rep_array(-999, N_CV, MAX_P_CV);
  array[N_CV, MAX_P_CV] int cv_last_visit_week      = rep_array(-999, N_CV, MAX_P_CV);
  array[N_CV, MAX_P_CV] int cv_last_cal_day         = rep_array(-999, N_CV, MAX_P_CV);
  array[N_CV, MAX_P_CV] int cv_week_observed        = rep_array(-999, N_CV, MAX_P_CV);
  array[N_CV, MAX_P_CV] int cv_cutoff_last_idx      = rep_array(-999, N_CV, MAX_P_CV);

  for (c in 1:N_CV) {
    int np = n_patients_cv[c];
    int nv = n_visits_cv[c];
    array[np] int ld, lw, lcd, lwo, clvi;
    (ld, lw, lcd, lwo, clvi) = cutoff_visits(
      cutoff_cal_day_cv[c],
      patient_cal_day_cv[c, 1:np],
      t_visits_cv[c, 1:nv],
      t_visits_day_cv[c, 1:nv],
      visit_pos_cv[c, 1:(np+1)]
    );
    for (i in 1:np) {
      cv_last_visit_day[c, i]  = ld[i];
      cv_last_visit_week[c, i] = lw[i];
      cv_last_cal_day[c, i]    = lcd[i];
      cv_week_observed[c, i]   = lwo[i];
      cv_cutoff_last_idx[c, i] = clvi[i];
    }
  }

  // fine_cutoff_visits outputs
  array[N_FCV, MAX_P_FCV] int fcv_last_visit_day  = rep_array(-999, N_FCV, MAX_P_FCV);
  array[N_FCV, MAX_P_FCV] int fcv_last_visit_week = rep_array(-999, N_FCV, MAX_P_FCV);
  array[N_FCV, MAX_P_FCV] int fcv_last_cal_day    = rep_array(-999, N_FCV, MAX_P_FCV);

  for (c in 1:N_FCV) {
    int np = n_patients_fcv[c];
    int nm = n_measures_fcv[c];
    array[np] int ld, lw, lcd;
    (ld, lw, lcd) = fine_cutoff_visits(
      cutoff_cal_day_fcv[c],
      patient_cal_day_fcv[c, 1:np],
      t_measure_fcv[c, 1:nm],
      t_day_measure_fcv[c, 1:nm],
      measure_pos_fcv[c, 1:(np+1)]
    );
    for (i in 1:np) {
      fcv_last_visit_day[c, i]  = ld[i];
      fcv_last_visit_week[c, i] = lw[i];
      fcv_last_cal_day[c, i]    = lcd[i];
    }
  }
}
```

- [ ] **Step 2: Write the R-Stan integration test**

Create `tests/testthat/test-stan-lfo-cutoff-visits.R`:

```r
library(testthat)
library(cmdstanr)
library(posterior)

test_that("cutoff_visits and fine_cutoff_visits: Stan matches R reference", {

  # ---- cutoff_visits cases ----
  cv_cases <- list(
    list(
      name             = "single_patient_one_visit_before",
      n_patients       = 1L,
      cutoff           = 120L,
      patient_cal      = c(100L),
      t_visits         = c(2L),
      t_visits_day     = c(14L),
      visit_pos        = c(1L, 2L)
    ),
    list(
      name             = "single_patient_no_visits_before",
      n_patients       = 1L,
      cutoff           = 110L,
      patient_cal      = c(100L),
      t_visits         = c(2L),
      t_visits_day     = c(14L),
      visit_pos        = c(1L, 2L)
    ),
    list(
      name             = "exact_boundary_included",
      n_patients       = 1L,
      cutoff           = 113L,   # study_day = 14 = day of visit
      patient_cal      = c(100L),
      t_visits         = c(2L),
      t_visits_day     = c(14L),
      visit_pos        = c(1L, 2L)
    ),
    list(
      name             = "multiple_visits_partial",
      n_patients       = 1L,
      cutoff           = 130L,
      patient_cal      = c(100L),
      t_visits         = c(1L, 2L, 4L, 6L),
      t_visits_day     = c(7L, 14L, 28L, 42L),
      visit_pos        = c(1L, 5L)
    ),
    list(
      name             = "two_patients",
      n_patients       = 2L,
      cutoff           = 125L,
      patient_cal      = c(100L, 110L),
      t_visits         = c(1L, 3L, 2L, 4L),
      t_visits_day     = c(7L, 21L, 14L, 28L),
      visit_pos        = c(1L, 3L, 5L)
    )
  )

  # ---- fine_cutoff_visits cases ----
  fcv_cases <- list(
    list(
      name           = "single_patient_one_measure_before",
      n_patients     = 1L,
      cutoff         = 120L,
      patient_cal    = c(100L),
      t_measure      = c(2L),
      t_day_measure  = c(14L),
      measure_pos    = c(1L, 2L)
    ),
    list(
      name           = "single_patient_no_measures_before",
      n_patients     = 1L,
      cutoff         = 110L,
      patient_cal    = c(100L),
      t_measure      = c(2L),
      t_day_measure  = c(14L),
      measure_pos    = c(1L, 2L)
    ),
    list(
      name           = "exact_boundary_included",
      n_patients     = 1L,
      cutoff         = 113L,
      patient_cal    = c(100L),
      t_measure      = c(2L),
      t_day_measure  = c(14L),
      measure_pos    = c(1L, 2L)
    ),
    list(
      name           = "multiple_measures_partial",
      n_patients     = 1L,
      cutoff         = 130L,
      patient_cal    = c(100L),
      t_measure      = c(1L, 2L, 4L, 6L),
      t_day_measure  = c(7L, 14L, 28L, 42L),
      measure_pos    = c(1L, 5L)
    ),
    list(
      name           = "two_patients",
      n_patients     = 2L,
      cutoff         = 125L,
      patient_cal    = c(100L, 110L),
      t_measure      = c(1L, 3L, 2L, 4L),
      t_day_measure  = c(7L, 21L, 14L, 28L),
      measure_pos    = c(1L, 3L, 5L)
    )
  )

  # Compute R expected values
  cv_expected  <- lapply(cv_cases,  function(cc)
    r_cutoff_visits(cc$cutoff, cc$patient_cal, cc$t_visits, cc$t_visits_day, cc$visit_pos))
  fcv_expected <- lapply(fcv_cases, function(cc)
    r_fine_cutoff_visits(cc$cutoff, cc$patient_cal, cc$t_measure, cc$t_day_measure, cc$measure_pos))

  # Build padded Stan data
  N_CV  <- length(cv_cases)
  N_FCV <- length(fcv_cases)

  MAX_P_CV <- max(sapply(cv_cases, function(x) x$n_patients))
  MAX_V_CV <- max(sapply(cv_cases, function(x) length(x$t_visits)))

  MAX_P_FCV <- max(sapply(fcv_cases, function(x) x$n_patients))
  MAX_M_FCV <- max(sapply(fcv_cases, function(x) length(x$t_measure)))

  pad <- function(x, n, fill = 0L) {
    if (length(x) >= n) return(x[seq_len(n)])
    c(x, rep(fill, n - length(x)))
  }

  cv_p_mat  <- matrix(0L, N_CV, MAX_P_CV)
  cv_v_mat  <- matrix(0L, N_CV, MAX_V_CV)
  cv_vd_mat <- matrix(0L, N_CV, MAX_V_CV)
  cv_pos_mat <- matrix(1L, N_CV, MAX_P_CV + 1)
  n_p_cv <- integer(N_CV); n_v_cv <- integer(N_CV); cut_cv <- integer(N_CV)

  for (i in seq_len(N_CV)) {
    cc <- cv_cases[[i]]
    np <- cc$n_patients; nv <- length(cc$t_visits)
    n_p_cv[i] <- np; n_v_cv[i] <- nv; cut_cv[i] <- cc$cutoff
    cv_p_mat[i,  seq_len(np)]   <- cc$patient_cal
    cv_v_mat[i,  seq_len(nv)]   <- cc$t_visits
    cv_vd_mat[i, seq_len(nv)]   <- cc$t_visits_day
    cv_pos_mat[i, seq_len(np+1)] <- cc$visit_pos
  }

  fcv_p_mat  <- matrix(0L, N_FCV, MAX_P_FCV)
  fcv_m_mat  <- matrix(0L, N_FCV, MAX_M_FCV)
  fcv_md_mat <- matrix(0L, N_FCV, MAX_M_FCV)
  fcv_pos_mat <- matrix(1L, N_FCV, MAX_P_FCV + 1)
  n_p_fcv <- integer(N_FCV); n_m_fcv <- integer(N_FCV); cut_fcv <- integer(N_FCV)

  for (i in seq_len(N_FCV)) {
    cc <- fcv_cases[[i]]
    np <- cc$n_patients; nm <- length(cc$t_measure)
    n_p_fcv[i] <- np; n_m_fcv[i] <- nm; cut_fcv[i] <- cc$cutoff
    fcv_p_mat[i,   seq_len(np)]   <- cc$patient_cal
    fcv_m_mat[i,   seq_len(nm)]   <- cc$t_measure
    fcv_md_mat[i,  seq_len(nm)]   <- cc$t_day_measure
    fcv_pos_mat[i, seq_len(np+1)] <- cc$measure_pos
  }

  stan_data <- list(
    N_CV = N_CV, MAX_P_CV = MAX_P_CV, MAX_V_CV = MAX_V_CV,
    n_patients_cv = n_p_cv, n_visits_cv = n_v_cv,
    cutoff_cal_day_cv = cut_cv,
    patient_cal_day_cv = cv_p_mat,
    t_visits_cv = cv_v_mat, t_visits_day_cv = cv_vd_mat,
    visit_pos_cv = cv_pos_mat,
    N_FCV = N_FCV, MAX_P_FCV = MAX_P_FCV, MAX_M_FCV = MAX_M_FCV,
    n_patients_fcv = n_p_fcv, n_measures_fcv = n_m_fcv,
    cutoff_cal_day_fcv = cut_fcv,
    patient_cal_day_fcv = fcv_p_mat,
    t_measure_fcv = fcv_m_mat, t_day_measure_fcv = fcv_md_mat,
    measure_pos_fcv = fcv_pos_mat
  )

  fit      <- test_stan_function(
    "tests/testthat/stan/test_lfo_cutoff_visits_all.stan",
    data = stan_data
  )
  draws_df <- posterior::as_draws_df(fit$draws())
  get_val  <- function(var, ...) get_stan_val(draws_df, var, ...)

  # Assert cutoff_visits
  for (i in seq_len(N_CV)) {
    cc  <- cv_cases[[i]]
    exp <- cv_expected[[i]]
    for (p in seq_len(cc$n_patients)) {
      lbl <- paste0("cv case ", i, " (", cc$name, ") patient ", p)
      expect_equal(get_val("cv_last_visit_day",  i, p), exp$last_visit_day[p],           label = paste(lbl, "last_visit_day"))
      expect_equal(get_val("cv_last_visit_week", i, p), exp$last_visit_week[p],          label = paste(lbl, "last_visit_week"))
      expect_equal(get_val("cv_last_cal_day",    i, p), exp$last_visit_calendar_day[p],  label = paste(lbl, "last_cal_day"))
      expect_equal(get_val("cv_week_observed",   i, p), exp$last_visit_week_observed[p], label = paste(lbl, "week_observed"))
      expect_equal(get_val("cv_cutoff_last_idx", i, p), exp$cutoff_last_visit_idx[p],    label = paste(lbl, "cutoff_last_idx"))
    }
  }

  # Assert fine_cutoff_visits
  for (i in seq_len(N_FCV)) {
    cc  <- fcv_cases[[i]]
    exp <- fcv_expected[[i]]
    for (p in seq_len(cc$n_patients)) {
      lbl <- paste0("fcv case ", i, " (", cc$name, ") patient ", p)
      expect_equal(get_val("fcv_last_visit_day",  i, p), exp$last_visit_day[p],          label = paste(lbl, "last_visit_day"))
      expect_equal(get_val("fcv_last_visit_week", i, p), exp$last_visit_week[p],         label = paste(lbl, "last_visit_week"))
      expect_equal(get_val("fcv_last_cal_day",    i, p), exp$last_visit_calendar_day[p], label = paste(lbl, "last_cal_day"))
    }
  }
})
```

- [ ] **Step 3: Run integration test**

```bash
Rscript -e 'testthat::test_file("tests/testthat/test-stan-lfo-cutoff-visits.R")'
```

Expected: PASS — 40 assertions (5 cases × 1 patient + 1 case × 2 patients = 7 patient-slots × 5 outputs = 35 for cv; 5 cases × ~1.2 avg patients × 3 outputs = 18 for fcv; total ~53 assertions).

- [ ] **Step 4: Run full test suite**

```bash
Rscript -e 'testthat::test_dir("tests/testthat")'
```

Expected: PASS all, count ≥ 811 (783 + ~28 new).

- [ ] **Step 5: Commit**

```bash
git add tests/testthat/stan/test_lfo_cutoff_visits_all.stan \
        tests/testthat/test-stan-lfo-cutoff-visits.R
git commit -m "test: add cutoff_visits + fine_cutoff_visits Stan tests — LFO coverage complete"
```

---

## Summary

| Task | New files | Assertions |
|---|---|---|
| 1. R references + unit tests | 1 modified, 1 created | ~28 |
| 2. Stan harness + integration | 2 created | ~53 |
| **Total** | **3 files** | **~81** |

LFO pipeline test coverage is now complete: `get_oos_patients_idx`, `get_testing_visit_week_bounds`, `cutoff_visits`, `fine_cutoff_visits` all tested.
