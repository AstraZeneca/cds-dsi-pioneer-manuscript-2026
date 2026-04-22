library(testthat)
library(dplyr)
library(here)

source(here("r", "pioneer", "prepare_analysis_data.R"))

# =============================================================================
# Test helpers
# =============================================================================

# Minimal raw trial patient row for prepare_trial_analysis_data() tests.
# Provides sensible defaults; override specific fields per test.
make_raw_trial_patient <- function(
  usubjid                    = "PT-001",
  studyid                    = "FPI-TEST",
  death                      = FALSE,
  death_week                 = NA_integer_,
  # Overall PFS (radiographic composite — NOT PSA-specific)
  pfs                        = 20L,
  right_censored             = TRUE,
  interval_censored          = 0L,
  progress_week              = NA_integer_,
  progressor                 = FALSE,
  # PIONEER PSA-specific columns (what the fixes should use)
  pioneer_progressor         = 0L,
  deterministic_PFS_PIONEER  = 20L,
  patient_max_t              = 20L,
  calendar_week              = 100L,
  calendar_day               = 700L
) {
  tibble(
    studyid                   = studyid,
    usubjid                   = usubjid,
    death                     = death,
    death_week                = as.integer(death_week),
    pfs                       = as.integer(pfs),
    right_censored            = right_censored,
    interval_censored         = as.integer(interval_censored),
    progress_week             = as.integer(progress_week),
    progressor                = progressor,
    pioneer_progressor        = pioneer_progressor,
    deterministic_PFS_PIONEER = as.integer(deterministic_PFS_PIONEER),
    patient_max_t             = as.integer(patient_max_t),
    calendar_week             = as.integer(calendar_week),
    calendar_day              = as.integer(calendar_day),
    # Required by validate_flatiron_analysis_data() but not under test here
    trtsdt = as.Date("2024-01-01"),
    age = 65L, sex = "M", country = "USA", ecogbl = 1L, bmibl = 25.0,
    race = "White", arm = "Pluvicto monotherapy",
    visceral_mets = FALSE, lymph_mets = FALSE, bone_mets = FALSE,
    liver_mets = FALSE,
    baseline_albumin = 4.0, baseline_ALP = 100.0, baseline_AST = 30.0,
    baseline_chloride = 100.0, baseline_creatinine = 1.0,
    baseline_hematocrit = 40.0, baseline_hemoglobin = 13.0,
    baseline_ldh = 200.0, baseline_lymph = 1.5, baseline_monocytes = 0.5,
    baseline_neutrophils = 5.0, baseline_nlr = 3.0,
    prior_taxanes = FALSE, imputed = FALSE,
    pfs_cnsr = as.integer(right_censored),
    progression_before_death = NA,
    pioneer_pfs_cnsr = 1L,
    pd_reason = NA_character_,
    os_time = as.integer(patient_max_t),
    patient_min_t = 0L, patient_first_visit = 1L, patient_last_visit = patient_max_t,
    patient_t_width = patient_max_t, trtedt = as.Date("2024-06-01"),
    treatment_end_week = patient_max_t, treatment_end_day = patient_max_t * 7L,
    age_group = "65+", pioneer_progress_week = NA_integer_
  )
}

make_raw_trial_visit <- function(studyid = "FPI-TEST", usubjid = "PT-001",
                                 week = 4L, measurement_value = 10.0,
                                 psa_pd_confirmed = 0L) {
  tibble(
    studyid          = studyid,
    usubjid          = usubjid,
    visitnum         = 1L,
    ady              = as.integer(week * 7),
    day              = as.integer(week * 7),
    week             = as.integer(week),
    measurement_value = measurement_value,
    response         = NA_character_,
    psa_nadir = measurement_value, psa_change_from_nadir = 0.0,
    psa_pct_change_from_nadir = 0.0,
    psa_pd = 0L, psa_pd_confirmed = as.integer(psa_pd_confirmed),
    psa50 = 0L, psa90 = 0L, bone_pd = 0L, bone_pd_confirmed = 0L,
    new_lesion = 0L, combined_response = "SD", combined_pd = 0L,
    radiographic_pd = 0L
  )
}

# Build a multi-visit tibble for a single patient from parallel vectors
make_visits <- function(usubjid = "PT-001", weeks, psas,
                        psa_pd_confirmed = rep(0L, length(weeks))) {
  map2(weeks, psas, \(w, p) make_raw_trial_visit(
    usubjid = usubjid, week = w, measurement_value = p,
    psa_pd_confirmed = psa_pd_confirmed[match(w, weeks)]
  )) |>
    bind_rows()
}

harmonize <- function(patient_rows, visit_rows = NULL) {
  if (is.null(visit_rows)) {
    visit_rows <- bind_rows(
      lapply(seq_len(nrow(patient_rows)), function(i) {
        make_raw_trial_visit(
          studyid = patient_rows$studyid[[i]],
          usubjid = patient_rows$usubjid[[i]]
        )
      })
    )
  }
  prepare_trial_analysis_data(patient_rows, visit_rows)
}

# =============================================================================
# Bug 1: psa_right_censored must come from pioneer_progressor, not right_censored
# =============================================================================

test_that("dead patient without PSA progression → psa_right_censored = TRUE", {
  # right_censored=FALSE (had an event=death) but pioneer_progressor=0 → no PSA PD
  # Bug: current code maps right_censored → psa_right_censored, giving FALSE (wrong)
  # Fix: map pioneer_progressor → psa_right_censored, giving TRUE (correct)
  pd <- make_raw_trial_patient(
    death              = TRUE,
    death_week         = 30L,
    right_censored     = FALSE,
    pioneer_progressor = 0L,
    pfs                = 29L,
    patient_max_t      = 20L
  )
  result <- harmonize(pd)
  expect_true(result$psa_right_censored[[1]])
})

test_that("patient with PSA progression → psa_right_censored = FALSE", {
  pd <- make_raw_trial_patient(
    pioneer_progressor        = 1L,
    deterministic_PFS_PIONEER = 15L,
    right_censored            = FALSE,  # had an event
    pfs                       = 15L,
    patient_max_t             = 20L
  )
  result <- harmonize(pd)
  expect_false(result$psa_right_censored[[1]])
})

test_that("right-censored patient with no events → psa_right_censored = TRUE", {
  pd <- make_raw_trial_patient(
    death              = FALSE,
    right_censored     = TRUE,
    pioneer_progressor = 0L
  )
  result <- harmonize(pd)
  expect_true(result$psa_right_censored[[1]])
})

# =============================================================================
# Bug 2: psa_pfs must come from deterministic_PFS_PIONEER, not pfs
# =============================================================================

test_that("PSA progressor: psa_pfs = deterministic_PFS_PIONEER, not pfs", {
  # Bug: pfs = 52 (death_week - 1); deterministic_PFS_PIONEER = 15 (actual PSA PD)
  # Current bug gives psa_pfs = 52; fix gives psa_pfs = 15
  pd <- make_raw_trial_patient(
    pioneer_progressor        = 1L,
    deterministic_PFS_PIONEER = 15L,
    pfs                       = 52L,   # composite PFS (different from PSA PD time)
    patient_max_t             = 52L
  )
  result <- harmonize(pd)
  expect_equal(result$psa_pfs[[1]], 15L)
})

test_that("non-progressor: psa_pfs = deterministic_PFS_PIONEER (= patient_max_t)", {
  pd <- make_raw_trial_patient(
    pioneer_progressor        = 0L,
    deterministic_PFS_PIONEER = 20L,
    pfs                       = 20L,
    patient_max_t             = 20L
  )
  result <- harmonize(pd)
  expect_equal(result$psa_pfs[[1]], 20L)
})

test_that("dead patient without PSA prog: psa_pfs = deterministic_PFS_PIONEER, not death_week-1", {
  # Bug: pfs = death_week - 1 = 29; fix: deterministic_PFS_PIONEER = 20 (last clean week)
  pd <- make_raw_trial_patient(
    death                     = TRUE,
    death_week                = 30L,
    right_censored            = FALSE,
    pioneer_progressor        = 0L,
    pfs                       = 29L,
    deterministic_PFS_PIONEER = 20L,
    patient_max_t             = 20L
  )
  result <- harmonize(pd)
  expect_equal(result$psa_pfs[[1]], 20L)
})

# =============================================================================
# Bug 3: other_events_pfs must capture RECIST progression, not just death
# =============================================================================

test_that("RECIST progressor alive: other_events_pfs = progress_week", {
  # Bug: current code only uses death → other_events_pfs = patient_max_t (wrong)
  # Fix: progressor=TRUE → other_events_pfs = progress_week
  pd <- make_raw_trial_patient(
    progressor    = TRUE,
    progress_week = 18L,
    right_censored = FALSE,
    death         = FALSE,
    patient_max_t = 30L,
    pfs           = 18L
  )
  result <- harmonize(pd)
  expect_equal(result$other_events_pfs[[1]], 18L)
  expect_false(result$other_events_right_censored[[1]])
})

test_that("RECIST progressor who then died: other_events_pfs = progress_week", {
  pd <- make_raw_trial_patient(
    progressor    = TRUE,
    progress_week = 18L,
    right_censored = FALSE,
    death         = TRUE,
    death_week    = 30L,
    patient_max_t = 30L,
    pfs           = 18L
  )
  result <- harmonize(pd)
  expect_equal(result$other_events_pfs[[1]], 18L)
  expect_false(result$other_events_right_censored[[1]])
})

test_that("dead patient without RECIST progression: other_events_pfs = death_week", {
  pd <- make_raw_trial_patient(
    progressor     = FALSE,
    death          = TRUE,
    death_week     = 30L,
    right_censored = FALSE,
    patient_max_t  = 20L,
    pfs            = 29L
  )
  result <- harmonize(pd)
  expect_equal(result$other_events_pfs[[1]], 30L)
  expect_false(result$other_events_right_censored[[1]])
})

test_that("right-censored patient without events: other_events_pfs = patient_max_t", {
  pd <- make_raw_trial_patient(
    progressor     = FALSE,
    death          = FALSE,
    right_censored = TRUE,
    patient_max_t  = 20L
  )
  result <- harmonize(pd)
  expect_equal(result$other_events_pfs[[1]], 20L)
  expect_true(result$other_events_right_censored[[1]])
})

# =============================================================================
# compute_pcwg3_from_psa: unit tests
# =============================================================================

test_that("PSA50 achiever gets PR", {
  # baseline (week -1) = 100, post PSA = 40 → 60% reduction → PR
  week <- c(-1L, 4L, 8L)
  psa  <- c(100, 40, 35)
  res  <- compute_pcwg3_from_psa(week, psa, rep(0L, 3))
  expect_equal(res[week > 0], c("PR", "PR"))
  expect_true(is.na(res[week < 0]))
})

test_that("undetectable PSA gets CR", {
  week <- c(-1L, 4L)
  psa  <- c(50, 0.05)
  res  <- compute_pcwg3_from_psa(week, psa, c(0L, 0L))
  expect_equal(res[[2]], "CR")
})

test_that("PSA below 50% reduction gets SD", {
  week <- c(-1L, 4L)
  psa  <- c(100, 60)   # only 40% reduction
  res  <- compute_pcwg3_from_psa(week, psa, c(0L, 0L))
  expect_equal(res[[2]], "SD")
})

test_that("confirmed PSA-PD gets PD", {
  week <- c(-1L, 4L)
  psa  <- c(100, 60)
  res  <- compute_pcwg3_from_psa(week, psa, c(0L, 1L))
  expect_equal(res[[2]], "PD")
})

test_that("missing post-baseline PSA gets NE", {
  week <- c(-1L, 4L)
  psa  <- c(100, NA_real_)
  res  <- compute_pcwg3_from_psa(week, psa, c(0L, 0L))
  expect_equal(res[[2]], "NE")
})

test_that("screening visits return NA", {
  week <- c(-4L, -1L, 4L)
  psa  <- c(120, 100, 40)
  res  <- compute_pcwg3_from_psa(week, psa, rep(0L, 3))
  expect_true(all(is.na(res[week <= 0])))
})

test_that("week 0 is excluded from baseline (regression: 115-203 pattern)", {
  # week -2 PSA = 53.3 (last pre-treatment), week 0 PSA = 38.1 (day-of-treatment)
  # Post PSA = 20.2 → PSA50 only when baseline = 53.3 (53.3/2 = 26.65 > 20.2)
  # Old bug: used week <= 0 → baseline = 38.1 → 38.1/2 = 19.05 < 20.2 → SD (wrong)
  week <- c(-2L, 0L, 11L)
  psa  <- c(53.3, 38.1, 20.2)
  res  <- compute_pcwg3_from_psa(week, psa, c(0L, 0L, 0L))
  expect_equal(res[[3]], "PR")   # must be PR, not SD
})

test_that("no pre-treatment PSA (only week 0) returns all NA", {
  # week 0 is day-of-dose; with week < 0 baseline, no valid baseline → all NA
  week <- c(0L, 4L)
  psa  <- c(100, 40)
  res  <- compute_pcwg3_from_psa(week, psa, c(0L, 0L))
  expect_true(all(is.na(res)))
})

# =============================================================================
# compute_pcwg3_from_psa: integration via prepare_trial_analysis_data
# =============================================================================

test_that("det_response in visit_data is derived from raw PSA, not response column", {
  # response column is NA at all visits (as in real data) — det_response must
  # still correctly classify PSA50 from measurement_value
  pd <- make_raw_trial_patient()
  visits <- make_visits(
    weeks = c(-1L, 4L, 8L),
    psas  = c(100, 40, 35)   # 60% and 65% reductions → PR
  )
  result <- harmonize(pd, visits)
  vd <- result$visit_data[[1]]
  expect_equal(vd$det_response[vd$week > 0], c("PR", "PR"))
})
