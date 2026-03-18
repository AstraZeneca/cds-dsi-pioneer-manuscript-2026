library(testthat)
source(here::here("testthat/helper-lfo.R"))

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
