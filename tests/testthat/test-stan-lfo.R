library(testthat)
library(cmdstanr)
source(here::here("tests/testthat/helper-stan.R"))

# Test cases for get_testing_visit_week_bounds function
get_testing_visit_week_bounds_cases <- list(
  # Case 1: Simple case with 2 patients, 2 cutoffs
  list(
    case_name = "simple_2_patients_2_cutoffs",
    n_patients = 2L,
    n_cutoffs = 2L,
    n_visits = 6L,
    oos_patient_idx = c(1L, 1L),
    last_visit_calendar_day_sort_idx = c(1L, 2L),
    cutoff_calendar_day = c(130L, 150L),
    patient_calendar_day = c(100L, 110L),
    t_patient_visits = c(1L, 3L, 5L, 2L, 4L, 6L),
    t_patient_visits_day = c(7L, 21L, 35L, 14L, 28L, 42L),
    patient_visit_pos = c(1L, 4L, 7L),
    expected_first_testing_visit_week = matrix(
      c(5L, 0L, 4L, 0L),
      nrow = 2,
      ncol = 2,
      byrow = TRUE
    ),
    expected_testing_start_idx = matrix(
      c(3L, 0L, 5L, 0L),
      nrow = 2,
      ncol = 2,
      byrow = TRUE
    )
  )
)


# Inline all test cases for cutoff_visits
# (already fixed above, no duplicate or stray lists)

# Inline all test cases for fine_cutoff_visits
# fine_cutoff_visits is like cutoff_visits but operates on tumor measurements
# (t_measure / t_day_measure / patient_tumor_measure_pos) instead of visit arrays.
# Expected values derived by: for each patient i, find last measure j where
#   patient_calendar_day[i] + t_day_measure[j] <= cutoff_calendar_day
fine_cutoff_visits_cases <- list(
  # Case 1: single patient, 3 measures all before cutoff
  list(
    cutoff_calendar_day = 30,
    patient_calendar_day = c(0),
    t_measure = c(1, 2, 3),
    t_day_measure = c(7, 14, 21),        # absolute days: 7, 14, 21 — all <= 30
    patient_tumor_measure_pos = c(1, 4),
    expected_last_visit_day = c(21),
    expected_last_visit_week = c(3)
  ),
  # Case 2: single patient, no measures before cutoff
  list(
    cutoff_calendar_day = 5,
    patient_calendar_day = c(0),
    t_measure = c(1, 2),
    t_day_measure = c(10, 20),           # absolute days: 10, 20 — none <= 5
    patient_tumor_measure_pos = c(1, 3),
    expected_last_visit_day = c(0),
    expected_last_visit_week = c(0)
  ),
  # Case 3: single patient, last measure exactly at cutoff boundary
  list(
    cutoff_calendar_day = 10,
    patient_calendar_day = c(0),
    t_measure = c(1, 2),
    t_day_measure = c(7, 10),            # day 10 exactly at cutoff — should be included
    patient_tumor_measure_pos = c(1, 3),
    expected_last_visit_day = c(10),
    expected_last_visit_week = c(2)
  ),
  # Case 4: two patients with different entry dates
  list(
    cutoff_calendar_day = 20,
    patient_calendar_day = c(0, 5),
    # Patient 1 (entry=0): relative cutoff=20, measures at days 7,14 — both qualify, last=14 wk2
    # Patient 2 (entry=5): relative cutoff=15, measures at days 8,18 — day 8 qualifies, last=8 wk1
    t_measure = c(1, 2, 1, 2),
    t_day_measure = c(7, 14, 8, 18),
    patient_tumor_measure_pos = c(1, 3, 5),
    expected_last_visit_day = c(14, 8),
    expected_last_visit_week = c(2, 1)
  )
)

# Inline all test cases for get_oos_patients_idx
get_oos_patients_idx_cases <- list(
  # Minimal dummy OOS test case
  list(
    n_patients = 1L,
    n_cutoffs = 1L,
    sorted_last_visit_calendar_day = matrix(0L, 1, 1),
    cutoff_calendar_day = matrix(0L, 1, 1)
    # No expected output needed for now
  )
)

# --- Prepare data for cutoff_visits ---
# test-stan-lfo.R
# All-in-one test for lfo.stan functions (cutoff_visits, etc.)
# Follows the same pattern as pos.stan tests

# Ensure test helpers (e.g., stan_test_helper) are available
source(here::here("tests/testthat/helper-stan.R"))

test_that("cutoff_visits and fine_cutoff_visits handle all edge cases and boundaries in a single Stan run", {
  # Load fine_cutoff_visits cases

  fine_cases <- fine_cutoff_visits_cases

  # Inline cutoff_visits cases
  cutoff_cases <- list(
    # 1. All patients have visits before the cutoff, different entry dates
    list(
      cutoff_calendar_day = 30,
      patient_calendar_day = c(0, 5, 10),
      t_patient_visits = c(1, 2, 3, 1, 2, 1, 2, 3, 4),
      t_patient_visits_day = c(5, 10, 15, 8, 18, 2, 4, 6, 8),
      patient_visit_pos = c(1, 4, 6, 10),
      expected_last_visit_day = c(15, 18, 8),
      expected_last_visit_week = c(3, 2, 4),
      expected_cutoff_last_visit_idx = c(3, 5, 9) # Patient 1: visits 1,2,3 -> last is idx 3; Patient 2: visits 4,5 -> last is idx 5; Patient 3: visits 6,7,8,9 -> last is idx 9
    ),
    # 2. Some patients have no visits before the cutoff
    list(
      cutoff_calendar_day = 5,
      patient_calendar_day = c(0, 0, 0),
      t_patient_visits = c(1, 2, 1, 2, 1, 2),
      t_patient_visits_day = c(10, 20, 10, 20, 10, 20),
      patient_visit_pos = c(1, 3, 5, 7),
      expected_last_visit_day = c(0, 0, 0),
      expected_last_visit_week = c(0, 0, 0),
      expected_cutoff_last_visit_idx = c(0, 0, 0) # No visits before cutoff
    ),
    # 3. Some patients have visits exactly at the cutoff
    list(
      cutoff_calendar_day = 10,
      patient_calendar_day = c(0, 5),
      t_patient_visits = c(1, 2, 1, 2),
      t_patient_visits_day = c(10, 20, 10, 20),
      patient_visit_pos = c(1, 3, 5),
      expected_last_visit_day = c(10, 0),
      expected_last_visit_week = c(1, 0),
      expected_cutoff_last_visit_idx = c(1, 0) # Patient 1: visit 1 has week 1 <= 1; Patient 2: no visits before cutoff
    ),
    # 4. Patient with visits at irregular intervals, cutoff between visits
    list(
      cutoff_calendar_day = 17,
      patient_calendar_day = c(0),
      t_patient_visits = c(1, 2, 3, 4),
      t_patient_visits_day = c(5, 10, 20, 30),
      patient_visit_pos = c(1, 5),
      expected_last_visit_day = c(10),
      expected_last_visit_week = c(2),
      expected_cutoff_last_visit_idx = c(2) # Visits 1,2 have weeks 1,2 <= 2
    ),
    # 5. Patient with negative entry date (enrolled before study "start")
    list(
      cutoff_calendar_day = 10,
      patient_calendar_day = c(-5),
      t_patient_visits = c(1, 2),
      t_patient_visits_day = c(3, 8),
      patient_visit_pos = c(1, 3),
      expected_last_visit_day = c(8),
      expected_last_visit_week = c(2),
      expected_cutoff_last_visit_idx = c(2) # Visits 1,2 have weeks 1,2 <= 2
    ),
    # 6. Patient with visits on non-monotonic days (should be sorted)
    list(
      cutoff_calendar_day = 15,
      patient_calendar_day = c(0),
      t_patient_visits = c(1, 2, 3),
      t_patient_visits_day = c(10, 5, 15),
      patient_visit_pos = c(1, 4),
      expected_last_visit_day = c(15),
      expected_last_visit_week = c(3),
      expected_cutoff_last_visit_idx = c(3) # Visits 1,2,3 have weeks 1,2,3 <= 3
    ),
    # 7. Patient with duplicate visit days
    list(
      cutoff_calendar_day = 10,
      patient_calendar_day = c(0),
      t_patient_visits = c(1, 2, 3),
      t_patient_visits_day = c(5, 5, 10),
      patient_visit_pos = c(1, 4),
      expected_last_visit_day = c(10),
      expected_last_visit_week = c(3),
      expected_cutoff_last_visit_idx = c(3) # Visits 1,2,3 have weeks 1,2,3 <= 3
    ),
    # 8. Patient with only one visit, after cutoff
    list(
      cutoff_calendar_day = 5,
      patient_calendar_day = c(0),
      t_patient_visits = c(1),
      t_patient_visits_day = c(10),
      patient_visit_pos = c(1, 2),
      expected_last_visit_day = c(0),
      expected_last_visit_week = c(0),
      expected_cutoff_last_visit_idx = c(0) # No visits before cutoff
    )
  )

  # Combine all cases for both functions
  all_cases <- c(
    lapply(fine_cases, function(x) c(x, list(.type = "fine"))),
    lapply(cutoff_cases, function(x) c(x, list(.type = "cutoff")))
  )
  N_cases <- length(all_cases)
  max_n_patients <- max(sapply(all_cases, function(x) {
    length(x$patient_calendar_day)
  }))
  max_n_measures <- max(
    1,
    max(sapply(all_cases, function(x) {
      if (!is.null(x$t_measure)) length(x$t_measure) else 0
    }))
  )
  max_n_visits <- max(sapply(all_cases, function(x) {
    if (!is.null(x$t_patient_visits)) length(x$t_patient_visits) else 0
  }))

  # Rectangularize all arrays
  cutoff_calendar_day <- integer(N_cases)
  patient_calendar_day <- matrix(0, N_cases, max_n_patients)
  t_measure <- matrix(0, N_cases, max_n_measures)
  t_day_measure <- matrix(0, N_cases, max_n_measures)
  patient_tumor_measure_pos <- matrix(1, N_cases, max_n_patients + 1)
  n_measures <- rep(0L, N_cases)
  t_patient_visits <- matrix(0, N_cases, max_n_visits)
  t_patient_visits_day <- matrix(0, N_cases, max_n_visits)
  patient_visit_pos <- matrix(1, N_cases, max_n_patients + 1)
  n_patients <- integer(N_cases)
  n_visits <- integer(N_cases)

  # Expected outputs for both functions
  expected_fine_last_visit_day <- matrix(NA_integer_, N_cases, max_n_patients)
  expected_fine_last_visit_week <- matrix(NA_integer_, N_cases, max_n_patients)
  expected_cutoff_last_visit_day <- matrix(NA_integer_, N_cases, max_n_patients)
  expected_cutoff_last_visit_week <- matrix(
    NA_integer_,
    N_cases,
    max_n_patients
  )
  expected_cutoff_last_visit_idx <- matrix(NA_integer_, N_cases, max_n_patients)

  for (i in seq_len(N_cases)) {
    case <- all_cases[[i]]
    np <- length(case$patient_calendar_day)
    cutoff_calendar_day[i] <- case$cutoff_calendar_day
    patient_calendar_day[i, seq_len(np)] <- case$patient_calendar_day
    n_patients[i] <- np

    if (!is.null(case$t_measure) && !is.null(case$t_day_measure)) {
      # fine_cutoff_visits case
      nm <- length(case$t_measure)
      if (nm > 0) {
        t_measure[i, seq_len(nm)] <- case$t_measure
        t_day_measure[i, seq_len(nm)] <- case$t_day_measure
      }
      # pad with zeros if nm == 0
      if (nm == 0) {
        t_measure[i, 1] <- 0
        t_day_measure[i, 1] <- 0
      }
      patient_tumor_measure_pos[
        i,
        seq_len(np + 1)
      ] <- case$patient_tumor_measure_pos
      n_measures[i] <- nm
      expected_fine_last_visit_day[i, 1:np] <- case$expected_last_visit_day
      expected_fine_last_visit_week[i, 1:np] <- case$expected_last_visit_week
    } else {
      # cutoff_visits-only case
      nv <- length(case$t_patient_visits)
      if (nv > 0) {
        t_patient_visits[i, seq_len(nv)] <- case$t_patient_visits
        t_patient_visits_day[i, seq_len(nv)] <- case$t_patient_visits_day
      }
      patient_visit_pos[i, seq_len(np + 1)] <- case$patient_visit_pos
      n_visits[i] <- nv
      # For fine_cutoff_visits, force one dummy measure for these patients
      n_measures[i] <- 0
      t_measure[i, 1] <- 0
      t_day_measure[i, 1] <- 0
      patient_tumor_measure_pos[i, seq_len(np + 1)] <- rep(1L, np + 1)
      expected_cutoff_last_visit_day[i, 1:np] <- case$expected_last_visit_day
      expected_cutoff_last_visit_week[i, 1:np] <- case$expected_last_visit_week
      expected_cutoff_last_visit_idx[
        i,
        1:np
      ] <- case$expected_cutoff_last_visit_idx
      expected_fine_last_visit_day[i, 1:np] <- 0
      expected_fine_last_visit_week[i, 1:np] <- 0
    }
  }

  # Prepare OOS (get_oos_patients_idx) data from dummy case
  # Prepare testing bounds cases data
  n_testing_bounds_cases <- length(get_testing_visit_week_bounds_cases)
  max_n_patients_testing <- if (n_testing_bounds_cases > 0) {
    max(sapply(get_testing_visit_week_bounds_cases, function(x) x$n_patients))
  } else {
    1L
  }
  max_n_cutoffs_testing <- if (n_testing_bounds_cases > 0) {
    max(sapply(get_testing_visit_week_bounds_cases, function(x) x$n_cutoffs))
  } else {
    1L
  }
  max_n_visits_testing <- if (n_testing_bounds_cases > 0) {
    max(sapply(get_testing_visit_week_bounds_cases, function(x) x$n_visits))
  } else {
    1L
  }

  # Initialize testing bounds arrays
  n_patients_testing <- integer(max(n_testing_bounds_cases, 1))
  n_cutoffs_testing <- integer(max(n_testing_bounds_cases, 1))
  n_visits_testing <- integer(max(n_testing_bounds_cases, 1))
  oos_patient_idx_testing <- matrix(
    1L,
    max(n_testing_bounds_cases, 1),
    max_n_cutoffs_testing
  )
  last_visit_calendar_day_sort_idx_testing <- matrix(
    1L,
    max(n_testing_bounds_cases, 1),
    max_n_patients_testing
  )
  cutoff_calendar_day_testing <- matrix(
    0L,
    max(n_testing_bounds_cases, 1),
    max_n_cutoffs_testing
  )
  patient_calendar_day_testing <- matrix(
    0L,
    max(n_testing_bounds_cases, 1),
    max_n_patients_testing
  )
  t_patient_visits_testing <- matrix(
    0L,
    max(n_testing_bounds_cases, 1),
    max_n_visits_testing
  )
  t_patient_visits_day_testing <- matrix(
    0L,
    max(n_testing_bounds_cases, 1),
    max_n_visits_testing
  )
  patient_visit_pos_testing <- matrix(
    1L,
    max(n_testing_bounds_cases, 1),
    max_n_patients_testing + 1
  )

  # Fill testing bounds data
  if (n_testing_bounds_cases > 0) {
    for (i in seq_len(n_testing_bounds_cases)) {
      case <- get_testing_visit_week_bounds_cases[[i]]
      n_patients_testing[i] <- case$n_patients
      n_cutoffs_testing[i] <- case$n_cutoffs
      n_visits_testing[i] <- case$n_visits
      oos_patient_idx_testing[i, 1:case$n_cutoffs] <- case$oos_patient_idx
      last_visit_calendar_day_sort_idx_testing[
        i,
        1:case$n_patients
      ] <- case$last_visit_calendar_day_sort_idx
      cutoff_calendar_day_testing[
        i,
        1:case$n_cutoffs
      ] <- case$cutoff_calendar_day
      patient_calendar_day_testing[
        i,
        1:case$n_patients
      ] <- case$patient_calendar_day
      t_patient_visits_testing[i, 1:case$n_visits] <- case$t_patient_visits
      t_patient_visits_day_testing[
        i,
        1:case$n_visits
      ] <- case$t_patient_visits_day
      patient_visit_pos_testing[
        i,
        1:(case$n_patients + 1)
      ] <- case$patient_visit_pos
    }
  }
  stan_data <- list(
    N_cases = N_cases,
    max_n_patients = max_n_patients,
    max_n_measures = max_n_measures,
    max_n_visits = max_n_visits,
    cutoff_calendar_day = cutoff_calendar_day,
    patient_calendar_day = patient_calendar_day,
    t_measure = t_measure,
    t_day_measure = t_day_measure,
    patient_tumor_measure_pos = patient_tumor_measure_pos,
    n_measures = n_measures,
    t_patient_visits = t_patient_visits,
    t_patient_visits_day = t_patient_visits_day,
    patient_visit_pos = patient_visit_pos,
    n_patients = n_patients,
    n_visits = n_visits,
    # Minimal valid OOS data
    n_oos_cases = 1L,
    max_n_patients_oos = 1L,
    max_n_cutoffs_oos = 1L,
    n_patients_oos = 1L,
    n_cutoffs_oos = 1L,
    sorted_last_visit_calendar_day_oos = matrix(0L, 1, 1),
    cutoff_calendar_day_oos = matrix(0L, 1, 1),
    # Testing bounds data
    n_testing_bounds_cases = n_testing_bounds_cases,
    max_n_patients_testing = max_n_patients_testing,
    max_n_cutoffs_testing = max_n_cutoffs_testing,
    max_n_visits_testing = max_n_visits_testing,
    n_patients_testing = n_patients_testing,
    n_cutoffs_testing = n_cutoffs_testing,
    n_visits_testing = n_visits_testing,
    oos_patient_idx_testing = oos_patient_idx_testing,
    last_visit_calendar_day_sort_idx_testing = last_visit_calendar_day_sort_idx_testing,
    cutoff_calendar_day_testing = cutoff_calendar_day_testing,
    patient_calendar_day_testing = patient_calendar_day_testing,
    t_patient_visits_testing = t_patient_visits_testing,
    t_patient_visits_day_testing = t_patient_visits_day_testing,
    patient_visit_pos_testing = patient_visit_pos_testing
  )

  fit <- test_stan_function("stan/test_lfo_all.stan", data = stan_data)
  library(tidybayes)
  draws_df <- spread_draws(
    fit$draws(),
    last_visit_day[case, patient],
    last_visit_week[case, patient],
    fine_last_visit_day[case, patient],
    fine_last_visit_week[case, patient],
    cutoff_last_visit_idx[case, patient]
  )

  for (i in seq_len(N_cases)) {
    np <- n_patients[i]
    df <- draws_df[
      draws_df$.iteration == 1 &
        draws_df$.chain == 1 &
        draws_df$case == i &
        draws_df$patient <= np,
    ]
    # Check cutoff_visits outputs if expected values are not NA
    if (!all(is.na(expected_cutoff_last_visit_day[i, 1:np]))) {
      expect_equal(
        as.integer(df$last_visit_day),
        expected_cutoff_last_visit_day[i, 1:np],
        label = paste("last_visit_day, case", i)
      )
      expect_equal(
        as.integer(df$last_visit_week),
        expected_cutoff_last_visit_week[i, 1:np],
        label = paste("last_visit_week, case", i)
      )
      expect_equal(
        as.integer(df$cutoff_last_visit_idx),
        expected_cutoff_last_visit_idx[i, 1:np],
        label = paste("cutoff_last_visit_idx, case", i)
      )
    }
    # Check fine_cutoff_visits outputs if expected values are not NA
    if (!all(is.na(expected_fine_last_visit_day[i, 1:np]))) {
      expect_equal(
        as.integer(df$fine_last_visit_day),
        expected_fine_last_visit_day[i, 1:np],
        label = paste("fine_last_visit_day, case", i)
      )
      expect_equal(
        as.integer(df$fine_last_visit_week),
        expected_fine_last_visit_week[i, 1:np],
        label = paste("fine_last_visit_week, case", i)
      )
    }
  }

  # get_testing_visit_week_bounds is exercised by the Stan run above but not
  # asserted here — see test-stan-get_testing_visit_week_bounds.R for dedicated
  # coverage of that function.
})
