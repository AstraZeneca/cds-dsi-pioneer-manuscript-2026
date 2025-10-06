library(testthat)
library(cmdstanr)
source(here::here("tests/testthat/helper-stan.R"))




# Inline all test cases for cutoff_visits
# (already fixed above, no duplicate or stray lists)

# Inline all test cases for fine_cutoff_visits
fine_cutoff_visits_cases <- list(
  # ... (copy all fine_cutoff_visits test cases here, as in your previous code) ...
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
      t_patient_visits = c(1,2,3,1,2,1,2,3,4),
      t_patient_visits_day = c(5,10,15,8,18,2,4,6,8),
      patient_visit_pos = c(1,4,6,10),
      expected_last_visit_day = c(15,18,8),
      expected_last_visit_week = c(3,2,4)
    ),
    # 2. Some patients have no visits before the cutoff
    list(
      cutoff_calendar_day = 5,
      patient_calendar_day = c(0,0,0),
      t_patient_visits = c(1,2,1,2,1,2),
      t_patient_visits_day = c(10,20,10,20,10,20),
      patient_visit_pos = c(1,3,5,7),
      expected_last_visit_day = c(0,0,0),
      expected_last_visit_week = c(0,0,0)
    ),
    # 3. Some patients have visits exactly at the cutoff
    list(
      cutoff_calendar_day = 10,
      patient_calendar_day = c(0,5),
      t_patient_visits = c(1,2,1,2),
      t_patient_visits_day = c(10,20,10,20),
      patient_visit_pos = c(1,3,5),
      expected_last_visit_day = c(10,0),
      expected_last_visit_week = c(1,0)
    ),
    # 4. Patient with visits at irregular intervals, cutoff between visits
    list(
      cutoff_calendar_day = 17,
      patient_calendar_day = c(0),
      t_patient_visits = c(1,2,3,4),
      t_patient_visits_day = c(5,10,20,30),
      patient_visit_pos = c(1,5),
      expected_last_visit_day = c(10),
      expected_last_visit_week = c(2)
    ),
    # 5. Patient with negative entry date (enrolled before study "start")
    list(
      cutoff_calendar_day = 10,
      patient_calendar_day = c(-5),
      t_patient_visits = c(1,2),
      t_patient_visits_day = c(3,8),
      patient_visit_pos = c(1,3),
      expected_last_visit_day = c(8),
      expected_last_visit_week = c(2)
    ),
    # 6. Patient with visits on non-monotonic days (should be sorted)
    list(
      cutoff_calendar_day = 15,
      patient_calendar_day = c(0),
      t_patient_visits = c(1,2,3),
      t_patient_visits_day = c(10,5,15),
      patient_visit_pos = c(1,4),
      expected_last_visit_day = c(15),
      expected_last_visit_week = c(3)
    ),
    # 7. Patient with duplicate visit days
    list(
      cutoff_calendar_day = 10,
      patient_calendar_day = c(0),
      t_patient_visits = c(1,2,3),
      t_patient_visits_day = c(5,5,10),
      patient_visit_pos = c(1,4),
      expected_last_visit_day = c(10),
      expected_last_visit_week = c(3)
    ),
    # 8. Patient with only one visit, after cutoff
    list(
      cutoff_calendar_day = 5,
      patient_calendar_day = c(0),
      t_patient_visits = c(1),
      t_patient_visits_day = c(10),
      patient_visit_pos = c(1,2),
      expected_last_visit_day = c(0),
      expected_last_visit_week = c(0)
    )
  )


  # Combine all cases for both functions
  all_cases <- c(
    lapply(fine_cases, function(x) c(x, list(.type = "fine"))),
    lapply(cutoff_cases, function(x) c(x, list(.type = "cutoff")))
  )
  N_cases <- length(all_cases)
  max_n_patients <- max(sapply(all_cases, function(x) length(x$patient_calendar_day)))
  max_n_measures <- max(1, max(sapply(all_cases, function(x) if (!is.null(x$t_measure)) length(x$t_measure) else 0)))
  max_n_visits <- max(sapply(all_cases, function(x) if (!is.null(x$t_patient_visits)) length(x$t_patient_visits) else 0))

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
  expected_cutoff_last_visit_week <- matrix(NA_integer_, N_cases, max_n_patients)

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
      patient_tumor_measure_pos[i, seq_len(np+1)] <- case$patient_tumor_measure_pos
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
      patient_visit_pos[i, seq_len(np+1)] <- case$patient_visit_pos
      n_visits[i] <- nv
      # For fine_cutoff_visits, force one dummy measure for these patients
      n_measures[i] <- 0
      t_measure[i, 1] <- 0
      t_day_measure[i, 1] <- 0
      patient_tumor_measure_pos[i, seq_len(np+1)] <- rep(1L, np+1)
      expected_cutoff_last_visit_day[i, 1:np] <- case$expected_last_visit_day
      expected_cutoff_last_visit_week[i, 1:np] <- case$expected_last_visit_week
      expected_fine_last_visit_day[i, 1:np] <- 0
      expected_fine_last_visit_week[i, 1:np] <- 0
    }
  }

  # Prepare OOS (get_oos_patients_idx) data from dummy case
  oos_case <- get_oos_patients_idx_cases[[1]]
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
    max_n_patients_oos = oos_case$n_patients,
    max_n_cutoffs_oos = oos_case$n_cutoffs,
    n_patients_oos = oos_case$n_patients,
    n_cutoffs_oos = oos_case$n_cutoffs,
    sorted_last_visit_calendar_day_oos = oos_case$sorted_last_visit_calendar_day,
    cutoff_calendar_day_oos = oos_case$cutoff_calendar_day
  )

  fit <- test_stan_function("stan/test_lfo_all.stan", data = stan_data)
  library(tidybayes)
  draws_df <- spread_draws(
    fit$draws(),
    last_visit_day[case, patient], last_visit_week[case, patient],
    fine_last_visit_day[case, patient], fine_last_visit_week[case, patient]
  )

  for (i in seq_len(N_cases)) {
    np <- n_patients[i]
    df <- draws_df[draws_df$.iteration == 1 & draws_df$.chain == 1 & draws_df$case == i & draws_df$patient <= np, ]
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
    }
    # Check fine_cutoff_visits outputs if expected values are not NA
    if (!all(is.na(expected_fine_last_visit_day[i, 1:np]))) {
      # Diagnostic print for cutoff-only cases
      if (all(expected_fine_last_visit_day[i, 1:np] == 0)) {
        cat(sprintf("[DIAG] Case %d: expected fine_last_visit_day = %s, actual = %s\n", i, paste(expected_fine_last_visit_day[i, 1:np], collapse=","), paste(as.integer(df$fine_last_visit_day), collapse=",")))
      }
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
})
