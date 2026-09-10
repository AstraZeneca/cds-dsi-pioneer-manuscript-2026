library(testthat)
library(cmdstanr)
library(posterior)
source(here::here("tests/testthat/helper-stan.R"))
source(here::here("tests/testthat/helper-lfo.R"))

test_that("get_testing_visit_week_bounds - multi-patient comprehensive test", {
  # Multi-patient test cases to validate the function handles multiple patients correctly
  test_cases <- list(
    # Case 1: Two patients, single cutoff - basic multi-patient scenario
    list(
      case_name = "two_patients_single_cutoff",
      n_patients = 2L,
      n_cutoffs = 1L,
      n_visits = 6L, # Total visits across both patients
      oos_patient_idx = c(1L), # Include all patients starting from patient 1
      last_visit_calendar_day_sort_idx = c(1L, 2L), # Both patients ordered by ID
      cutoff_calendar_day = c(115L), # Single cutoff at day 115
      patient_calendar_day = c(100L, 105L), # Patient 1 starts day 100, Patient 2 starts day 105
      # Patient 1: visits at weeks 1,2,3 (days 7,14,21 -> absolute 107,114,121)
      # Patient 2: visits at weeks 1,2,4 (days 7,14,28 -> absolute 112,119,133)
      t_patient_visits = c(1L, 2L, 3L, 1L, 2L, 4L),
      t_patient_visits_day = c(7L, 14L, 21L, 7L, 14L, 28L),
      patient_visit_pos = c(1L, 4L, 7L), # Patient 1: visits 1-3, Patient 2: visits 4-6
      should_error = FALSE
    ),

    # Case 2: Two patients, multiple cutoffs - complex scenario
    list(
      case_name = "two_patients_multiple_cutoffs",
      n_patients = 2L,
      n_cutoffs = 2L,
      n_visits = 8L,
      oos_patient_idx = c(1L, 1L), # Both cutoffs include all patients
      last_visit_calendar_day_sort_idx = c(1L, 2L),
      cutoff_calendar_day = c(115L, 125L), # Two cutoffs at days 115, 125
      patient_calendar_day = c(100L, 105L),
      # Patient 1: visits at weeks 1,2,3,5 (days 7,14,21,35 -> absolute 107,114,121,135)
      # Patient 2: visits at weeks 1,2,3,4 (days 7,14,21,28 -> absolute 112,119,126,133)
      t_patient_visits = c(1L, 2L, 3L, 5L, 1L, 2L, 3L, 4L),
      t_patient_visits_day = c(7L, 14L, 21L, 35L, 7L, 14L, 21L, 28L),
      patient_visit_pos = c(1L, 5L, 9L), # Patient 1: visits 1-4, Patient 2: visits 5-8
      should_error = FALSE
    ),

    # Case 3: Three patients, single cutoff - stress test with more patients
    list(
      case_name = "three_patients_single_cutoff",
      n_patients = 3L,
      n_cutoffs = 1L,
      n_visits = 9L,
      oos_patient_idx = c(1L),
      last_visit_calendar_day_sort_idx = c(1L, 2L, 3L),
      cutoff_calendar_day = c(120L),
      patient_calendar_day = c(100L, 110L, 105L), # Different start days
      # Patient 1: weeks 2,4 (days 14,28 -> absolute 114,128)
      # Patient 2: weeks 1,3,5 (days 7,21,35 -> absolute 117,131,145)
      # Patient 3: weeks 1,2,3,6 (days 7,14,21,42 -> absolute 112,119,126,147)
      t_patient_visits = c(2L, 4L, 1L, 3L, 5L, 1L, 2L, 3L, 6L),
      t_patient_visits_day = c(14L, 28L, 7L, 21L, 35L, 7L, 14L, 21L, 42L),
      patient_visit_pos = c(1L, 3L, 6L, 10L), # Patient 1: 1-2, Patient 2: 3-5, Patient 3: 6-9
      should_error = FALSE
    )
  )

  cat(
    "Running multi-patient comprehensive stress test with",
    length(test_cases),
    "test cases\n"
  )
  cat(
    "Test scenarios: 2 patients + single cutoff, 2 patients + multiple cutoffs, 3 patients + single cutoff\n"
  )

  # Process each test case individually
  for (case_idx in seq_along(test_cases)) {
    case <- test_cases[[case_idx]]

    cat("\n=== Testing Case", case_idx, ":", case$case_name, " ===\n")
    cat(
      "Patients:",
      case$n_patients,
      "Cutoffs:",
      case$n_cutoffs,
      "Total visits:",
      case$n_visits,
      "\n"
    )

    # Set up Stan data for this case
    stan_data <- list(
      N_cases = 1L,
      max_n_patients = case$n_patients,
      max_n_cutoffs = case$n_cutoffs,
      max_n_visits = case$n_visits,
      n_patients = array(case$n_patients),
      n_cutoffs = array(case$n_cutoffs),
      n_visits = array(case$n_visits),
      expect_error = array(0L), # All cases should succeed
      oos_patient_idx = array(case$oos_patient_idx, dim = c(1, case$n_cutoffs)),
      last_visit_calendar_day_sort_idx = array(
        case$last_visit_calendar_day_sort_idx,
        dim = c(1, case$n_patients)
      ),
      cutoff_calendar_day = array(
        case$cutoff_calendar_day,
        dim = c(1, case$n_cutoffs)
      ),
      patient_calendar_day = array(
        case$patient_calendar_day,
        dim = c(1, case$n_patients)
      ),
      t_patient_visits = array(
        case$t_patient_visits,
        dim = c(1, case$n_visits)
      ),
      t_patient_visits_day = array(
        case$t_patient_visits_day,
        dim = c(1, case$n_visits)
      ),
      patient_visit_pos = array(
        case$patient_visit_pos,
        dim = c(1, case$n_patients + 1)
      )
    )

    cat("Running case", case_idx, ":", case$case_name, "\n")

    # Run and validate
    fit <- test_stan_function(
      "tests/testthat/stan/test_get_testing_visit_week_bounds_all.stan",
      data = stan_data
    )

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
    expect_equal(case_status, 0, label = str_glue("Case {case_idx} should not error"))

    for (n in 1:case$n_cutoffs) {
      for (i in 1:case$n_patients) {
        expect_equal(
          get_val("first_testing_visit_week", 1L, n, i),
          expected$first_testing_visit_week[n, i],
          label = str_glue("Case {case_idx} {case$case_name} first_testing_visit_week n= {n} i= {i}")
        )
        expect_equal(
          get_val("testing_start_idx", 1L, n, i),
          expected$testing_start_idx[n, i],
          label = str_glue("Case {case_idx} {case$case_name} testing_start_idx n= {n} i= {i}")
        )
        for (m in 1:case$n_cutoffs) {
          expect_equal(
            get_val("last_testing_visit_week", 1L, n, m, i),
            expected$last_testing_visit_week[n, m, i],
            label = str_glue("Case {case_idx} {case$case_name} last_testing_visit_week n= {n} m= {m} i= {i}")
          )
          expect_equal(
            get_val("testing_end_idx", 1L, n, m, i),
            expected$testing_end_idx[n, m, i],
            label = str_glue("Case {case_idx} {case$case_name} testing_end_idx n= {n} m= {m} i= {i}")
          )
        }
      }
    }

    cat("Case", case_idx, "completed successfully!\n")
  }

  cat("\n*** MULTI-PATIENT COMPREHENSIVE STRESS TESTS PASSED! ***\n")
  cat("Successfully validated:\n")
  cat("- 2 patients with single cutoff ✓\n")
  cat("- 2 patients with multiple cutoffs ✓\n")
  cat("- 3 patients with single cutoff ✓\n")
  cat(
    "Function get_testing_visit_week_bounds correctly handles multiple patients and cutoffs!\n"
  )
})
