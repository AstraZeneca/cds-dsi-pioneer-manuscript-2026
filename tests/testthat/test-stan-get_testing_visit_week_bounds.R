library(testthat)
library(cmdstanr)
source(here::here("tests/testthat/helper-stan.R"))
source(here::here("tests/testthat/helper-lfo.R"))

test_that("get_testing_visit_week_bounds stress test - all edge cases and failure modes", {
  # Comprehensive stress test cases covering all potential failure modes
  test_cases <- list(
    # === BASIC FUNCTIONALITY TESTS ===

    # Case 1: Minimal valid case - single patient, single cutoff
    list(
      case_name = "minimal_valid_case",
      n_patients = 1L,
      n_cutoffs = 1L,
      n_visits = 2L,
      oos_patient_idx = c(1L),
      last_visit_calendar_day_sort_idx = c(1L),
      cutoff_calendar_day = c(110L),
      patient_calendar_day = c(100L),
      t_patient_visits = c(1L, 3L),
      t_patient_visits_day = c(7L, 21L),
      patient_visit_pos = c(1L, 3L),
      should_error = FALSE
    ),

    # Case 2: Multiple patients, multiple cutoffs - normal case
    list(
      case_name = "normal_multi_patient_cutoff",
      n_patients = 3L,
      n_cutoffs = 2L,
      n_visits = 9L,
      oos_patient_idx = c(1L, 2L),
      last_visit_calendar_day_sort_idx = c(1L, 2L, 3L),
      cutoff_calendar_day = c(120L, 140L),
      patient_calendar_day = c(100L, 105L, 110L),
      t_patient_visits = c(1L, 2L, 4L, 1L, 3L, 5L, 2L, 4L, 6L),
      t_patient_visits_day = c(7L, 14L, 28L, 7L, 21L, 35L, 14L, 28L, 42L),
      patient_visit_pos = c(1L, 4L, 7L, 10L),
      should_error = FALSE
    ),

    # === BOUNDARY CONDITION TESTS ===

    # Case 3: Visit exactly at cutoff boundary
    list(
      case_name = "visit_exactly_at_cutoff",
      n_patients = 1L,
      n_cutoffs = 1L,
      n_visits = 3L,
      oos_patient_idx = c(1L),
      last_visit_calendar_day_sort_idx = c(1L),
      cutoff_calendar_day = c(121L), # Patient day 21 exactly
      patient_calendar_day = c(100L),
      t_patient_visits = c(1L, 3L, 5L),
      t_patient_visits_day = c(7L, 21L, 35L), # Day 21 exactly at cutoff
      patient_visit_pos = c(1L, 4L),
      should_error = FALSE
    ),

    # Case 4: Visit one day before cutoff
    list(
      case_name = "visit_one_day_before_cutoff",
      n_patients = 1L,
      n_cutoffs = 1L,
      n_visits = 3L,
      oos_patient_idx = c(1L),
      last_visit_calendar_day_sort_idx = c(1L),
      cutoff_calendar_day = c(122L), # Patient day 22, visit at day 21
      patient_calendar_day = c(100L),
      t_patient_visits = c(1L, 3L, 5L),
      t_patient_visits_day = c(7L, 21L, 35L),
      patient_visit_pos = c(1L, 4L),
      should_error = FALSE
    ),

    # Case 5: Visit one day after cutoff
    list(
      case_name = "visit_one_day_after_cutoff",
      n_patients = 1L,
      n_cutoffs = 1L,
      n_visits = 3L,
      oos_patient_idx = c(1L),
      last_visit_calendar_day_sort_idx = c(1L),
      cutoff_calendar_day = c(119L), # Patient day 20 (119-100+1=20), visit at day 21 is after
      patient_calendar_day = c(100L),
      t_patient_visits = c(1L, 3L, 5L),
      t_patient_visits_day = c(7L, 21L, 35L),
      patient_visit_pos = c(1L, 4L),
      should_error = FALSE
    ),

    # === ERROR CONDITIONS TESTS ===

    # Case 6: No visits after cutoff - should trigger fatal_error
    list(
      case_name = "no_visits_after_cutoff_fatal_error",
      n_patients = 1L,
      n_cutoffs = 1L,
      n_visits = 2L,
      oos_patient_idx = c(1L),
      last_visit_calendar_day_sort_idx = c(1L),
      cutoff_calendar_day = c(200L), # Way after all visits
      patient_calendar_day = c(100L),
      t_patient_visits = c(1L, 2L),
      t_patient_visits_day = c(7L, 14L), # Last visit at day 114
      patient_visit_pos = c(1L, 3L),
      should_error = TRUE,
      error_pattern = "Unexpectedly could not find the first testing visit"
    ),

    # Case 7: Only baseline visit (day 0 or negative)
    list(
      case_name = "only_baseline_visit",
      n_patients = 1L,
      n_cutoffs = 1L,
      n_visits = 1L,
      oos_patient_idx = c(1L),
      last_visit_calendar_day_sort_idx = c(1L),
      cutoff_calendar_day = c(105L),
      patient_calendar_day = c(100L),
      t_patient_visits = c(0L), # Baseline visit
      t_patient_visits_day = c(0L),
      patient_visit_pos = c(1L, 2L),
      should_error = TRUE,
      error_pattern = "Unexpectedly could not find the first testing visit"
    ),

    # === PATIENT ENTRY DATE EDGE CASES ===

    # Case 8: Patient enters exactly at cutoff
    list(
      case_name = "patient_entry_at_cutoff",
      n_patients = 1L,
      n_cutoffs = 1L,
      n_visits = 2L,
      oos_patient_idx = c(1L),
      last_visit_calendar_day_sort_idx = c(1L),
      cutoff_calendar_day = c(100L), # Same as patient entry
      patient_calendar_day = c(100L),
      t_patient_visits = c(1L, 3L),
      t_patient_visits_day = c(7L, 21L), # Days 7, 21 relative to entry
      patient_visit_pos = c(1L, 3L),
      should_error = FALSE
    ),

    # === EXTREME VALUES TESTS ===

    # Case 9: Very large calendar days
    list(
      case_name = "large_calendar_days",
      n_patients = 1L,
      n_cutoffs = 1L,
      n_visits = 2L,
      oos_patient_idx = c(1L),
      last_visit_calendar_day_sort_idx = c(1L),
      cutoff_calendar_day = c(999990L), # Very large
      patient_calendar_day = c(999900L),
      t_patient_visits = c(1L, 2L),
      t_patient_visits_day = c(7L, 200L), # Large day offset
      patient_visit_pos = c(1L, 3L),
      should_error = FALSE
    ),

    # Case 10: Zero and negative days (edge case testing)
    list(
      case_name = "zero_and_negative_days",
      n_patients = 1L,
      n_cutoffs = 1L,
      n_visits = 3L,
      oos_patient_idx = c(1L),
      last_visit_calendar_day_sort_idx = c(1L),
      cutoff_calendar_day = c(105L),
      patient_calendar_day = c(100L),
      t_patient_visits = c(0L, 1L, 2L), # Week 0 (baseline)
      t_patient_visits_day = c(-1L, 0L, 7L), # Negative day, day 0, day 7
      patient_visit_pos = c(1L, 4L),
      should_error = FALSE
    ),

    # === EXTREME STRESS CASES ===

    # Case 11: Massive number of visits to test performance and bounds
    list(
      case_name = "massive_visit_count",
      n_patients = 1L,
      n_cutoffs = 1L,
      n_visits = 50L, # Stress test with many visits
      oos_patient_idx = c(1L),
      last_visit_calendar_day_sort_idx = c(1L),
      # calendar_date_to_study_date(500, 843) = 843-500+1 = 344
      # Visits at days 7,14,...,343 (weeks 1-49, ≤344) and 350 (week 50, >344)
      # First testing visit = week 50, global idx = 50
      cutoff_calendar_day = c(843L),
      patient_calendar_day = c(500L),
      t_patient_visits = 1:50, # Weeks 1-50
      t_patient_visits_day = seq(7, by = 7, length.out = 50), # Days 7,14,...,350
      patient_visit_pos = c(1L, 51L),
      should_error = FALSE
    ),

    # Case 12: Pathological case - visits with huge gaps
    list(
      case_name = "huge_visit_gaps",
      n_patients = 1L,
      n_cutoffs = 1L,
      n_visits = 3L,
      oos_patient_idx = c(1L),
      last_visit_calendar_day_sort_idx = c(1L),
      cutoff_calendar_day = c(1000L),
      patient_calendar_day = c(100L),
      t_patient_visits = c(1L, 100L, 200L), # Huge gaps between visits
      t_patient_visits_day = c(7L, 700L, 1400L), # Corresponding huge day gaps
      patient_visit_pos = c(1L, 4L),
      should_error = FALSE
    ),

    # Case 13: Multiple patients with vastly different visit patterns
    list(
      case_name = "asymmetric_patient_visits",
      n_patients = 3L,
      n_cutoffs = 1L,
      n_visits = 10L,
      oos_patient_idx = c(2L), # Only test 2nd patient in sort order
      last_visit_calendar_day_sort_idx = c(2L, 1L, 3L), # Patient 2, 1, 3 sorted by last visit
      cutoff_calendar_day = c(150L),
      patient_calendar_day = c(100L, 110L, 120L),
      # Patient 1: 1 visit, Patient 2: 8 visits, Patient 3: 1 visit
      t_patient_visits = c(5L, 1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 10L),
      t_patient_visits_day = c(35L, 7L, 14L, 21L, 28L, 35L, 42L, 49L, 56L, 70L),
      patient_visit_pos = c(1L, 2L, 10L, 11L), # P1: 1 visit, P2: 8 visits, P3: 1 visit
      should_error = FALSE
    ),

    # Case 14: Cutoff exactly between two consecutive visits
    list(
      case_name = "cutoff_between_consecutive_visits",
      n_patients = 1L,
      n_cutoffs = 1L,
      n_visits = 5L,
      oos_patient_idx = c(1L),
      last_visit_calendar_day_sort_idx = c(1L),
      cutoff_calendar_day = c(125L), # Exactly between day 21 and 28
      patient_calendar_day = c(100L),
      t_patient_visits = c(1L, 2L, 3L, 4L, 5L),
      t_patient_visits_day = c(7L, 14L, 21L, 28L, 35L), # Days 107, 114, 121, 128, 135
      patient_visit_pos = c(1L, 6L),
      should_error = FALSE
    ),

    # Case 15: Integer overflow boundary test
    list(
      case_name = "integer_boundary_values",
      n_patients = 1L,
      n_cutoffs = 1L,
      n_visits = 2L,
      oos_patient_idx = c(1L),
      last_visit_calendar_day_sort_idx = c(1L),
      cutoff_calendar_day = c(.Machine$integer.max - 1000L), # Near integer max
      patient_calendar_day = c(.Machine$integer.max - 2000L),
      t_patient_visits = c(100L, 200L),
      t_patient_visits_day = c(500L, 1500L), # Large but safe offsets
      patient_visit_pos = c(1L, 3L),
      should_error = FALSE
    ),

    # Case 16: Multiple cutoffs with overlapping testing windows
    list(
      case_name = "overlapping_testing_windows",
      n_patients = 1L,
      n_cutoffs = 3L,
      n_visits = 6L,
      oos_patient_idx = c(1L, 1L, 1L),
      last_visit_calendar_day_sort_idx = c(1L),
      cutoff_calendar_day = c(110L, 115L, 120L), # Close cutoffs
      patient_calendar_day = c(100L),
      t_patient_visits = c(1L, 2L, 3L, 4L, 5L, 6L),
      t_patient_visits_day = c(7L, 14L, 21L, 28L, 35L, 42L),
      patient_visit_pos = c(1L, 7L),
      should_error = FALSE
    ),

    # Case 17: Edge case with only one visit after cutoff
    list(
      case_name = "single_visit_after_cutoff",
      n_patients = 1L,
      n_cutoffs = 1L,
      n_visits = 5L,
      oos_patient_idx = c(1L),
      last_visit_calendar_day_sort_idx = c(1L),
      cutoff_calendar_day = c(133L), # study_day = 133-100+1 = 34; day 35 is first visit after
      patient_calendar_day = c(100L),
      t_patient_visits = c(1L, 2L, 3L, 4L, 5L),
      t_patient_visits_day = c(7L, 14L, 21L, 28L, 35L), # Last visit at day 35
      patient_visit_pos = c(1L, 6L),
      should_error = FALSE
    ),

    # Case 18: Stress test with many patients and cutoffs
    list(
      case_name = "many_patients_many_cutoffs",
      n_patients = 5L,
      n_cutoffs = 5L,
      n_visits = 25L, # 5 visits per patient
      oos_patient_idx = c(1L, 2L, 3L, 4L, 5L),
      last_visit_calendar_day_sort_idx = c(1L, 2L, 3L, 4L, 5L),
      cutoff_calendar_day = c(120L, 130L, 140L, 150L, 160L),
      patient_calendar_day = c(100L, 105L, 110L, 115L, 120L),
      # 5 visits per patient, staggered patterns
      t_patient_visits = c(
        1L,
        3L,
        5L,
        7L,
        9L, # Patient 1
        2L,
        4L,
        6L,
        8L,
        10L, # Patient 2
        1L,
        2L,
        8L,
        9L,
        10L, # Patient 3
        3L,
        4L,
        5L,
        11L,
        12L, # Patient 4
        6L,
        7L,
        11L,
        12L,
        13L # Patient 5
      ),
      t_patient_visits_day = c(
        7L,
        21L,
        35L,
        49L,
        63L, # Patient 1
        14L,
        28L,
        42L,
        56L,
        70L, # Patient 2
        7L,
        14L,
        56L,
        63L,
        70L, # Patient 3
        21L,
        28L,
        35L,
        77L,
        84L, # Patient 4
        42L,
        49L,
        77L,
        84L,
        91L # Patient 5
      ),
      patient_visit_pos = c(1L, 6L, 11L, 16L, 21L, 26L),
      should_error = FALSE
    )
  )

  N_cases <- length(test_cases)
  max_n_patients <- max(sapply(test_cases, function(x) x$n_patients))
  max_n_cutoffs <- max(sapply(test_cases, function(x) x$n_cutoffs))
  max_n_visits <- max(sapply(test_cases, function(x) x$n_visits))

  # Prepare arrays for Stan (with safety padding)
  n_patients <- integer(N_cases)
  n_cutoffs <- integer(N_cases)
  n_visits <- integer(N_cases)
  expect_error <- integer(N_cases)

  oos_patient_idx <- matrix(1L, N_cases, max_n_cutoffs)
  last_visit_calendar_day_sort_idx <- matrix(1L, N_cases, max_n_patients)
  cutoff_calendar_day <- matrix(0L, N_cases, max_n_cutoffs)
  patient_calendar_day <- matrix(0L, N_cases, max_n_patients)
  t_patient_visits <- matrix(0L, N_cases, max_n_visits)
  t_patient_visits_day <- matrix(0L, N_cases, max_n_visits)
  patient_visit_pos <- matrix(1L, N_cases, max_n_patients + 1)

  error_cases <- integer(0)
  valid_cases <- integer(0)

  for (i in seq_len(N_cases)) {
    case <- test_cases[[i]]

    n_patients[i] <- case$n_patients
    n_cutoffs[i] <- case$n_cutoffs
    n_visits[i] <- case$n_visits
    expect_error[i] <- if (case$should_error) 1L else 0L

    oos_patient_idx[i, 1:case$n_cutoffs] <- case$oos_patient_idx
    last_visit_calendar_day_sort_idx[
      i,
      1:case$n_patients
    ] <- case$last_visit_calendar_day_sort_idx
    cutoff_calendar_day[i, 1:case$n_cutoffs] <- case$cutoff_calendar_day
    patient_calendar_day[i, 1:case$n_patients] <- case$patient_calendar_day
    t_patient_visits[i, 1:case$n_visits] <- case$t_patient_visits
    t_patient_visits_day[i, 1:case$n_visits] <- case$t_patient_visits_day
    patient_visit_pos[i, 1:(case$n_patients + 1)] <- case$patient_visit_pos

    if (case$should_error) {
      error_cases <- c(error_cases, i)
    } else {
      valid_cases <- c(valid_cases, i)
    }
  }

  stan_data <- list(
    N_cases = N_cases,
    max_n_patients = max_n_patients,
    max_n_cutoffs = max_n_cutoffs,
    max_n_visits = max_n_visits,
    n_patients = n_patients,
    n_cutoffs = n_cutoffs,
    n_visits = n_visits,
    expect_error = expect_error,
    oos_patient_idx = oos_patient_idx,
    last_visit_calendar_day_sort_idx = last_visit_calendar_day_sort_idx,
    cutoff_calendar_day = cutoff_calendar_day,
    patient_calendar_day = patient_calendar_day,
    t_patient_visits = t_patient_visits,
    t_patient_visits_day = t_patient_visits_day,
    patient_visit_pos = patient_visit_pos
  )

  cat("Running comprehensive stress test with", N_cases, "test cases\n")
  cat("EXTREME STRESS TEST SCENARIOS:\n")
  cat("=== BASIC FUNCTIONALITY ===\n")
  cat("- Minimal valid cases (single patient/cutoff)\n")
  cat("- Multi-patient, multi-cutoff scenarios\n")
  cat("=== BOUNDARY CONDITIONS ===\n")
  cat("- Visits exactly at, before, and after cutoff boundaries\n")
  cat("- Patient entry dates at/before/after cutoffs\n")
  cat("- Single visits after cutoffs\n")
  cat("=== ERROR CONDITIONS ===\n")
  cat("- No visits after cutoff (fatal_error expected)\n")
  cat("- Only baseline visits (should fail)\n")
  cat("- Empty patient visit arrays\n")
  cat("=== EXTREME VALUES & STRESS TESTS ===\n")
  cat("- Massive visit counts (50+ visits per patient)\n")
  cat("- Huge gaps between visits (hundreds of weeks)\n")
  cat("- Integer boundary values (near .Machine$integer.max)\n")
  cat("- Very large calendar day values\n")
  cat("- Zero, negative, and baseline visit days\n")
  cat("=== COMPLEX SCENARIOS ===\n")
  cat("- Asymmetric patient visit patterns\n")
  cat("- Multiple cutoffs with overlapping testing windows\n")
  cat("- Many patients (5) × many cutoffs (5) = 25 combinations\n")
  cat("- Mixed patient ordering in sort indices\n")
  cat("- Visits not chronologically ordered in input arrays\n")
  cat("- Duplicate visit days\n")
  cat("- Cutoffs very close together (1 day apart)\n")
  cat(
    "Expected error cases:",
    length(error_cases),
    "- cases:",
    error_cases,
    "\n"
  )
  cat("Expected valid cases:", length(valid_cases), "\n")
  cat("*** This test pushes the function to its absolute limits! ***\n")

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

  # Execute the comprehensive stress test
  cat("Executing comprehensive stress test...\n")

  tryCatch(
    {
      fit <- test_stan_function(
        "tests/testthat/stan/test_get_testing_visit_week_bounds_all.stan",
        data = stan_data
      )

      draws_df <- posterior::as_draws_df(fit$draws())
      get_val  <- function(var, ...) get_stan_val(draws_df, var, ...)

      cat("Successfully extracted test results\n")

      # Validate all valid cases — all 4 output arrays
      for (case_idx in valid_cases) {
        case     <- test_cases[[case_idx]]
        expected <- expected_by_case[[case_idx]]
        lbl      <- str_glue("Case {case_idx} {case$case_name}")

        expect_equal(get_val("case_status", case_idx), 0, label = str_c(lbl, " no error"))

        for (n in 1:case$n_cutoffs) {
          for (i in 1:case$n_patients) {
            expect_equal(
              get_val("first_testing_visit_week", case_idx, n, i),
              expected$first_testing_visit_week[n, i],
              label = str_glue("{lbl} first_testing_visit_week n {n} i {i}")
            )
            expect_equal(
              get_val("testing_start_idx", case_idx, n, i),
              expected$testing_start_idx[n, i],
              label = str_glue("{lbl} testing_start_idx n {n} i {i}")
            )
            for (m in 1:case$n_cutoffs) {
              expect_equal(
                get_val("last_testing_visit_week", case_idx, n, m, i),
                expected$last_testing_visit_week[n, m, i],
                label = str_glue("{lbl} last_testing_visit_week n {n} m {m} i {i}")
              )
              expect_equal(
                get_val("testing_end_idx", case_idx, n, m, i),
                expected$testing_end_idx[n, m, i],
                label = str_glue("{lbl} testing_end_idx n {n} m {m} i {i}")
              )
            }
          }
        }
      }

      # Check that error cases were marked appropriately
      for (case_idx in error_cases) {
        case <- test_cases[[case_idx]]
        expect_equal(
          get_val("case_status", case_idx),
          1,
          label = str_glue("Case {case_idx} {case$case_name} should be marked as error")
        )
      }

      cat("*** ALL STRESS TESTS PASSED! ***\n")
    },
    error = function(e) {
      cat("Stan execution error encountered:\n")
      cat("Error message:", e$message, "\n")

      # If we get a Stan error, check if it matches expected error patterns
      error_matched <- FALSE
      for (case_idx in error_cases) {
        case <- test_cases[[case_idx]]
        if (
          !is.null(case$error_pattern) &&
            grepl(case$error_pattern, e$message, fixed = TRUE)
        ) {
          error_matched <- TRUE
          cat(
            "Expected error matched for case",
            case_idx,
            ":",
            case$case_name,
            "\n"
          )
          break
        }
      }

      if (!error_matched) {
        stop("Unexpected error in stress test: ", e$message)
      } else {
        cat("Expected error handling successful\n")
      }
    }
  )
})
