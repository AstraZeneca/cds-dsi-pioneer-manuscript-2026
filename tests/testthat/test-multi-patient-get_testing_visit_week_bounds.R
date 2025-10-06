library(testthat)
library(cmdstanr)
source(here::here("tests/testthat/helper-stan.R"))

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
      should_error = FALSE,
      # Cutoff at day 115: Patient 1 cutoff at study day 16, Patient 2 cutoff at study day 11
      # Patient 1: first visit after day 16 should be day 21 (week 3) at global index 3
      # Patient 2: first visit after day 11 should be day 14 (week 2) at global index 5
      expected_first_testing_visit_week = matrix(c(3L, 2L), nrow=1, ncol=2),
      expected_testing_start_idx = matrix(c(3L, 5L), nrow=1, ncol=2)
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
      should_error = FALSE,
      # Cutoff 1 (day 115): P1 study day 16 -> week 3, P2 study day 11 -> week 2
      # Cutoff 2 (day 125): P1 study day 26 -> week 5, P2 study day 21 -> week 4
      expected_first_testing_visit_week = matrix(c(3L, 2L, 5L, 4L), nrow=2, ncol=2, byrow=TRUE),
      expected_testing_start_idx = matrix(c(3L, 6L, 4L, 8L), nrow=2, ncol=2, byrow=TRUE)
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
      should_error = FALSE,
      # Cutoff at day 120: P1 study day 21 -> week 4, P2 study day 11 -> week 3, P3 study day 16 -> week 3
      expected_first_testing_visit_week = matrix(c(4L, 3L, 3L), nrow=1, ncol=3),
      expected_testing_start_idx = matrix(c(2L, 4L, 8L), nrow=1, ncol=3)
    )
  )
  
  cat("Running multi-patient comprehensive stress test with", length(test_cases), "test cases\n")
  cat("Test scenarios: 2 patients + single cutoff, 2 patients + multiple cutoffs, 3 patients + single cutoff\n")
  
  # Process each test case individually
  for (case_idx in seq_along(test_cases)) {
    case <- test_cases[[case_idx]]
    
    cat("\n=== Testing Case", case_idx, ":", case$case_name, " ===\n")
    cat("Patients:", case$n_patients, "Cutoffs:", case$n_cutoffs, "Total visits:", case$n_visits, "\n")
    
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
      last_visit_calendar_day_sort_idx = array(case$last_visit_calendar_day_sort_idx, dim = c(1, case$n_patients)),
      cutoff_calendar_day = array(case$cutoff_calendar_day, dim = c(1, case$n_cutoffs)),
      patient_calendar_day = array(case$patient_calendar_day, dim = c(1, case$n_patients)),
      t_patient_visits = array(case$t_patient_visits, dim = c(1, case$n_visits)),
      t_patient_visits_day = array(case$t_patient_visits_day, dim = c(1, case$n_visits)),
      patient_visit_pos = array(case$patient_visit_pos, dim = c(1, case$n_patients + 1))
    )
    
    cat("Running case", case_idx, ":", case$case_name, "\n")
    
    # Run and validate
    fit <- test_stan_function("tests/testthat/stan/test_get_testing_visit_week_bounds_all.stan", data = stan_data)
    
    # Extract results
    first_testing_visit_week_draws <- fit$draws("first_testing_visit_week", format = "draws_array")
    testing_start_idx_draws <- fit$draws("testing_start_idx", format = "draws_array")
    case_status_draws <- fit$draws("case_status", format = "draws_array")
    
    cat("Stan execution successful!\n")
    cat("Result dimensions - first_testing_visit_week:", dim(first_testing_visit_week_draws), "\n")
    
    # Extract the actual results based on actual dimensions
    case_status <- as.numeric(case_status_draws[1,1,1])
    expect_equal(case_status, 0, label = paste("Case", case_idx, "should not error"))
    
    # Validate specific results for each cutoff and patient
    for (n in 1:case$n_cutoffs) {
      for (i in 1:case$n_patients) {
        # For multi-patient and/or multi-cutoff, extract based on dimensions
        if (case$n_cutoffs > 1 && case$n_patients > 1) {
          # Both dimensions present: [draws, chains, cutoffs, patients]
          first_week_actual <- as.numeric(first_testing_visit_week_draws[1,1,n,i])
          start_idx_actual <- as.numeric(testing_start_idx_draws[1,1,n,i])
        } else if (case$n_cutoffs > 1) {
          # Multiple cutoffs, single patient: [draws, chains, cutoffs]
          first_week_actual <- as.numeric(first_testing_visit_week_draws[1,1,n])
          start_idx_actual <- as.numeric(testing_start_idx_draws[1,1,n])
        } else {
          # Single cutoff, multiple patients: [draws, chains, patients]
          first_week_actual <- as.numeric(first_testing_visit_week_draws[1,1,i])
          start_idx_actual <- as.numeric(testing_start_idx_draws[1,1,i])
        }
        
        expect_equal(
          first_week_actual, 
          case$expected_first_testing_visit_week[n, i],
          label = paste("Case", case_idx, case$case_name, "first_testing_visit_week n=", n, "i=", i)
        )
        
        expect_equal(
          start_idx_actual, 
          case$expected_testing_start_idx[n, i],
          label = paste("Case", case_idx, case$case_name, "testing_start_idx n=", n, "i=", i)
        )
        
        cat("✓ Cutoff", n, "Patient", i, ": first_week =", first_week_actual, ", start_idx =", start_idx_actual, "\n")
      }
    }
    
    cat("Case", case_idx, "completed successfully!\n")
  }
  
  cat("\n*** MULTI-PATIENT COMPREHENSIVE STRESS TESTS PASSED! ***\n")
  cat("Successfully validated:\n")
  cat("- 2 patients with single cutoff ✓\n")
  cat("- 2 patients with multiple cutoffs ✓\n") 
  cat("- 3 patients with single cutoff ✓\n")
  cat("Function get_testing_visit_week_bounds correctly handles multiple patients and cutoffs!\n")
})
