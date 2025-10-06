# tests/testthat/test-stan-estimate_kaplan_meier-array_int-array_int-int-int.R
# 
# IMPORTANT: Stan's estimate_kaplan_meier function uses a specific output convention:
# - The returned survival vector contains S(0), S(1), S(2), ..., S(max_t) 
# - S(0) accounts for immediate events: S(0) = (n_patients - immediate_events) / n_patients
# - If no immediate events (event_time=0), then S(0) = 1.0
# - This differs from R's survfit which typically reports survival only at event times
#
# PARAMETER CONVENTION:
# - event_time represents the TIME AT WHICH events occur (not survival duration)
# - Example: event_time=3 means event occurs at time 3, not after surviving 3 time units
#
library(here)
library(testthat)

# Helper function to run Stan test and extract results
run_km_test <- function(event_time, right_censored, max_t, description = "") {
  data <- list(
    n_patients = length(event_time),
    event_time = event_time,
    right_censored = right_censored,
    max_t = max_t,
    pfs_offset = 0
  )
  
  fit <- test_stan_function(
    here::here("tests", "testthat", "stan", "test_estimate_kaplan_meier-array_int-array_int-int-int.stan"),
    data
  )
  
  list(
    survival = as.numeric(fit$draws("km_survival", format = "matrix")[1, ]),
    at_risk = as.numeric(fit$draws("at_risk", format = "matrix")[1, ]),
    n_right_censored = as.numeric(fit$draws("n_right_censored", format = "matrix")[1, ]),
    n_exited = as.numeric(fit$draws("n_exited", format = "matrix")[1, ]),
    data = data
  )
}

# Helper function to validate basic KM properties
validate_km_properties <- function(result, description = "") {
  survival <- result$survival
  at_risk <- result$at_risk
  n_exited <- result$n_exited
  n_right_censored <- result$n_right_censored
  data <- result$data
  
  # Basic survival properties
  expect_length(survival, data$max_t + 1)
  expect_true(all(survival >= 0 & survival <= 1), 
              info = paste(description, "- Survival probabilities should be between 0 and 1"))
  expect_true(all(diff(survival) <= 0), 
              info = paste(description, "- Survival should be non-increasing"))
  
  # Array dimensions
  expect_length(at_risk, data$max_t + 1)
  expect_length(n_exited, data$max_t + 1)
  expect_length(n_right_censored, data$max_t + 1)
  
  # Non-negative counts
  expect_true(all(at_risk >= 0))
  expect_true(all(n_exited >= 0))
  expect_true(all(n_right_censored >= 0))
  
  # Initial at-risk should equal total patients
  expect_equal(at_risk[1], data$n_patients, 
               info = paste(description, "- All patients should be at risk initially"))
}

test_that("estimate_kaplan_meier basic functionality", {
  skip_if_not_installed("cmdstanr")
  
  # Test with random data
  data <- create_mock_km_data(n_patients = 20, max_t = 50, pfs_offset = 0)
  result <- run_km_test(data$event_time, data$right_censored, data$max_t, "Random data")
  validate_km_properties(result, "Random data")
  
  # Check immediate events handling
  immediate_events <- sum(data$event_time == 0 & data$right_censored == 0)
  if (immediate_events > 0) {
    expected_s0 <- (data$n_patients - immediate_events) / data$n_patients
    expect_equal(result$survival[1], expected_s0, tolerance = 1e-6, 
                 info = "S(0) should account for immediate events")
  } else {
    expect_equal(result$survival[1], 1.0, tolerance = 1e-6,
                 info = "S(0) should be 1.0 when no immediate events")
  }
})

test_that("estimate_kaplan_meier edge cases", {
  skip_if_not_installed("cmdstanr")
  
  # Test 1: Immediate events (event_time=0)
  result1 <- run_km_test(c(0, 0, 1, 2, 3), c(0, 0, 0, 1, 0), 5, "Immediate events")
  validate_km_properties(result1, "Immediate events")
  expect_true(result1$survival[1] < 1.0, info = "S(0) should be < 1 with immediate events")
  
  # Test 2: All events at same time
  result2 <- run_km_test(c(2, 2, 2, 2), c(0, 0, 0, 0), 5, "Simultaneous events")
  validate_km_properties(result2, "Simultaneous events")
  
  # Test 3: All censored
  result3 <- run_km_test(c(1, 3, 5), c(1, 1, 1), 6, "All censored")
  validate_km_properties(result3, "All censored")
  expect_true(all(result3$survival == 1), info = "Survival should stay at 1 with all censored")
  expect_true(all(result3$n_exited == 0), info = "No events with all censored")
  
  # Test 4: Single patient
  result4 <- run_km_test(c(3), c(0), 5, "Single patient")
  validate_km_properties(result4, "Single patient")
  expected_survival4 <- c(1, 1, 1, 0, 0, 0)  # Event at time 3
  expect_equal(result4$survival, expected_survival4, tolerance = 1e-6,
               info = "Single patient: survival=1 until event at time 3")
})

test_that("estimate_kaplan_meier algorithmic correctness", {
  skip_if_not_installed("cmdstanr")
  
  # Test with simple known case: 3 patients, events at times 1, 2, 3
  result <- run_km_test(c(1, 2, 3), c(0, 0, 0), 5, "Sequential events")
  
  # Manual calculation:
  # S(0) = 1.0 (no immediate events)
  # S(1) = 1.0 * (3-1)/3 = 2/3 = 0.667
  # S(2) = (2/3) * (2-1)/2 = 1/3 = 0.333  
  # S(3) = (1/3) * (1-1)/1 = 0.0
  expected_survival <- c(1.0, 2/3, 1/3, 0.0, 0.0, 0.0)
  expect_equal(result$survival, expected_survival, tolerance = 1e-6,
               info = "Sequential events should follow standard KM formula")
  
  # Test at-risk counts
  expected_at_risk <- c(3, 3, 2, 1, 0, 0)
  expect_equal(result$at_risk, expected_at_risk,
               info = "At-risk counts should decrease with events")
  
  # Test event counts
  expected_n_exited <- c(0, 1, 1, 1, 0, 0)
  expect_equal(result$n_exited, expected_n_exited,
               info = "Event counts should match input data")
})

test_that("estimate_kaplan_meier handles non-monotonic and complex data", {
  skip_if_not_installed("cmdstanr")
  
  # Non-monotonic event time values
  result1 <- run_km_test(c(5, 1, 8, 2, 6, 3), c(0, 0, 1, 0, 0, 0), 10, "Non-monotonic")
  validate_km_properties(result1, "Non-monotonic")
  
  total_events <- sum(1 - c(0, 0, 1, 0, 0, 0))
  expect_equal(sum(result1$n_exited), total_events,
               info = "Total events should match regardless of event time ordering")
  
  # Mixed events and censoring with immediate events
  result2 <- run_km_test(c(0, 0, 1, 2, 3, 5, 7, 10), 
                        c(0, 0, 0, 1, 0, 0, 1, 0), 12, "Complex mixed")
  validate_km_properties(result2, "Complex mixed")
  expect_true(result2$survival[1] < 1.0, info = "S(0) should account for immediate events")
})

test_that("estimate_kaplan_meier comparison with R survfit", {
  skip_if_not_installed("cmdstanr")
  skip_if_not_installed("ggsurvfit") 
  skip_if_not_installed("survival")
  
  # Simple test case for comparison
  event_time <- c(1, 2, 3, 4)
  right_censored <- c(0, 0, 0, 0)
  
  # Stan result
  stan_result <- run_km_test(event_time, right_censored, 6, "R comparison")
  
  # R result
  r_data <- data.frame(time = event_time, event = 1 - right_censored)
  r_survfit <- ggsurvfit::survfit2(survival::Surv(time, event) ~ 1, data = r_data)
  r_summary <- summary(r_survfit)
  
  # Compare at event times (account for indexing difference)
  # R gives survival AT event times, Stan gives survival AFTER processing each time
  r_event_times <- r_summary$time
  r_survival_values <- r_summary$surv
  
  for (i in seq_along(r_event_times)) {
    event_time <- r_event_times[i]
    # Stan survival after processing events at time t is in position t+1
    stan_survival_at_time <- stan_result$survival[event_time + 1]
    expect_equal(stan_survival_at_time, r_survival_values[i], tolerance = 1e-6,
                 info = sprintf("Stan and R should match at event time %d", event_time))
  }
})

test_that("estimate_kaplan_meier mathematical consistency", {
  skip_if_not_installed("cmdstanr")
  
  # Test mathematical relationships between tuple components
  test_scenarios <- list(
    list(event_time = c(0, 1, 3, 8, 12, 15), right_censored = c(0, 0, 1, 0, 1, 0), max_t = 20),
    list(event_time = c(2, 4, 6, 8, 10), right_censored = c(1, 1, 1, 1, 0), max_t = 15),
    list(event_time = c(3, 3, 3, 7, 7, 10, 12), right_censored = c(0, 0, 1, 0, 0, 1, 0), max_t = 15)
  )
  
  for (i in seq_along(test_scenarios)) {
    scenario <- test_scenarios[[i]]
    result <- run_km_test(scenario$event_time, scenario$right_censored, scenario$max_t, 
                         paste("Consistency scenario", i))
    
    # Validate relationships
    n <- length(scenario$event_time)
    
    # At start, all patients at risk
    expect_equal(result$at_risk[1], n)
    
    # Total events should match input
    total_events <- sum(1 - scenario$right_censored)
    expect_equal(sum(result$n_exited), total_events)
    
    # Total censored should match input  
    total_censored <- sum(scenario$right_censored)
    expect_equal(sum(result$n_right_censored), total_censored)
    
    # At-risk should decrease monotonically (when events/censoring occur)
    risk_changes <- which(result$n_exited > 0 | result$n_right_censored > 0)
    if (length(risk_changes) > 1) {
      for (j in 2:length(risk_changes)) {
        expect_lte(result$at_risk[risk_changes[j]], result$at_risk[risk_changes[j-1]])
      }
    }
  }
})
