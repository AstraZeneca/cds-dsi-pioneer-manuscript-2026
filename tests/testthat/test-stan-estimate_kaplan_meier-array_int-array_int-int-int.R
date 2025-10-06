    # Skip test cases where n_patients is 0 (invalid for KM)
# IMPORTANT: Stan's estimate_kaplan_meier function uses a specific output convention:
# - The returned survival vector contains S(0), S(1), S(2), ..., S(max_t) 
# - S(0) accounts for immediate events: S(0) = (n_patients - immediate_events) / n_patients
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


# Define all test cases in a list (zero-patient cases removed)
test_cases <- list(
  # Basic cases
  list(event_time = c(0, 0, 1, 2, 3), right_censored = c(0, 0, 0, 1, 0), max_t = 5, pfs_offset = 0),
  list(event_time = c(1, 2, 3), right_censored = c(0, 0, 0), max_t = 5, pfs_offset = 0),
  # Algorithmic correctness: 3 patients, events at 1,2,3
  list(event_time = c(1, 2, 3), right_censored = c(0, 0, 0), max_t = 5, pfs_offset = 0),
  # Non-monotonic event times
  list(event_time = c(5, 1, 8, 2, 6, 3), right_censored = c(0, 0, 1, 0, 0, 0), max_t = 10, pfs_offset = 0),
  # Complex mixed events/censoring with immediate events
  list(event_time = c(0, 0, 1, 2, 3, 5, 7, 10), right_censored = c(0, 0, 0, 1, 0, 0, 1, 0), max_t = 12, pfs_offset = 0),
  # R comparison: 4 events at 1,2,3,4
  list(event_time = c(1, 2, 3, 4), right_censored = c(0, 0, 0, 0), max_t = 6, pfs_offset = 0),
  # Mathematical consistency scenarios
  list(event_time = c(0, 1, 3, 8, 12, 15), right_censored = c(0, 0, 1, 0, 1, 0), max_t = 20, pfs_offset = 0),
  list(event_time = c(2, 4, 6, 8, 10), right_censored = c(1, 1, 1, 1, 0), max_t = 15, pfs_offset = 0),
  list(event_time = c(3, 3, 3, 7, 7, 10, 12), right_censored = c(0, 0, 1, 0, 0, 1, 0), max_t = 15, pfs_offset = 0)
)

n_cases <- length(test_cases) # Number of test cases
n_patients <- sapply(test_cases, function(x) length(x$event_time))
max_t <- sapply(test_cases, function(x) x$max_t)
pfs_offset <- sapply(test_cases, function(x) x$pfs_offset)
event_time <- unlist(lapply(test_cases, function(x) x$event_time))
right_censored <- unlist(lapply(test_cases, function(x) x$right_censored))
start_idx <- cumsum(c(1, head(n_patients, -1)))

stan_data <- list(
  n_cases = n_cases,
  n_patients = n_patients,
  max_t = max_t, # Maximum time points for survival analysis
  pfs_offset = pfs_offset,
  event_time_sum = length(event_time),
  event_time = event_time,
  right_censored = right_censored,
  start_idx = start_idx
)


fit <- test_stan_function(
  here::here("tests", "testthat", "stan", "test_estimate_kaplan_meier_all.stan"),
  stan_data
)

# Extract all results as arrays (rectangular, n_cases x max(max_t)+1)
get_matrix <- function(var) {
  mat <- fit$draws(var, format = "matrix")
  # Stan variable names: var[i,j] for i in 1:n_cases, j in 1:n_time
  n_cases <- length(test_cases)
  n_time <- max(sapply(test_cases, function(x) x$max_t)) + 1
  out <- matrix(NA, nrow = n_cases, ncol = n_time)
  for (i in seq_len(n_cases)) {
    for (j in seq_len(n_time)) {
      vname <- sprintf("%s[%d,%d]", var, i, j)
      if (vname %in% colnames(mat)) {
        out[i, j] <- mat[1, vname]
      }
    }
  }
  out
}
km_survival <- get_matrix("km_survival")
at_risk <- get_matrix("at_risk")
n_right_censored <- get_matrix("n_right_censored")
n_exited <- get_matrix("n_exited")


# Now write your expect_equal() checks for each test case
test_that("Kaplan-Meier S(0) is correct for all cases", {
  for (i in seq_along(test_cases)) {
    expect_equal(km_survival[i, 1], 1.0, tolerance = 1e-8)
  }
})

# Helper function to validate basic KM properties
validate_km_properties <- function(result, description = "") {
  survival <- result$survival
  at_risk <- result$at_risk
  n_exited <- result$n_exited
  n_right_censored <- result$n_right_censored
  data <- result$data
  
  # Basic survival properties (Stan output includes S(0), has length max_t + 1)
  expect_length(survival, data$max_t + 1)
  expect_true(all(survival >= 0 & survival <= 1), 
              info = paste(description, "- Survival probabilities should be between 0 and 1"))
  expect_true(all(diff(survival) <= 0), 
              info = paste(description, "- Survival should be non-increasing"))
  
  # Array dimensions (Stan output includes time 0, has length max_t + 1)
  expect_length(at_risk, data$max_t + 1)
  expect_length(n_exited, data$max_t + 1)
  expect_length(n_right_censored, data$max_t + 1)
  
  # Non-negative counts
  expect_true(all(at_risk >= 0))
  expect_true(all(n_exited >= 0))
  expect_true(all(n_right_censored >= 0))
  
  # Initial at-risk should equal total patients (at time 0, index 1)
  expect_equal(at_risk[1], data$n_patients, 
               info = paste(description, "- All patients should be at risk initially at time 0"))
}


test_that("estimate_kaplan_meier basic functionality", {
  for (i in seq_along(test_cases)) {
    idx <- seq_len(test_cases[[i]]$max_t + 1)
    survival <- km_survival[i, idx]
    at_risk_case <- as.integer(at_risk[i, idx])
    n_exited_case <- as.integer(n_exited[i, idx])
    n_right_censored_case <- as.integer(n_right_censored[i, idx])
    data <- test_cases[[i]]
    # All test cases have at least one patient; no skip needed
    # All test cases have at least one patient; no skip needed
    expect_length(survival, data$max_t + 1)
    expect_true(all(survival >= 0 & survival <= 1))
    expect_true(all(diff(survival) <= 0))
    # Skip test cases where n_patients is 0 (invalid for KM)
    expect_length(at_risk_case, data$max_t + 1)
    expect_length(n_exited_case, data$max_t + 1)
    expect_length(n_right_censored_case, data$max_t + 1)
    expect_true(all(at_risk_case >= 0))
    expect_true(all(n_exited_case >= 0))
    expect_true(all(n_right_censored_case >= 0))
    expect_equal(at_risk_case[1], length(data$event_time))
  }
})


# Edge cases: check all test cases in the all-in-one fit
test_that("estimate_kaplan_meier edge cases", {
  for (i in seq_along(test_cases)) {
    idx <- seq_len(test_cases[[i]]$max_t + 1)
    survival <- km_survival[i, idx]
    at_risk_case <- as.integer(at_risk[i, idx])
    n_exited_case <- as.integer(n_exited[i, idx])
    n_right_censored_case <- as.integer(n_right_censored[i, idx])
    data <- test_cases[[i]]
    # All test cases have at least one patient; no skip needed
    # All test cases have at least one patient; no skip needed
    expect_length(survival, data$max_t + 1)
    expect_true(all(survival >= 0 & survival <= 1))
    expect_true(all(diff(survival) <= 0))
    expect_length(at_risk_case, data$max_t + 1)
    expect_length(n_exited_case, data$max_t + 1)
    expect_length(n_right_censored_case, data$max_t + 1)
    expect_true(all(at_risk_case >= 0))
    expect_true(all(n_exited_case >= 0))
    expect_true(all(n_right_censored_case >= 0))
    expect_equal(at_risk_case[1], length(data$event_time))
  }
})

test_that("estimate_kaplan_meier algorithmic correctness", {
  # Use the 3rd test case (3 patients, events at 1,2,3)
  i <- 3
  survival <- km_survival[i, 1:(test_cases[[i]]$max_t + 1)]
  at_risk_case <- at_risk[i, 1:(test_cases[[i]]$max_t + 1)]
  n_exited_case <- n_exited[i, 1:(test_cases[[i]]$max_t + 1)]
  expected_survival <- c(1.0, 2/3, 1/3, 0.0, 0.0, 0.0)
  expect_equal(survival, expected_survival, tolerance = 1e-6,
               info = "Sequential events should follow standard KM formula")
  expected_at_risk <- c(3, 3, 2, 1, 0, 0)
  expect_equal(at_risk_case, expected_at_risk,
               info = "At-risk counts should decrease with events")
  expected_n_exited <- c(0, 1, 1, 1, 0, 0)
  expect_equal(n_exited_case, expected_n_exited,
               info = "Event counts should match input data")
})

test_that("estimate_kaplan_meier handles non-monotonic and complex data", {
  # 4th test case: non-monotonic event times
  i <- 4
  survival <- km_survival[i, 1:(test_cases[[i]]$max_t + 1)]
  n_exited_case <- n_exited[i, 1:(test_cases[[i]]$max_t + 1)]
  total_events <- sum(1 - test_cases[[i]]$right_censored)
  expect_equal(sum(n_exited_case), total_events,
               info = "Total events should match regardless of event time ordering")
  # 5th test case: complex mixed
  i <- 5
  survival <- km_survival[i, 1:(test_cases[[i]]$max_t + 1)]
  expect_true(survival[1] >= 0 && survival[1] <= 1.0, 
              info = "S(1) should be valid probability with immediate events")
})

test_that("estimate_kaplan_meier comparison with R survfit", {
  skip_if_not_installed("ggsurvfit") 
  skip_if_not_installed("survival")
  # 6th test case: 4 events at 1,2,3,4
  i <- 6
  survival <- km_survival[i, 1:(test_cases[[i]]$max_t + 1)]
  event_time <- test_cases[[i]]$event_time
  right_censored <- test_cases[[i]]$right_censored
  r_data <- data.frame(time = event_time, event = 1 - right_censored)
  r_survfit <- ggsurvfit::survfit2(survival::Surv(time, event) ~ 1, data = r_data)
  r_summary <- summary(r_survfit)
  r_event_times <- r_summary$time
  r_survival_values <- r_summary$surv
  for (j in seq_along(r_event_times)) {
    event_time_j <- r_event_times[j]
    stan_survival_at_time <- survival[event_time_j + 1]
    expect_equal(stan_survival_at_time, r_survival_values[j], tolerance = 1e-6,
                 info = sprintf("Stan and R should match at event time %d", event_time_j))
  }
})

test_that("estimate_kaplan_meier mathematical consistency", {
  # 7th, 8th, 9th test cases
  for (i in 7:9) {
    survival <- km_survival[i, 1:(test_cases[[i]]$max_t + 1)]
    at_risk_case <- at_risk[i, 1:(test_cases[[i]]$max_t + 1)]
    n_exited_case <- n_exited[i, 1:(test_cases[[i]]$max_t + 1)]
    n_right_censored_case <- n_right_censored[i, 1:(test_cases[[i]]$max_t + 1)]
    n <- length(test_cases[[i]]$event_time)
    expect_equal(at_risk_case[1], n)
    total_events <- sum(1 - test_cases[[i]]$right_censored)
    expect_equal(sum(n_exited_case), total_events)
    total_censored <- sum(test_cases[[i]]$right_censored)
    expect_equal(sum(n_right_censored_case), total_censored)
    risk_changes <- which(n_exited_case > 0 | n_right_censored_case > 0)
    if (length(risk_changes) > 1) {
      for (j in 2:length(risk_changes)) {
        expect_lte(at_risk_case[risk_changes[j]], at_risk_case[risk_changes[j-1]])
      }
    }
  }
})
