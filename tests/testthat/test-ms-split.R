library(testthat)
library(dplyr)
library(here)

source(here("r", "pioneer", "prepare_analysis_data.R"))
source(here("r", "multi_level_hierarchy.R"))

# =============================================================================
# ms_args_for_split
# =============================================================================

test_that("ms_args_for_split('all'): ms_split_level=0, empty target groups", {
  analysis_data <- tibble(trial = c("FPI-2265-202", "FLATIRON_PLUVICTO"))
  r <- ms_args_for_split("all", analysis_data)
  expect_equal(r$ms_split_level, 0L)
  expect_equal(r$ms_target_groups, integer(0))
})

test_that("ms_args_for_split('trial'): ms_split_level=1, RWD group excluded", {
  # FPI-2265-202 is the actual trial studyid; FLATIRON_PLUVICTO is RWD.
  # The trial split must include only the trial group ID.
  analysis_data <- tibble(trial = c("FPI-2265-202", "FPI-2265-202", "FLATIRON_PLUVICTO"))
  r <- ms_args_for_split("trial", analysis_data)
  expect_equal(r$ms_split_level, 1L)
  expect_length(r$ms_target_groups, 1L)
  rwd_id <- as.integer(as_factor(analysis_data$trial))[analysis_data$trial == "FLATIRON_PLUVICTO"][1]
  expect_false(rwd_id %in% r$ms_target_groups)
})

test_that("ms_args_for_split('trial'): target group IDs are valid factor levels", {
  analysis_data <- tibble(trial = c("FPI-2265-202", "FPI-2265-202", "FLATIRON_PLUVICTO"))
  r <- ms_args_for_split("trial", analysis_data)
  n_groups <- length(unique(analysis_data$trial))
  expect_true(all(r$ms_target_groups >= 1L & r$ms_target_groups <= n_groups))
})

test_that("ms_args_for_split: unknown split value errors", {
  expect_error(
    ms_args_for_split("rwd", tibble(trial = "FPI-2265-202")),
    "unknown ms_split"
  )
})

test_that("ms_args_for_split('trial'): single-trial dataset returns that one ID", {
  analysis_data <- tibble(trial = c("FPI-2265-202", "FPI-2265-202"))
  r <- ms_args_for_split("trial", analysis_data)
  expect_equal(r$ms_split_level, 1L)
  expect_length(r$ms_target_groups, 1L)
})

# =============================================================================
# compute_ms_patient_idx
# =============================================================================

# Build a minimal stan_data stub for compute_ms_patient_idx.
# n_patients total; forecast_patient_idx is the subset that are forecast.
# patient_level_groups is n_patients x n_levels, with level 1 = trial group.
make_ms_stan_data <- function(
    forecast_patient_idx,
    patient_trial_groups,   # integer vector, length = n_patients
    ms_split_level,
    ms_target_groups) {
  n_patients <- length(patient_trial_groups)
  list(
    forecast_patient_idx = forecast_patient_idx,
    patient_level_groups = matrix(patient_trial_groups, ncol = 1),
    ms_split_level       = ms_split_level,
    ms_target_groups     = ms_target_groups
  )
}

test_that("compute_ms_patient_idx: ms_split_level=0 returns all forecast patients", {
  sd <- make_ms_stan_data(
    forecast_patient_idx = 1:5,
    patient_trial_groups = c(1L, 1L, 1L, 2L, 2L),  # 3 trial, 2 RWD
    ms_split_level       = 0L,
    ms_target_groups     = integer(0)
  )
  expect_equal(compute_ms_patient_idx(sd), 1:5)
})

test_that("compute_ms_patient_idx: ms_split_level=1 keeps only target-group patients", {
  # Patients 1-3 are in group 1 (trial), 4-5 in group 2 (RWD).
  # Target = trial group 1.
  sd <- make_ms_stan_data(
    forecast_patient_idx = 1:5,
    patient_trial_groups = c(1L, 1L, 1L, 2L, 2L),
    ms_split_level       = 1L,
    ms_target_groups     = 1L
  )
  expect_equal(compute_ms_patient_idx(sd), 1:3)
})

test_that("compute_ms_patient_idx: target group matching all patients = same as split_level=0", {
  sd <- make_ms_stan_data(
    forecast_patient_idx = 1:4,
    patient_trial_groups = c(1L, 1L, 1L, 1L),
    ms_split_level       = 1L,
    ms_target_groups     = 1L
  )
  expect_equal(compute_ms_patient_idx(sd), 1:4)
})

test_that("compute_ms_patient_idx: no patients in target group returns empty integer", {
  # Target group 99 does not exist — no patients match.
  sd <- make_ms_stan_data(
    forecast_patient_idx = 1:3,
    patient_trial_groups = c(1L, 1L, 2L),
    ms_split_level       = 1L,
    ms_target_groups     = 99L
  )
  expect_equal(compute_ms_patient_idx(sd), integer(0))
})

test_that("compute_ms_patient_idx: subset of forecast_patient_idx is respected", {
  # All 6 patients exist, but only 3-6 are forecast.
  # Of those, 3-4 are trial, 5-6 are RWD.
  sd <- make_ms_stan_data(
    forecast_patient_idx = 3:6,
    patient_trial_groups = c(2L, 2L, 1L, 1L, 2L, 2L),
    ms_split_level       = 1L,
    ms_target_groups     = 1L
  )
  expect_equal(compute_ms_patient_idx(sd), c(3L, 4L))
})
