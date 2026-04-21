library(testthat)
library(dplyr)
library(here)

source(here("r", "multistate.R"))

# =============================================================================
# Test helpers
# =============================================================================

# Build a single-row patient tibble (all fields needed by validate_ms_inputs)
make_patient <- function(
  original_right_censored  = FALSE,
  progression_before_death = FALSE,
  death                    = FALSE,
  death_week               = NA_integer_,
  patient_max_t            = 20L,
  potential_followup       = 80.0
) {
  tibble(
    original_right_censored  = original_right_censored,
    progression_before_death = progression_before_death,
    death                    = death,
    death_week               = death_week,
    patient_max_t            = patient_max_t,
    potential_followup       = potential_followup
  )
}

# One representative patient per pattern (6 rows total)
make_all_patterns <- function() {
  bind_rows(
    # A: admin_censored
    make_patient(original_right_censored = TRUE,  patient_max_t = 50L, potential_followup = 55.0),
    # B: true_dropout
    make_patient(original_right_censored = TRUE,  patient_max_t = 20L, potential_followup = 80.0),
    # C: progressed_alive
    make_patient(progression_before_death = TRUE, death = FALSE, patient_max_t = 30L, potential_followup = 80.0),
    # D: progressed_died
    make_patient(progression_before_death = TRUE, death = TRUE, death_week = 35L, patient_max_t = 30L, potential_followup = 80.0),
    # E: died_on_trial
    make_patient(death = TRUE, death_week = 23L, patient_max_t = 20L, potential_followup = 80.0),
    # F: died_off_trial
    make_patient(death = TRUE, death_week = 40L, patient_max_t = 20L, potential_followup = 80.0)
  )
}

# Build analysis_data suitable for derive_ms_fields (adds extra fields beyond patient data)
make_analysis_data <- function(patient_data, pfs = NULL, interval_censored = NULL) {
  n <- nrow(patient_data)
  patient_data |>
    mutate(
      # Default pfs: detection-adjusted event time (= patient_max_t for censored)
      pfs = if (!is.null(pfs)) pfs else as.integer(patient_max_t),
      interval_censored    = if (!is.null(interval_censored)) interval_censored else rep(0L, n),
      ms_prog_deterministic = rep(0L, n),
      ms_pattern = classify_ms_patients(pick(everything()))
    )
}


# =============================================================================
# validate_ms_inputs
# =============================================================================

test_that("validate_ms_inputs passes for all 6 valid patterns", {
  expect_true(validate_ms_inputs(make_all_patterns()))
})

test_that("I1: non-logical original_right_censored fails", {
  d <- make_patient()
  d$original_right_censored <- 0L
  expect_error(validate_ms_inputs(d), "I1")
})

test_that("I1: NA in original_right_censored fails", {
  d <- make_patient()
  d$original_right_censored <- NA
  expect_error(validate_ms_inputs(d), "I1")
})

test_that("I2: NA in progression_before_death fails", {
  d <- make_patient(death = TRUE, death_week = 20L)
  d$progression_before_death <- NA
  expect_error(validate_ms_inputs(d), "I2")
})

test_that("I3: NA in death fails", {
  d <- make_patient(original_right_censored = TRUE)
  d$death <- NA
  expect_error(validate_ms_inputs(d), "I3")
})

test_that("I4: NA death_week for dead patient fails", {
  d <- make_patient(death = TRUE, death_week = NA_integer_)
  expect_error(validate_ms_inputs(d), "I4")
})

test_that("I4: negative death_week fails", {
  d <- make_patient(death = TRUE, death_week = -1L)
  expect_error(validate_ms_inputs(d), "I4")
})

test_that("I5: patient_max_t = 0 fails", {
  d <- make_patient(original_right_censored = TRUE, patient_max_t = 0L,
                    potential_followup = 10.0)
  expect_error(validate_ms_inputs(d), "I5")
})

test_that("I6: potential_followup < patient_max_t fails", {
  d <- make_patient(original_right_censored = TRUE, patient_max_t = 50L,
                    potential_followup = 40.0)
  expect_error(validate_ms_inputs(d), "I6")
})

test_that("I7: uncensored with no progression and no death fails", {
  d <- make_patient(
    original_right_censored  = FALSE,
    progression_before_death = FALSE,
    death                    = FALSE
  )
  expect_error(validate_ms_inputs(d), "I7")
})

test_that("I8: right-censored AND progression_before_death warns (data anomaly)", {
  d <- make_patient(
    original_right_censored  = TRUE,
    progression_before_death = TRUE
  )
  # Downgraded to warning: these patients still route correctly as admin_censored/true_dropout
  expect_warning(validate_ms_inputs(d), "I8")
})


# =============================================================================
# classify_ms_patients
# =============================================================================

test_that("Pattern A: admin_censored when rc and near DCO", {
  d <- make_patient(original_right_censored = TRUE, patient_max_t = 50L, potential_followup = 54.0)
  expect_equal(as.character(classify_ms_patients(d)), "admin_censored")
})

test_that("Pattern B: true_dropout when rc and far from DCO", {
  d <- make_patient(original_right_censored = TRUE, patient_max_t = 20L, potential_followup = 80.0)
  expect_equal(as.character(classify_ms_patients(d)), "true_dropout")
})

test_that("Pattern C: progressed_alive when progressed and not dead", {
  d <- make_patient(progression_before_death = TRUE, death = FALSE, patient_max_t = 25L,
                    potential_followup = 80.0)
  expect_equal(as.character(classify_ms_patients(d)), "progressed_alive")
})

test_that("Pattern D: progressed_died when progressed and dead", {
  d <- make_patient(progression_before_death = TRUE, death = TRUE, death_week = 35L,
                    patient_max_t = 25L, potential_followup = 80.0)
  expect_equal(as.character(classify_ms_patients(d)), "progressed_died")
})

test_that("Pattern E: died_on_trial when death within buffer of last visit", {
  d <- make_patient(death = TRUE, death_week = 23L, patient_max_t = 20L, potential_followup = 80.0)
  expect_equal(as.character(classify_ms_patients(d)), "died_on_trial")
})

test_that("Pattern F: died_off_trial when death far after last visit", {
  d <- make_patient(death = TRUE, death_week = 40L, patient_max_t = 20L, potential_followup = 80.0)
  expect_equal(as.character(classify_ms_patients(d)), "died_off_trial")
})

test_that("All 6 patterns get assigned (no NAs) on a full dataset", {
  result <- classify_ms_patients(make_all_patterns())
  expect_false(anyNA(result))
  expect_setequal(
    as.character(result),
    c("admin_censored", "true_dropout", "progressed_alive",
      "progressed_died", "died_on_trial", "died_off_trial")
  )
})

test_that("Factor has exactly 6 levels in canonical order", {
  result <- classify_ms_patients(make_all_patterns())
  expect_equal(
    levels(result),
    c("admin_censored", "true_dropout", "progressed_alive",
      "progressed_died", "died_on_trial", "died_off_trial")
  )
})

test_that("Edge case: death_week == patient_max_t + 6 → died_on_trial (on boundary)", {
  # death_week - patient_max_t == 6 == buffer → should be on_trial (<=)
  d <- make_patient(death = TRUE, death_week = 26L, patient_max_t = 20L, potential_followup = 80.0)
  expect_equal(as.character(classify_ms_patients(d)), "died_on_trial")
})

test_that("Edge case: death_week == patient_max_t + 7 → died_off_trial (just over boundary)", {
  d <- make_patient(death = TRUE, death_week = 27L, patient_max_t = 20L, potential_followup = 80.0)
  expect_equal(as.character(classify_ms_patients(d)), "died_off_trial")
})


# =============================================================================
# validate_ms_routing
# =============================================================================

# Build minimal stan_fields for state-0 patients
make_stan_fields_state0 <- function(n, patient_max_t = 20L) {
  list(
    ms_final_state        = rep(0L, n),
    ms_time_01            = rep(patient_max_t, n),
    ms_time_02            = rep(patient_max_t, n),
    ms_time_12            = rep(0L, n),
    ms_censored_01        = rep(1L, n),
    ms_censored_02        = rep(1L, n),
    ms_censored_12        = rep(1L, n),
    ms_time_03            = rep(patient_max_t, n),
    ms_time_32            = rep(0L, n),
    ms_censored_32        = rep(1L, n),
    ms_prog_deterministic = rep(0L, n),
    ms_max_sojourn_t      = 11L,
    ms_max_sojourn_t_32   = 11L,
    ms_gp_grid_step       = 4L
  )
}

make_ad_row <- function(patient_max_t = 20L, ms_pattern = "admin_censored") {
  tibble(patient_max_t = patient_max_t, ms_pattern = ms_pattern)
}

test_that("validate_ms_routing passes for none mode all-zeros", {
  ad <- make_ad_row(ms_pattern = "admin_censored")
  sf <- make_stan_fields_state0(1)
  expect_true(validate_ms_routing(ad, sf, "none"))
})

test_that("validate_ms_routing passes for full mode state 3", {
  ad <- tibble(patient_max_t = 20L, ms_pattern = "true_dropout")
  sf <- make_stan_fields_state0(1)
  sf$ms_final_state <- 3L
  expect_true(validate_ms_routing(ad, sf, "full"))
})

test_that("R1: NA ms_pattern fails", {
  ad <- tibble(patient_max_t = 20L, ms_pattern = NA_character_)
  sf <- make_stan_fields_state0(1)
  expect_error(validate_ms_routing(ad, sf, "none"), "R1")
})

test_that("R2: state outside valid range for mode fails", {
  ad <- make_ad_row(ms_pattern = "admin_censored")
  sf <- make_stan_fields_state0(1)
  sf$ms_final_state <- 2L  # pfs only allows 0–1
  expect_error(validate_ms_routing(ad, sf, "pfs"), "R2")
})

test_that("R3: state 0 with ms_censored_01 != 1 fails", {
  ad <- make_ad_row(ms_pattern = "admin_censored")
  sf <- make_stan_fields_state0(1)
  sf$ms_censored_01 <- 0L
  expect_error(validate_ms_routing(ad, sf, "none"), "R3")
})

test_that("R4: state 1 with ms_time_01 = 0 fails", {
  ad <- make_ad_row(ms_pattern = "progressed_alive")
  sf <- make_stan_fields_state0(1)
  sf$ms_final_state <- 1L
  sf$ms_censored_01 <- 0L
  sf$ms_time_01 <- 0L
  expect_error(validate_ms_routing(ad, sf, "pfs"), "R4")
})

test_that("R8: state 3 in illness_death mode fails (caught by R2 or R8)", {
  # R2 fires first (state 3 > max_state 2), which is correct behaviour
  ad <- tibble(patient_max_t = 20L, ms_pattern = "true_dropout")
  sf <- make_stan_fields_state0(1)
  sf$ms_final_state <- 3L
  expect_error(validate_ms_routing(ad, sf, "illness_death"))
})

test_that("R9: progressed-then-died with ms_time_12 = 0 fails", {
  ad <- make_ad_row(ms_pattern = "progressed_died")
  sf <- make_stan_fields_state0(1)
  sf$ms_final_state <- 2L
  sf$ms_censored_01 <- 0L
  sf$ms_time_01 <- 15L
  sf$ms_time_12 <- 0L
  expect_error(validate_ms_routing(ad, sf, "full"), "R9")
})


# =============================================================================
# derive_ms_fields
# =============================================================================

test_that("derive_ms_fields: 'none' mode gives all state 0", {
  ad <- make_analysis_data(make_all_patterns())
  result <- derive_ms_fields(ad, "none")
  expect_true(all(result$ms_final_state == 0L))
  expect_true(all(result$ms_censored_01 == 1L))
})

test_that("derive_ms_fields: 'pfs' mode routes C/D/E → state 1, A/B/F → state 0", {
  patients <- make_all_patterns()
  ad <- make_analysis_data(patients, pfs = c(50L, 20L, 30L, 30L, 23L, 20L))
  result <- derive_ms_fields(ad, "pfs")
  patterns <- as.character(ad$ms_pattern)
  expect_equal(result$ms_final_state[patterns == "admin_censored"], 0L)
  expect_equal(result$ms_final_state[patterns == "true_dropout"], 0L)
  expect_equal(result$ms_final_state[patterns == "progressed_alive"], 1L)
  expect_equal(result$ms_final_state[patterns == "progressed_died"], 1L)
  expect_equal(result$ms_final_state[patterns == "died_on_trial"], 1L)
  expect_equal(result$ms_final_state[patterns == "died_off_trial"], 0L)
})

test_that("derive_ms_fields: 'illness_death' mode routes correctly", {
  patients <- make_all_patterns()
  ad <- make_analysis_data(patients, pfs = c(50L, 20L, 25L, 25L, 20L, 20L))
  result <- derive_ms_fields(ad, "illness_death")
  patterns <- as.character(ad$ms_pattern)
  expect_equal(result$ms_final_state[patterns == "admin_censored"], 0L)
  expect_equal(result$ms_final_state[patterns == "true_dropout"], 0L)
  expect_equal(result$ms_final_state[patterns == "progressed_alive"], 1L)
  expect_equal(result$ms_final_state[patterns == "progressed_died"], 2L)
  expect_equal(result$ms_final_state[patterns == "died_on_trial"], 2L)
  expect_equal(result$ms_final_state[patterns == "died_off_trial"], 0L)
})

test_that("derive_ms_fields: 'full' mode routes B/F → state 3", {
  patients <- make_all_patterns()
  ad <- make_analysis_data(patients, pfs = c(50L, 20L, 25L, 25L, 20L, 20L))
  result <- derive_ms_fields(ad, "full")
  patterns <- as.character(ad$ms_pattern)
  expect_equal(result$ms_final_state[patterns == "true_dropout"], 3L)
  expect_equal(result$ms_final_state[patterns == "died_off_trial"], 3L)
  expect_equal(result$ms_final_state[patterns == "admin_censored"], 0L)
})

test_that("derive_ms_fields: died_on_trial has ms_time_01 = 0 (sentinel) in illness_death", {
  d <- make_patient(death = TRUE, death_week = 23L, patient_max_t = 20L, potential_followup = 80.0)
  ad <- make_analysis_data(d, pfs = 20L)
  result <- derive_ms_fields(ad, "illness_death")
  expect_equal(result$ms_time_01, 0L)
  expect_equal(result$ms_final_state, 2L)
  expect_equal(result$ms_time_02, 23L)  # raw death_week, no adjustment
})

test_that("derive_ms_fields: progressed_died sojourn uses pmax(1, death_week - pfs)", {
  # death_week = 35, pfs (adjusted) = 30 → sojourn = 5
  d <- make_patient(progression_before_death = TRUE, death = TRUE, death_week = 35L,
                    patient_max_t = 30L, potential_followup = 80.0)
  ad <- make_analysis_data(d, pfs = 30L)
  result <- derive_ms_fields(ad, "full")
  expect_equal(result$ms_time_12, 5L)
  expect_equal(result$ms_censored_12, 0L)
})

test_that("derive_ms_fields: sojourn clamped to 1 when death same week as adjusted pfs", {
  # death_week = 25, pfs_adjusted = 26 → raw sojourn = -1 → pmax(1, -1) = 1
  d <- make_patient(progression_before_death = TRUE, death = TRUE, death_week = 25L,
                    patient_max_t = 20L, potential_followup = 80.0)
  ad <- make_analysis_data(d, pfs = 26L)
  result <- derive_ms_fields(ad, "full")
  expect_equal(result$ms_time_12, 1L)
})

test_that("derive_ms_fields: died_off_trial has correct ms_time_32", {
  # patient_max_t = 20, death_week = 40 → ms_time_32 = 20
  d <- make_patient(death = TRUE, death_week = 40L, patient_max_t = 20L, potential_followup = 80.0)
  ad <- make_analysis_data(d, pfs = 20L)
  result <- derive_ms_fields(ad, "full")
  expect_equal(result$ms_time_32, 20L)
  expect_equal(result$ms_censored_32, 0L)
  expect_equal(result$ms_final_state, 3L)
})

test_that("derive_ms_fields returns all required Stan field names", {
  ad <- make_analysis_data(make_all_patterns())
  result <- derive_ms_fields(ad, "full")
  expected_names <- c(
    "ms_final_state", "ms_time_01", "ms_time_02", "ms_time_12",
    "ms_censored_01", "ms_censored_02", "ms_censored_12",
    "ms_time_03", "ms_time_32", "ms_censored_32",
    "ms_prog_deterministic", "ms_max_sojourn_t", "ms_max_sojourn_t_32",
    "ms_gp_grid_step", "ms_os_event_12", "interval_censored"
  )
  expect_true(all(expected_names %in% names(result)))
})
