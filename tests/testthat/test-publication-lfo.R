library(testthat)
library(targets)
library(tidyverse)

source(here::here("r/util.R"))
source(here::here("r/multistate.R"))
source(here::here("r/sclc/priors.R"))
source(here::here("r/sclc/prepare_analysis_data.R"))
source(here::here("r/publication/prepare_analysis_data.R"))

make_minimal_pub_data <- function() {
  target_patient <- tibble(
    studyid = "S1", usubjid = "P1", trial = "lilly_cxcr4",
    calendar_day = 100L, calendar_week = 15L,
    pfs = 10L, death = FALSE, death_week = NA_integer_,
    progression_before_death = NA, right_censored = TRUE,
    interval_censored = 0L, age = 60, sex = 0L, ecog = 1L,
    hgb = 120, ldh_log = 5, albumin = 40,
    race = "White", stage = "IV", prev_lines = 1L, smoker = TRUE,
    patient_min_t = 1L, patient_max_t = 10L, patient_t_width = 10L
  )
  historical_patient <- target_patient |>
    mutate(usubjid = "P2", trial = "amgen_darbe", calendar_day = 1L, calendar_week = 1L)

  target_visit <- tibble(
    studyid = "S1", usubjid = "P1",
    visitnum = 1L, ady = 1L, week = 1L, mmsumdiam = 50, response = "SD"
  )
  historical_visit <- target_visit |>
    mutate(usubjid = "P2")

  list(
    target_patient = target_patient,
    historical_patient = historical_patient,
    target_visit = target_visit |> determine_visit_data_response(),
    historical_visit = historical_visit |> determine_visit_data_response()
  )
}

test_that("prepare_publication_analysis_data adds visit_calendar_day to nested visit data", {
  d <- make_minimal_pub_data()
  result <- prepare_publication_analysis_data(
    d$target_patient, d$historical_patient,
    d$target_visit, d$historical_visit
  )
  visit_cols <- result |> pull(visit_data) |> purrr::map(colnames) |> purrr::reduce(intersect)
  expect_true("visit_calendar_day" %in% visit_cols)
})

test_that("visit_calendar_day equals calendar_day + ady - 1", {
  d <- make_minimal_pub_data()
  result <- prepare_publication_analysis_data(
    d$target_patient, d$historical_patient,
    d$target_visit, d$historical_visit
  )
  # Target patient: calendar_day=100, ady=1 => visit_calendar_day = 100 + 1 - 1 = 100
  target_row <- result |> filter(fct_match(trial, "lilly_cxcr4"))
  expect_equal(
    target_row$visit_data[[1]]$visit_calendar_day,
    100L  # calendar_day=100, ady=1 => 100 + 1 - 1 = 100
  )
})

source(here::here("r/accuracy.R"))
source(here::here("r/sclc/accuracy.R"))

test_that("get_lfo_cutoffs works for a non-sclc target_trial", {
  # Two patients: target trial lilly_cxcr4 (calendar_day 100, 130),
  # historical amgen_darbe (calendar_day 1)
  d <- tibble(
    trial = as_factor(c("lilly_cxcr4", "lilly_cxcr4", "amgen_darbe")),
    calendar_day = c(100L, 130L, 1L),
    patient_max_t = c(10L, 5L, 20L),
    visit_data = list(
      tibble(ady = c(1L, 7L), week = c(1L, 2L),
             visit_calendar_day = 100L + c(1L, 7L) - 1L),
      tibble(ady = c(1L, 7L), week = c(1L, 2L),
             visit_calendar_day = 130L + c(1L, 7L) - 1L),
      tibble(ady = c(1L, 7L), week = c(1L, 2L),
             visit_calendar_day = 1L + c(1L, 7L) - 1L)
    )
  )
  result <- get_lfo_cutoffs(d, lfo_step = 30L,
                             target_trial = "lilly_cxcr4")
  expect_s3_class(result, "data.frame")
  expect_true(all(c("n", "cutoff_date", "cutoff_calendar_day") %in% names(result)))
  expect_true(nrow(result) >= 1L)
  # All cutoff_calendar_days should be within the target trial's calendar_day range
  expect_true(all(result$cutoff_calendar_day >= min(d$calendar_day[d$trial == "lilly_cxcr4"])))
})

test_that("get_lfo_cutoffs default target_trial='sclc' still works", {
  # Sclc data has trtsdt — origin recovery branch must fire
  d <- tibble(
    trial = as_factor(c("sclc", "historical")),
    trtsdt = as.Date(c("2020-01-01", "2019-06-01")),
    calendar_day = c(215L, 1L),
    patient_max_t = c(10L, 20L),
    visit_data = list(
      tibble(ady = c(1L, 7L), week = c(1L, 2L),
             visit_calendar_day = 215L + c(1L, 7L) - 1L),
      tibble(ady = c(1L, 7L), week = c(1L, 2L),
             visit_calendar_day = 1L + c(1L, 7L) - 1L)
    )
  )
  result <- get_lfo_cutoffs(d, lfo_step = 30L)
  expect_s3_class(result, "data.frame")
  expect_true(nrow(result) >= 1L)
  # cutoff_date should be a real Date
  expect_s3_class(result$cutoff_date, "Date")
})
