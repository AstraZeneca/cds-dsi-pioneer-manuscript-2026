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
