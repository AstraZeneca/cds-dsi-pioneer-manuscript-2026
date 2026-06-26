library(testthat)
library(targets)
library(tidyverse)

source(here::here("r/util.R"))
source(here::here("r/multistate.R"))
source(here::here("r/priors.R"))
source(here::here("r/prepare_analysis_data.R"))
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

# L1: get_lfo_cutoffs drops cutoffs where no target-trial patient has a post-cutoff visit
test_that("L1: get_lfo_cutoffs drops cutoffs with no target-trial future visits", {
  # Three patients: two target (lilly_cxcr4) enrolled late, one historical (amgen_darbe).
  # At the early cutoff (day 5), no target patient has a post-cutoff visit — only the
  # historical patient does. That cutoff must be dropped by n_target_future_observed > 0.
  #
  # target 1: calendar_day=100, visit at ady=1 (visit_calendar_day=100)
  # target 2: calendar_day=110, visit at ady=1 (visit_calendar_day=110)
  # historical: calendar_day=1,  visit at ady=1 and ady=10 (calendar days 1 and 10)
  # lfo_step=50 → cutoffs at days 100 and 150 (roughly)
  # At cutoff=5: historical has future visit (day 10 > 5) but no target patient has a
  # post-cutoff visit → should be dropped.
  # At cutoff=105: target 1 has visit at day 100 (≤ cutoff) but target 2 has visit at
  # day 110 (> cutoff) → kept.
  d <- tibble(
    trial = as_factor(c("lilly_cxcr4", "lilly_cxcr4", "amgen_darbe")),
    calendar_day = c(100L, 110L, 1L),
    patient_max_t = c(2L, 2L, 5L),
    visit_data = list(
      tibble(ady = 1L, week = 1L, visit_calendar_day = 100L),
      tibble(ady = 1L, week = 1L, visit_calendar_day = 110L),
      tibble(ady = c(1L, 10L), week = c(1L, 2L), visit_calendar_day = c(1L, 10L))
    )
  )
  result <- get_lfo_cutoffs(d, lfo_step = 50L, target_trial = "lilly_cxcr4")

  # All retained cutoffs must have at least one target patient with a future visit
  expect_true(
    all(result$n_target_future_observed > 0),
    label = "every retained cutoff must have >= 1 target patient with a post-cutoff visit"
  )

  # Specifically: no cutoff where only historical patients have future visits
  # (i.e., no cutoff where cutoff_calendar_day < 100, the first target enrollment)
  if (nrow(result) > 0) {
    expect_true(
      all(result$cutoff_calendar_day >= 100L),
      label = "no cutoff should be retained before any target patient has enrolled"
    )
  }
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

# H2: lfo_log_lik_rvar with exact=TRUE must not call loo::psis() and must not
# crash even when the log-ratio column contains NaN.
# The bug: in exact mode loo::psis() was called unconditionally, and a NaN in
# the input caused it to throw "NAs not allowed in input", aborting the whole
# lfo() recursion. The fix bypasses the psis_/k_/lwt_ mutates entirely in
# exact mode. Verified invariants: (1) no error thrown, (2) k = NA_real_.
test_that("H2: lfo_log_lik_rvar exact=TRUE does not crash on NaN log-ratio draws", {
  library(posterior)

  # Build the minimal tidy input that lfo_log_lik_rvar expects.
  # One cutoff (n=1), one forecast horizon (m=1), two patients.
  # One draw of patient 1 is NaN — this is what caused loo::psis() to crash
  # when the exact=TRUE path didn't skip the PSIS call.
  n_draws <- 10L

  draws_p1 <- c(NaN, rep(-1.5, n_draws - 1L))
  draws_p2 <- rep(-2.0, n_draws)
  mat <- cbind(draws_p1, draws_p2)

  log_lik_input <- tibble(
    n = 1L,
    m = 1L,
    patient_log_lik = list(mat)
  )

  # The fix: this must complete without error
  expect_no_error(
    result <- lfo_log_lik_rvar(log_lik_input, max_n = Inf, future_window = 1L, exact = TRUE)
  )

  # Must produce exactly one row
  expect_equal(nrow(result), 1L)

  # k must be NA_real_ in exact mode (PSIS was skipped, not computed)
  expect_true(
    is.na(result$k[[1]]),
    label = "k must be NA in exact mode — PSIS is bypassed"
  )

  # approx_E_patient_log_lik must exist (placeholder column) so that the
  # if_else() in clean_lfo_results resolves without a missing-column error
  expect_true(
    "approx_E_patient_log_lik" %in% names(result),
    label = "approx_E_patient_log_lik placeholder must exist for clean_lfo_results"
  )
})
