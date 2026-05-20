# tests/testthat/test-validate-sd-modes.R
library(testthat)

source(here::here("r/multi_level_hierarchy.R"))

test_that("empty tibble passes", {
  sd_modes <- tibble::tibble(location_level = character(), sub_level = character(), mode = character())
  expect_silent(validate_sd_modes(sd_modes, level_stack = c("trial", "arm", "patient")))
})

test_that("valid rows pass", {
  sd_modes <- tibble::tribble(
    ~location_level, ~sub_level, ~mode,
    "patient",       "arm",      "re_cp",
    "patient",       "trial",    "re",
    "arm",           "trial",    "fe"
  )
  expect_silent(validate_sd_modes(sd_modes, level_stack = c("trial", "arm", "patient")))
})

test_that("unknown location_level errors", {
  sd_modes <- tibble::tibble(location_level = "galaxy", sub_level = "arm", mode = "re_cp")
  expect_error(
    validate_sd_modes(sd_modes, level_stack = c("trial", "arm", "patient")),
    regexp = "galaxy|unknown|level_stack"
  )
})

test_that("unknown sub_level errors", {
  sd_modes <- tibble::tibble(location_level = "patient", sub_level = "region", mode = "re_cp")
  expect_error(
    validate_sd_modes(sd_modes, level_stack = c("trial", "arm", "patient")),
    regexp = "region|unknown|level_stack"
  )
})

test_that("sub_level not positionally before location_level errors", {
  sd_modes <- tibble::tibble(location_level = "trial", sub_level = "arm", mode = "re")
  expect_error(
    validate_sd_modes(sd_modes, level_stack = c("trial", "arm", "patient")),
    regexp = "position|before|nested"
  )
})

test_that("sub_level equal to location_level errors", {
  sd_modes <- tibble::tibble(location_level = "patient", sub_level = "patient", mode = "re_cp")
  expect_error(
    validate_sd_modes(sd_modes, level_stack = c("trial", "arm", "patient")),
    regexp = "position|before"
  )
})

test_that("disallowed mode 'gp' errors", {
  sd_modes <- tibble::tibble(location_level = "patient", sub_level = "arm", mode = "gp")
  expect_error(
    validate_sd_modes(sd_modes, level_stack = c("trial", "arm", "patient")),
    # Intentionally permissive: accept any of "mode", "gp", "allowed" in the message.
    # The implementer can word the error message naturally; the test asserts it names at least one of these.
    regexp = "mode|gp|allowed"
  )
})

test_that("typo mode 'random' errors", {
  sd_modes <- tibble::tibble(location_level = "patient", sub_level = "arm", mode = "random")
  expect_error(
    validate_sd_modes(sd_modes, level_stack = c("trial", "arm", "patient")),
    regexp = "mode|random|allowed"
  )
})

test_that("error message names the offending row", {
  sd_modes <- tibble::tibble(
    location_level = c("patient", "trial"),
    sub_level      = c("arm",     "arm"),
    mode           = c("re_cp",   "re")
  )
  expect_error(
    validate_sd_modes(sd_modes, level_stack = c("trial", "arm", "patient")),
    # Either the row number OR the (location, sub_level) pair must appear in the
    # error so the user can locate the offending row. Implementer can choose phrasing.
    regexp = "row 2|\\(trial, arm\\)"
  )
})
