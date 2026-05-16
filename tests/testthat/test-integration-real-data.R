# Integration tests using real data from Domino
# These tests are SKIPPED in CI (no access to /mnt/data)
# Run manually: Rscript -e 'Sys.setenv(RUN_INTEGRATION_TESTS="true"); testthat::test_file("tests/testthat/test-integration-real-data.R")'

library(testthat)
library(dplyr)
library(here)

# Skip these tests unless explicitly enabled
skip_if_no_integration <- function() {
  skip_if_not(
    Sys.getenv("RUN_INTEGRATION_TESTS") == "true",
    "Integration tests disabled. Set RUN_INTEGRATION_TESTS=true to run."
  )
  skip_if_not(
    dir.exists("/mnt/data"),
    "Real data not available (not on Domino)"
  )
}

test_that("sclc patient data loads correctly", {
  skip_if_no_integration()

  # Load real data from Domino
  patient_data <- read_csv(
    "/mnt/data/cooked/cooked_patient_data_200226.csv",
    show_col_types = FALSE
  )

  # Validate schema
  expect_true("pdl1_central" %in% names(patient_data))
  expect_gt(nrow(patient_data), 100)
})

test_that("pioneer multistate fields are valid on real data", {
  skip_if_no_integration()

  # Load from targets store
  tar_config_set(
    store = file.path(
      "/mnt/data/analysis-results",
      Sys.getenv("DOMINO_STARTING_USERNAME"),
      "pioneer/main/_targets"
    )
  )

  analysis_data <- tar_read(pioneer_analysis_data)

  # Validate multistate fields
  expect_true(all(analysis_data$ms_final_state %in% 0:3))
  expect_true(all(analysis_data$ms_time_01 >= 0))
})
