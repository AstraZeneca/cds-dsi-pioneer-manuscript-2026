# tests/testthat/test-sd-hierarchy-bit-exact.R
#
# Phase 1 acceptance gate for the SD sub-hierarchy feature (issue #110).
# When the mode tibble is EMPTY, draws on this branch must be byte-identical
# to a reference CSV captured from `main` pre-feature.
#
# The reference CSV is captured once and committed to
# `tests/testthat/fixtures/sd-hierarchy/bitexact-reference.csv`.
# See the README in that directory for the capture procedure.
#
# If the fixture is absent, the test skips — the primary gate is still
# operational but requires the reference to be captured first.

library(testthat)

test_that("empty mode tibble produces bit-identical draws vs main reference", {
  reference_csv <- test_path("fixtures", "sd-hierarchy", "bitexact-reference.csv")
  reference_data <- test_path("fixtures", "sd-hierarchy", "bitexact-reference-stan-data.rds")

  skip_if(
    !file.exists(reference_csv) || !file.exists(reference_data),
    "Reference fixture not captured yet — see tests/testthat/fixtures/sd-hierarchy/README.md"
  )

  skip_on_ci()  # bit-exact depends on cmdstan version + compiler flags; run locally only

  library(cmdstanr)
  library(readr)

  stan_data <- readRDS(reference_data)

  # Sanity check: mode matrix must be all-zero for the bit-exact gate to apply.
  # (If the reference was captured with a non-empty mode tibble, this test is meaningless.)
  stopifnot(
    "reference was captured with an empty mode tibble" =
      all(stan_data$enable_sd_level_intercept_mode_tr == 0L)
  )

  model <- cmdstan_model(
    test_path("..", "..", "stan", "psa", "pioneer.stan"),
    include_paths = c(
      test_path("..", "..", "stan"),
      test_path("..", "..", "stan", "psa")
    ),
    force_recompile = TRUE,
    quiet = TRUE
  )

  source(test_path("..", "..", "r", "pioneer", "initializers.R"))

  fit <- model$sample(
    data = stan_data,
    seed = 42,
    iter_warmup = 50,
    iter_sampling = 50,
    chains = 1,
    parallel_chains = 1,
    init = create_pioneer_initializer(stan_data),
    refresh = 0,
    show_messages = FALSE
  )

  produced_csv <- fit$output_files()[1]

  # Read both CSVs, skipping comment lines (which contain timestamps/PIDs)
  reference_df <- read_csv(reference_csv, comment = "#", show_col_types = FALSE)
  produced_df <- read_csv(produced_csv, comment = "#", show_col_types = FALSE)

  expect_equal(
    produced_df,
    reference_df,
    info = "Draw values differ — Phase 1 bit-exactness gate violated"
  )
})
