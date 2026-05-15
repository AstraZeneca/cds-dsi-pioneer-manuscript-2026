# tests/testthat/test-stan-split-sd-cp-ncp-pos.R
# Tests the split_sd_cp_ncp_pos helper in stan/hierarchy.stanfunctions.
# Pattern mirrors test-stan-split-cp-ncp-pos.R (the 1D split helper).

library(testthat)
library(here)
library(posterior)

# R oracle: route (L, ℓ) pairs to _raw_ if mode in {1 FE, 2 RE}; route to
# _cp_ if mode == 4 RE_CP. Skip entries with mode == 0 NONE.
# Reject mode == 3 RE_GP.
r_split_sd_cp_ncp_pos <- function(sd_mode, n_groups_per_level) {
  n_levels <- nrow(sd_mode)
  raw_counts <- matrix(0L, n_levels, n_levels)
  cp_counts  <- matrix(0L, n_levels, n_levels)

  for (lv in seq_len(n_levels)) {
    for (sub_lv in seq_len(n_levels)) {
      mode <- sd_mode[lv, sub_lv]

      # Validate: ℓ >= L should be mode 0 (upper-triangular + diagonal should be NONE)
      if (sub_lv >= lv && mode != 0L) {
        stop("Invalid sd_mode: sub_lv >= lv must be NONE (0)")
      }

      # Validate: RE_GP not allowed
      if (mode == 3L) {
        stop("RE_GP (mode 3) not allowed for sub-hierarchy")
      }

      if (mode == 4L) {
        # RE_CP → cp bucket
        cp_counts[lv, sub_lv] <- n_groups_per_level[sub_lv]
      } else if (mode %in% c(1L, 2L)) {
        # FE or RE → raw bucket
        raw_counts[lv, sub_lv] <- n_groups_per_level[sub_lv]
      }
      # mode == 0 → skip
    }
  }

  # Compute row totals and cumulative position arrays
  raw_total <- sum(raw_counts)
  cp_total  <- sum(cp_counts)

  # Position arrays: [L, ℓ+1] indexing (ℓ goes 0..n_levels)
  # Positions are GLOBAL offsets into the flat parameter vectors — each row

  # starts where the previous row left off (not reset to 1).
  raw_pos <- matrix(1L, n_levels, n_levels + 1)
  cp_pos  <- matrix(1L, n_levels, n_levels + 1)

  raw_cursor <- 0L
  cp_cursor  <- 0L
  for (lv in seq_len(n_levels)) {
    raw_pos[lv, 1] <- raw_cursor + 1L
    cp_pos[lv, 1]  <- cp_cursor + 1L
    for (sub_lv in seq_len(n_levels)) {
      raw_pos[lv, sub_lv + 1] <- raw_pos[lv, sub_lv] + raw_counts[lv, sub_lv]
      cp_pos[lv, sub_lv + 1]  <- cp_pos[lv, sub_lv]  + cp_counts[lv, sub_lv]
    }
    raw_cursor <- raw_cursor + sum(raw_counts[lv, ])
    cp_cursor  <- cp_cursor  + sum(cp_counts[lv, ])
  }

  list(
    raw_total = raw_total,
    cp_total  = cp_total,
    raw_pos   = raw_pos,
    cp_pos    = cp_pos
  )
}

run_split_sd_test <- function(n_levels, n_groups_per_level, sd_mode) {
  stan_data <- list(
    n_levels           = n_levels,
    n_groups_per_level = as.array(as.integer(n_groups_per_level)),
    sd_mode            = sd_mode
  )

  fit <- test_stan_function(
    here("tests", "testthat", "stan", "test_split_sd_cp_ncp_pos.stan"),
    data = stan_data
  )

  d <- as_draws_df(fit$draws())
  exp <- r_split_sd_cp_ncp_pos(sd_mode, n_groups_per_level)

  expect_equal(as.integer(d$raw_total[1]), exp$raw_total)
  expect_equal(as.integer(d$cp_total[1]),  exp$cp_total)
  for (lv in seq_len(n_levels)) {
    for (sub_lv in seq_len(n_levels + 1)) {
      expect_equal(
        as.integer(d[[paste0("raw_pos[", lv, ",", sub_lv, "]")]][1]),
        exp$raw_pos[lv, sub_lv]
      )
      expect_equal(
        as.integer(d[[paste0("cp_pos[", lv, ",", sub_lv, "]")]][1]),
        exp$cp_pos[lv, sub_lv]
      )
    }
  }
}

# For error-path tests: bypass R oracle and call Stan directly to test validation
run_split_sd_test_stan_only <- function(n_levels, n_groups_per_level, sd_mode) {
  stan_data <- list(
    n_levels           = n_levels,
    n_groups_per_level = as.array(as.integer(n_groups_per_level)),
    sd_mode            = sd_mode
  )

  test_stan_function(
    here("tests", "testthat", "stan", "test_split_sd_cp_ncp_pos.stan"),
    data = stan_data
  )
}

test_that("all-zero mode matrix yields zero totals and identity-position arrays", {
  sd_mode <- matrix(0L, 3, 3)
  run_split_sd_test(3, c(2L, 3L, 10L), sd_mode)
})

test_that("single RE_CP entry goes into cp bucket only", {
  sd_mode <- matrix(0L, 3, 3)
  sd_mode[3, 2] <- 4L  # (L=patient, ℓ=arm) = RE_CP
  run_split_sd_test(3, c(2L, 3L, 10L), sd_mode)
})

test_that("mixed FE/RE/RE_CP partitions groups correctly across buckets", {
  sd_mode <- matrix(0L, 3, 3)
  sd_mode[3, 1] <- 2L  # (patient, trial) = RE  → raw bucket, 2 groups
  sd_mode[3, 2] <- 4L  # (patient, arm)   = RE_CP → cp bucket, 3 groups
  sd_mode[2, 1] <- 1L  # (arm, trial)     = FE   → raw bucket, 2 groups
  run_split_sd_test(3, c(2L, 3L, 10L), sd_mode)
})

test_that("multi-row positions are globally non-overlapping (issue #119)", {
  # When two rows both have active sub-levels, their positions must not overlap.

  # Row 2 (arm) has trial as sub-level RE; Row 3 (patient) has trial as sub-level RE.
  # Both contribute to the raw bucket — positions must be sequential, not both starting at 1.
  sd_mode <- matrix(0L, 3, 3)
  sd_mode[2, 1] <- 2L  # (arm, trial) = RE → raw bucket, 2 groups
  sd_mode[3, 1] <- 2L  # (patient, trial) = RE → raw bucket, 2 groups

  exp <- r_split_sd_cp_ncp_pos(sd_mode, c(2L, 3L, 10L))

  # Row 2 should start at 1, Row 3 should start at 3 (after Row 2's 2 groups)
  expect_equal(exp$raw_pos[2, 1], 1L)
  expect_equal(exp$raw_pos[3, 1], 3L)
  expect_equal(exp$raw_total, 4L)

  # Also verify via Stan
  run_split_sd_test(3, c(2L, 3L, 10L), sd_mode)
})

test_that("RE_GP (mode=3) triggers fatal_error", {
  sd_mode <- matrix(0L, 3, 3)
  sd_mode[3, 2] <- 3L  # RE_GP not allowed for sub-hierarchy

  # cmdstanr emits Stan's fatal_error as warnings. Capture all output (stdout and stderr)
  # plus warnings to verify that the Stan validation logic fires.
  warnings <- character()
  all_output <- capture.output(
    {
      result <- withCallingHandlers(
        tryCatch(
          run_split_sd_test_stan_only(3, c(2L, 3L, 10L), sd_mode),
          error = function(e) {
            # Suppress the "Fitting failed. Unable to print." error
            NULL
          }
        ),
        warning = function(w) {
          warnings <<- c(warnings, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      )
      # Don't return the result, just NULL
      invisible(NULL)
    },
    type = "output"  # Capture stdout instead of just messages
  )

  # Verify that cmdstanr detected the Stan failure
  expect_true(
    any(grepl("Chain 1 finished unexpectedly|No chains finished successfully", warnings)),
    info = paste("Expected cmdstanr failure warning, got:", paste(warnings, collapse = "; "))
  )

  # Verify that the Stan error message mentions RE_GP rejection
  combined_output <- paste(c(all_output, warnings), collapse = "\n")
  expect_match(
    combined_output,
    "RE_GP.*not allowed|mode.*3.*not allowed",
    ignore.case = TRUE,
    info = paste("Expected Stan error about RE_GP rejection, got:", combined_output)
  )
})

test_that("upper-triangular entries (ℓ >= L) that are non-zero trigger fatal_error", {
  sd_mode <- matrix(0L, 3, 3)
  sd_mode[1, 2] <- 4L  # L=trial cannot have sub-level ℓ=arm

  # cmdstanr emits Stan's fatal_error as warnings. Capture all output (stdout and stderr)
  # plus warnings to verify that the Stan validation logic fires.
  warnings <- character()
  all_output <- capture.output(
    {
      result <- withCallingHandlers(
        tryCatch(
          run_split_sd_test_stan_only(3, c(2L, 3L, 10L), sd_mode),
          error = function(e) {
            # Suppress the "Fitting failed. Unable to print." error
            NULL
          }
        ),
        warning = function(w) {
          warnings <<- c(warnings, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      )
      # Don't return the result, just NULL
      invisible(NULL)
    },
    type = "output"  # Capture stdout instead of just messages
  )

  # Verify that cmdstanr detected the Stan failure
  expect_true(
    any(grepl("Chain 1 finished unexpectedly|No chains finished successfully", warnings)),
    info = paste("Expected cmdstanr failure warning, got:", paste(warnings, collapse = "; "))
  )

  # Verify that the Stan error message mentions the upper-triangular constraint
  combined_output <- paste(c(all_output, warnings), collapse = "\n")
  expect_match(
    combined_output,
    "sub_lv.*>=.*L|sub-levels must be positionally before",
    ignore.case = TRUE,
    info = paste("Expected Stan error about upper-triangular constraint, got:", combined_output)
  )
})
