library(testthat)
library(here)

# Helper: run the share_dead_gp_shape constraint harness.
run_share_dead_gp_shape <- function(share_dead_gp_shape, ms_time_scale_12,
                                    enable_ms_12, enable_ms_02) {
  data <- list(
    share_dead_gp_shape = as.integer(share_dead_gp_shape),
    ms_time_scale_12    = as.integer(ms_time_scale_12),
    enable_ms_12        = as.integer(enable_ms_12),
    enable_ms_02        = as.integer(enable_ms_02)
  )
  test_stan_function(
    here("tests", "testthat", "stan", "test_share_dead_gp_shape_constraints.stan"),
    data
  )
}

# Runs a configuration expected to trigger Stan's fatal_error in transformed data.
# Returns TRUE if the chain failed as expected.
check_stan_fatal <- function(share_dead_gp_shape, ms_time_scale_12,
                             enable_ms_12, enable_ms_02) {
  warnings <- character()
  withCallingHandlers(
    tryCatch(
      run_share_dead_gp_shape(
        share_dead_gp_shape, ms_time_scale_12, enable_ms_12, enable_ms_02
      ),
      error = function(e) NULL
    ),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  any(grepl("Chain 1 finished unexpectedly|No chains finished successfully", warnings))
}

# =============================================================================
# Valid: share_dead_gp_shape OFF — no constraints apply
# =============================================================================

test_that("share_dead_gp_shape=0: any flag combo passes silently", {
  expect_no_error(run_share_dead_gp_shape(0, ms_time_scale_12 = 1, enable_ms_12 = 1, enable_ms_02 = 0))
  expect_no_error(run_share_dead_gp_shape(0, ms_time_scale_12 = 0, enable_ms_12 = 0, enable_ms_02 = 0))
})

# =============================================================================
# Valid: share_dead_gp_shape ON with all requirements met
# =============================================================================

test_that("share_dead_gp_shape=1 with ms_time_scale_12=0, enable_ms_12=1, enable_ms_02=1: passes", {
  expect_no_error(
    run_share_dead_gp_shape(1, ms_time_scale_12 = 0, enable_ms_12 = 1, enable_ms_02 = 1)
  )
})

# =============================================================================
# Invalid: ms_time_scale_12 != 0
# =============================================================================

test_that("share_dead_gp_shape=1 with ms_time_scale_12=1 (semi-Markov): chain failure", {
  expect_true(
    check_stan_fatal(1, ms_time_scale_12 = 1, enable_ms_12 = 1, enable_ms_02 = 1),
    info = "Expected Stan fatal_error when ms_time_scale_12 != 0"
  )
})

test_that("share_dead_gp_shape=1 with ms_time_scale_12=2 (extended): chain failure", {
  expect_true(
    check_stan_fatal(1, ms_time_scale_12 = 2, enable_ms_12 = 1, enable_ms_02 = 1),
    info = "Expected Stan fatal_error when ms_time_scale_12 != 0"
  )
})

# =============================================================================
# Invalid: need_12_t_gp = 0 (enable_ms_12 = 0)
# =============================================================================

test_that("share_dead_gp_shape=1 with enable_ms_12=0 (need_12_t_gp=0): chain failure", {
  expect_true(
    check_stan_fatal(1, ms_time_scale_12 = 0, enable_ms_12 = 0, enable_ms_02 = 1),
    info = "Expected Stan fatal_error when need_12_t_gp=0"
  )
})

# =============================================================================
# Invalid: enable_ms_02 = 0
# =============================================================================

test_that("share_dead_gp_shape=1 with enable_ms_02=0: chain failure", {
  expect_true(
    check_stan_fatal(1, ms_time_scale_12 = 0, enable_ms_12 = 1, enable_ms_02 = 0),
    info = "Expected Stan fatal_error when enable_ms_02=0"
  )
})
