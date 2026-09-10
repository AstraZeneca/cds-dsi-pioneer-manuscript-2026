library(testthat)
library(tibble)
library(here)

source(here("r", "pioneer", "prepare_analysis_data.R"))
source(here("r", "multi_level_hierarchy.R"))

# Helpers -------------------------------------------------------------------

# level_stack and n_groups for the typical pioneer 3-level hierarchy:
#   trial (2 groups) -> arm (4 groups) -> patient (100 patients)
STACK_3  <- c("trial", "arm", "patient")
NGROUPS_3 <- c(trial = 2L, arm = 4L, patient = 100L)

serialize <- function(modes, stack = STACK_3, ngroups = NGROUPS_3) {
  serialize_sd_subhierarchy_modes(modes, stack, ngroups)
}

# =============================================================================
# NULL / empty input → all-NONE matrix
# =============================================================================

test_that("NULL input: mode matrix is all zeros", {
  r <- serialize(NULL)
  expect_true(all(r$enable_sd_level_intercept_mode_tr == 0L))
  expect_equal(dim(r$enable_sd_level_intercept_mode_tr), c(3L, 3L))
})

test_that("NULL input: raw and CP counts are zero", {
  r <- serialize(NULL)
  expect_equal(r$n_raw_groups_tr_log_sd_intercept, 0L)
  expect_equal(r$n_cp_groups_tr_log_sd_intercept,  0L)
})

test_that("empty tibble behaves identically to NULL", {
  empty <- tibble(location_level = character(), sub_level = character(), mode = character())
  r_null  <- serialize(NULL)
  r_empty <- serialize(empty)
  expect_equal(r_null, r_empty)
})

# =============================================================================
# Single-row inputs — matrix placement
# =============================================================================

test_that("single re_cp row: correct matrix cell set to 4", {
  modes <- tibble(location_level = "patient", sub_level = "arm", mode = "re_cp")
  r <- serialize(modes)
  # patient = level 3, arm = level 2 → cell [3, 2] should be 4 (RE_CP)
  expect_equal(r$enable_sd_level_intercept_mode_tr[3, 2], 4L)
  # All other cells zero
  mat <- r$enable_sd_level_intercept_mode_tr
  mat[3, 2] <- 0L
  expect_true(all(mat == 0L))
})

test_that("single re row: correct matrix cell set to 2", {
  modes <- tibble(location_level = "arm", sub_level = "trial", mode = "re")
  r <- serialize(modes)
  # arm = level 2, trial = level 1 → cell [2, 1] = 2
  expect_equal(r$enable_sd_level_intercept_mode_tr[2, 1], 2L)
})

test_that("single fe row: correct matrix cell set to 1", {
  modes <- tibble(location_level = "patient", sub_level = "trial", mode = "fe")
  r <- serialize(modes)
  expect_equal(r$enable_sd_level_intercept_mode_tr[3, 1], 1L)
})

# =============================================================================
# Bucket counts — must mirror Stan's split_sd_cp_ncp_pos logic
# =============================================================================

test_that("re mode at patient/arm: raw count = n_arm_groups, cp count = 0", {
  modes <- tibble(location_level = "patient", sub_level = "arm", mode = "re")
  r <- serialize(modes)
  # RE → raw bucket; arm has 4 groups
  expect_equal(r$n_raw_groups_tr_log_sd_intercept, 4L)
  expect_equal(r$n_cp_groups_tr_log_sd_intercept,  0L)
})

test_that("re_cp mode at patient/arm: cp count = n_arm_groups, raw count = 0", {
  modes <- tibble(location_level = "patient", sub_level = "arm", mode = "re_cp")
  r <- serialize(modes)
  expect_equal(r$n_raw_groups_tr_log_sd_intercept, 0L)
  expect_equal(r$n_cp_groups_tr_log_sd_intercept,  4L)
})

test_that("fe mode at patient/arm: raw count = n_arm_groups (FE routes to raw bucket)", {
  modes <- tibble(location_level = "patient", sub_level = "arm", mode = "fe")
  r <- serialize(modes)
  expect_equal(r$n_raw_groups_tr_log_sd_intercept, 4L)
  expect_equal(r$n_cp_groups_tr_log_sd_intercept,  0L)
})

test_that("mixed re + re_cp rows: raw and cp counts accumulate correctly", {
  modes <- tribble(
    ~location_level, ~sub_level, ~mode,
    "patient",       "arm",      "re",      # raw += 4 (n_arm_groups)
    "arm",           "trial",    "re_cp"    # cp  += 2 (n_trial_groups)
  )
  r <- serialize(modes)
  expect_equal(r$n_raw_groups_tr_log_sd_intercept, 4L)
  expect_equal(r$n_cp_groups_tr_log_sd_intercept,  2L)
})

test_that("none mode: contributes nothing to raw or cp counts", {
  modes <- tibble(location_level = "patient", sub_level = "arm", mode = "none")
  r <- serialize(modes)
  expect_equal(r$n_raw_groups_tr_log_sd_intercept, 0L)
  expect_equal(r$n_cp_groups_tr_log_sd_intercept,  0L)
  # Cell should still be set to 0
  expect_equal(r$enable_sd_level_intercept_mode_tr[3, 2], 0L)
})

test_that("bucket counts match between serialize and independent R oracle", {
  # Oracle: sum n_groups[sub_range] * (mode %in% raw_modes or == 4) for each L
  oracle_counts <- function(modes_tbl, stack, ngroups) {
    mode_codes <- c(none = 0L, fe = 1L, re = 2L, re_gp = 3L, re_cp = 4L)
    n <- length(stack)
    mat <- matrix(0L, n, n)
    for (i in seq_len(nrow(modes_tbl))) {
      L   <- which(stack == modes_tbl$location_level[i])
      sub <- which(stack == modes_tbl$sub_level[i])
      mat[L, sub] <- mode_codes[[modes_tbl$mode[i]]]
    }
    n_raw <- 0L
    n_cp  <- 0L
    for (L in seq_len(n)) {
      if (L == 1) next
      for (s in seq_len(L - 1)) {
        if (mat[L, s] %in% c(1L, 2L)) n_raw <- n_raw + ngroups[s]
        if (mat[L, s] == 4L)           n_cp  <- n_cp  + ngroups[s]
      }
    }
    list(n_raw = unname(n_raw), n_cp = unname(n_cp))
  }

  modes <- tribble(
    ~location_level, ~sub_level, ~mode,
    "patient",       "arm",      "re_cp",
    "patient",       "trial",    "re",
    "arm",           "trial",    "fe"
  )
  r   <- serialize(modes)
  exp <- oracle_counts(modes, STACK_3, NGROUPS_3)
  expect_equal(r$n_raw_groups_tr_log_sd_intercept, exp$n_raw)
  expect_equal(r$n_cp_groups_tr_log_sd_intercept,  exp$n_cp)
})

# =============================================================================
# Invalid inputs are caught by validate_sd_modes
# =============================================================================

test_that("unknown level in modes errors", {
  modes <- tibble(location_level = "galaxy", sub_level = "arm", mode = "re")
  expect_error(serialize(modes), regexp = "galaxy|unknown|level_stack")
})

test_that("sub_level not before location_level errors", {
  modes <- tibble(location_level = "trial", sub_level = "arm", mode = "re")
  expect_error(serialize(modes), regexp = "position|before|nested")
})
