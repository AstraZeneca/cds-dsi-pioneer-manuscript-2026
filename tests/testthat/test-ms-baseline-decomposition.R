# Tests for decompose_ms_level_baseline_hazard() and the per-slot init
# dimensioning that consumes the decomposed multistate baseline-hazard schema.
#
# The legacy single per-level flag `enable_ms_level_baseline_hazard` (0-4) was
# replaced by three orthogonal per-transition x per-level arrays:
#   enable_ms_level_gp            (0/1)
#   ms_level_intercept_mode       (0-3: none/FE/RE-NCP/RE-CP)
#   ms_level_intercept_corr_group (>=0)
# Phase 1 (decomposition only) must reproduce the legacy single-flag model
# bit-identically: corr_group all-zero, and every slot reconstructs exactly the
# legacy mode. These tests pin that contract on the R side.

library(testthat)
library(dplyr)

source(here::here("r", "multistate.R"))

# Mirror of the Stan helper ms_reconstruct_legacy_mode (slot active).
reconstruct_legacy_mode_r <- function(gp_row, mode_row) {
  dplyr::case_when(
    gp_row == 1L   ~ 3L,
    mode_row == 0L ~ 0L,
    mode_row == 1L ~ 1L,
    mode_row == 2L ~ 2L,
    TRUE           ~ 4L
  )
}

# ---------------------------------------------------------------------------
# DIMENSIONING
# ---------------------------------------------------------------------------

test_that("decompose returns 3 arrays sized N_TRANS x n_levels", {
  for (legacy in list(c(trial = 3L, patient = 0L), c(a = 0L), c(x = 4L, y = 2L, z = 1L))) {
    d <- decompose_ms_level_baseline_hazard(legacy)
    n_levels <- length(legacy)
    expect_equal(dim(d$enable_ms_level_gp), c(MS_N_INTERCEPT_SLOTS, n_levels))
    expect_equal(dim(d$ms_level_intercept_mode), c(MS_N_INTERCEPT_SLOTS, n_levels))
    expect_equal(dim(d$ms_level_intercept_corr_group), c(MS_N_INTERCEPT_SLOTS, n_levels))
    expect_type(d$enable_ms_level_gp, "integer")
    expect_type(d$ms_level_intercept_mode, "integer")
    expect_type(d$ms_level_intercept_corr_group, "integer")
  }
})

test_that("decompose keeps corr_group all-zero (Phase 1 = independent singletons)", {
  d <- decompose_ms_level_baseline_hazard(c(trial = 3L, patient = 2L))
  expect_true(all(d$ms_level_intercept_corr_group == 0L))
})

# ---------------------------------------------------------------------------
# DEGENERACY LINCHPIN — every slot reconstructs the legacy mode bit-identically
# ---------------------------------------------------------------------------

test_that("every slot round-trips to the exact legacy mode vector", {
  legacy_configs <- list(
    c(trial = 3L, patient = 0L),  # publication
    c(trial = 2L, patient = 0L),  # sclc joint / standalone
    c(0L, 0L),                    # all-off
    c(1L, 4L),                    # FE + RE-CP
    c(4L, 3L, 0L)                 # RE-CP + RE-GP + none, 3 levels
  )
  for (legacy in legacy_configs) {
    legacy <- as.integer(unname(legacy))
    d <- decompose_ms_level_baseline_hazard(legacy)
    for (k in seq_len(MS_N_INTERCEPT_SLOTS)) {
      rec <- reconstruct_legacy_mode_r(
        d$enable_ms_level_gp[k, ], d$ms_level_intercept_mode[k, ]
      )
      expect_equal(rec, legacy,
        info = paste("slot", k, "legacy", paste(legacy, collapse = ",")))
    }
  }
})

test_that("legacy GP (mode 3) maps to gp=1 + RE-NCP intercept", {
  d <- decompose_ms_level_baseline_hazard(c(trial = 3L, patient = 0L))
  # level 1 (trial) legacy 3 => gp 1, intercept mode 2 (RE-NCP)
  expect_true(all(d$enable_ms_level_gp[, 1] == 1L))
  expect_true(all(d$ms_level_intercept_mode[, 1] == 2L))
  # level 2 (patient) legacy 0 => gp 0, intercept mode 0
  expect_true(all(d$enable_ms_level_gp[, 2] == 0L))
  expect_true(all(d$ms_level_intercept_mode[, 2] == 0L))
})

test_that("legacy RE-CP (mode 4) maps to gp=0 + intercept mode 3", {
  d <- decompose_ms_level_baseline_hazard(c(trial = 4L, patient = 0L))
  expect_true(all(d$enable_ms_level_gp[, 1] == 0L))
  expect_true(all(d$ms_level_intercept_mode[, 1] == 3L))
})

# ---------------------------------------------------------------------------
# VALIDATION
# ---------------------------------------------------------------------------

test_that("decompose rejects out-of-range legacy modes", {
  expect_error(decompose_ms_level_baseline_hazard(c(trial = 5L)), "0-4")
  expect_error(decompose_ms_level_baseline_hazard(c(trial = -1L)), "0-4")
})

# ---------------------------------------------------------------------------
# INIT DIMENSIONING — per-slot group counts from the decomposed arrays
# ---------------------------------------------------------------------------

test_that("ms_init_values_fixed sizes the 0->1 RE intercept per-slot (trial RE)", {
  source(here::here("r", "initializers_fixed.R"))
  # trial-level RE (mode 2), patient none; 2 levels with 3 trials, 10 patients.
  d <- decompose_ms_level_baseline_hazard(c(trial = 2L, patient = 0L))
  env <- c(list(
    n_time_varying_covar = 0L, n_time_invariant_covar = 0L,
    n_levels = 2L, n_groups_per_level = c(trial = 3L, patient = 10L),
    n_patients = 10L,
    enable_ms_01 = 1L, enable_ms_02 = 0L, enable_ms_03 = 0L,
    enable_ms_12 = 0L, enable_ms_32 = 0L, ms_time_scale_12 = 1L,
    enable_ms_level_cov = c(trial = 0L, patient = 0L),
    enable_ms_pop_time_varying_cov = 0L, enable_ms_pop_time_invariant_cov = 0L,
    enable_ms_visit_gated_01 = 0L, enable_ms_visit_gated_latent_01 = 0L,
    enable_ms_02_time_varying_cov = 0L, enable_ms_03_time_invariant_cov = 0L,
    enable_ms_03_time_varying_cov = 0L, enable_ms_32_time_invariant_cov = 0L,
    enable_ms_12_entry_covar = 0L, enable_ms_32_entry_covar = 0L,
    enable_ms_baseline_trend_01 = 0L,
    share_dead_gp_shape = 0L,
    max_all_t = 50L, ms_max_sojourn_t = 50L, ms_max_sojourn_t_32 = 50L,
    ms_gp_grid_step = 4L, enable_states_full_grid = 0L
  ), d)
  init <- ms_init_values_fixed(env)
  # RE intercept at trial level => 3 enabled groups; no GP eta.
  expect_length(init$raw_log_lambda_gp_01_level_intercept, 3L)
  expect_length(init$log_lambda_gp_01_level_intercept_sd, 2L)  # any_re => n_levels
})

# ---------------------------------------------------------------------------
# CORRELATED INTERCEPT BLOCKS (Phase 2)
# ---------------------------------------------------------------------------

test_that("ms_corr_blocks: all-zero corr_group => no blocks (degeneracy)", {
  cg <- matrix(0L, nrow = MS_N_INTERCEPT_SLOTS, ncol = 2L)
  active <- c(1L, 0L, 1L, 0L, 0L, 0L)  # 01 + 03 active
  corr <- ms_corr_blocks(cg, active, c(trial = 3L, patient = 10L))
  expect_equal(corr$n_blocks, 0L)
  expect_equal(corr$dim, 0L)
  expect_equal(corr$n_groups, 0L)
  expect_length(ms_corr_block_inits(corr), 0L)
})

test_that("ms_corr_blocks: 01+03 patient-level group 1 => one d=2 block", {
  cg <- matrix(0L, nrow = MS_N_INTERCEPT_SLOTS, ncol = 2L)
  cg[1, 2] <- 1L  # slot 01, patient level
  cg[3, 2] <- 1L  # slot 03, patient level
  active <- c(1L, 0L, 1L, 0L, 0L, 0L)
  corr <- ms_corr_blocks(cg, active, c(trial = 3L, patient = 10L))
  expect_equal(corr$n_blocks, 1L)
  expect_equal(corr$dim, 2L)
  expect_equal(corr$n_groups, 10L)        # patient-level group count
  expect_equal(corr$blocks[[1]]$level, 2L)
  expect_equal(corr$blocks[[1]]$member_slots, c(1L, 3L))
})

test_that("ms_corr_blocks: singleton (one active member) collapses, no block", {
  cg <- matrix(0L, nrow = MS_N_INTERCEPT_SLOTS, ncol = 2L)
  cg[1, 2] <- 1L  # only slot 01 carries the code
  cg[3, 2] <- 1L  # slot 03 carries it too, but 03 is inactive below
  active <- c(1L, 0L, 0L, 0L, 0L, 0L)  # only 01 active => group 1 has 1 member
  corr <- ms_corr_blocks(cg, active, c(trial = 3L, patient = 10L))
  expect_equal(corr$n_blocks, 0L)
})

test_that("ms_corr_block_inits: identity Cholesky + sized z, length n_blocks", {
  cg <- matrix(0L, nrow = MS_N_INTERCEPT_SLOTS, ncol = 2L)
  cg[1, 2] <- 1L; cg[3, 2] <- 1L
  active <- c(1L, 0L, 1L, 0L, 0L, 0L)
  corr <- ms_corr_blocks(cg, active, c(trial = 3L, patient = 10L))
  set.seed(1)
  inits <- ms_corr_block_inits(corr, z_sd = 0.3)
  expect_length(inits$L_ms_intercept_corr, 1L)
  expect_length(inits$z_ms_intercept, 1L)
  expect_equal(inits$L_ms_intercept_corr[[1]], diag(2))         # identity
  expect_equal(dim(inits$z_ms_intercept[[1]]), c(2L, 10L))      # d x n_groups
})

test_that("ms_corr_blocks: mismatched dimensions across blocks errors", {
  # block A at patient level: 01+03 (d=2); block B at trial level: 02+12_s+12_t
  # (d=3) under a different group code -> non-uniform dimension -> stop().
  cg <- matrix(0L, nrow = MS_N_INTERCEPT_SLOTS, ncol = 2L)
  cg[1, 2] <- 1L; cg[3, 2] <- 1L           # patient block, d=2
  cg[2, 1] <- 1L; cg[4, 1] <- 1L; cg[5, 1] <- 1L  # trial block, d=3
  active <- c(1L, 1L, 1L, 1L, 1L, 0L)
  expect_error(
    ms_corr_blocks(cg, active, c(trial = 3L, patient = 10L)),
    "same member"
  )
})

test_that("ms_init_values_fixed: corr_group stays absent (no Phase-2 params yet)", {
  source(here::here("r", "initializers_fixed.R"))
  d <- decompose_ms_level_baseline_hazard(c(trial = 3L, patient = 0L))
  env <- c(list(
    n_time_varying_covar = 0L, n_time_invariant_covar = 0L,
    n_levels = 2L, n_groups_per_level = c(trial = 3L, patient = 10L),
    n_patients = 10L,
    enable_ms_01 = 1L, enable_ms_02 = 0L, enable_ms_03 = 0L,
    enable_ms_12 = 0L, enable_ms_32 = 0L, ms_time_scale_12 = 1L,
    enable_ms_level_cov = c(trial = 0L, patient = 0L),
    enable_ms_pop_time_varying_cov = 0L, enable_ms_pop_time_invariant_cov = 0L,
    enable_ms_visit_gated_01 = 0L, enable_ms_visit_gated_latent_01 = 0L,
    enable_ms_02_time_varying_cov = 0L, enable_ms_03_time_invariant_cov = 0L,
    enable_ms_03_time_varying_cov = 0L, enable_ms_32_time_invariant_cov = 0L,
    enable_ms_12_entry_covar = 0L, enable_ms_32_entry_covar = 0L,
    enable_ms_baseline_trend_01 = 0L,
    share_dead_gp_shape = 0L,
    max_all_t = 50L, ms_max_sojourn_t = 50L, ms_max_sojourn_t_32 = 50L,
    ms_gp_grid_step = 4L, enable_states_full_grid = 0L
  ), d)
  init <- ms_init_values_fixed(env)
  expect_null(init$L_ms_intercept_corr)
  expect_null(init$z_ms_intercept)
})
