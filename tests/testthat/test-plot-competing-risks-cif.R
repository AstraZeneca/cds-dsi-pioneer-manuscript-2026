library(testthat)
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(forcats)
library(ggplot2)
library(ggdist)
library(patchwork)
library(posterior)
library(tidybayes)
library(scales)
library(here)

# Minimal stubs that plot_functions.R references at call time — avoids loading
# util.R (which drags in targets). Assigned into globalenv so function lookups
# from within plot_functions.R closures resolve correctly.
local({
  weeks_to_months <- function(weeks) weeks * 7 * 12 / 365.25
  assign("weeks_to_months",       weeks_to_months,                             envir = globalenv())
  assign("months_to_weeks",       function(months) months / weeks_to_months(1), envir = globalenv())
  assign("label_weeks_to_months", scales::label_number(scale = weeks_to_months(1)), envir = globalenv())
  assign("AZ_navy",      "#003865",  envir = globalenv())
  assign("AZ_plum",      "#830051",  envir = globalenv())
  assign("AZ_gold",      "#F0AB00",  envir = globalenv())
  assign("AZ_turquoise", "#68D2DF",  envir = globalenv())
  assign("AZ_pink",      "#D0006F",  envir = globalenv())
  assign("AZ_platinum",  "#9DB0AC",  envir = globalenv())
  assign("AZ_palette",   c("#003865", "#F0AB00", "#68D2DF", "#830051", "#D0006F", "#9DB0AC"),
         envir = globalenv())
})

source(here("r", "plot_functions.R"))

# =============================================================================
# Helpers
# =============================================================================

# Minimal stan_data for plot_competing_risks_cif.
# Uses a mix of event types so cmprsk::cuminc() gets a valid competing-risks
# dataset (at least one event and one censoring per trial).
# max_all_t must be >= time_step so at least one time point exists in the grid.
make_cif_stan_data <- function(n_patients = 6, n_trials = 1) {
  per_trial <- n_patients %/% n_trials
  trial_factor <- factor(
    rep(seq_len(n_trials), each = per_trial),
    levels = seq_len(n_trials)
  )
  n <- n_trials * per_trial

  # Event pattern (recycled across patients):
  # progression (0→1), on-trial death (0→2), admin-censored (state 0),
  # dropout (state 3), progression+censor, on-trial death
  final_states  <- rep_len(c(1L, 2L, 0L, 3L, 1L, 2L), n)
  ms_time_01    <- ifelse(final_states == 1L, 4L, 0L)
  ms_censored_01 <- ifelse(final_states == 1L, 0L, 1L)
  ms_time_02    <- ifelse(final_states == 2L, 5L, 0L)
  ms_censored_02 <- ifelse(final_states == 2L, 0L, 1L)
  ms_time_03    <- ifelse(final_states %in% c(0L, 3L), 8L, 0L)

  list(
    patient_trial   = trial_factor,
    max_all_t       = 8L,
    ms_censored_01  = as.integer(ms_censored_01),
    ms_time_01      = as.integer(ms_time_01),
    ms_censored_02  = as.integer(ms_censored_02),
    ms_time_02      = as.integer(ms_time_02),
    ms_final_state  = as.integer(final_states),
    ms_time_03      = as.integer(ms_time_03)
  )
}

# Flat draws matrix: columns named `spop_cif_0{cause}[t_idx]`.
# time_grid = seq(time_step, max_all_t, by = time_step) → t_idx = t + 1.
# With time_step=4 and max_all_t=8: grid = c(4, 8) → t_idx = c(5, 9).
make_flat_draws <- function(n_draws = 100, cif_prefix = "spop",
                            time_grid = c(4L, 8L), causes = 1:3,
                            cif_value = 0.1) {
  t_idxs <- time_grid + 1L
  cols <- unlist(lapply(causes, function(cause) {
    paste0(cif_prefix, "_cif_0", cause, "[", t_idxs, "]")
  }))
  mat <- matrix(cif_value, nrow = n_draws, ncol = length(cols),
                dimnames = list(NULL, cols))
  as_draws_matrix(mat)
}

# 2D draws matrix: columns named `spop_cif_0{cause}[trial,t_idx]`.
make_2d_draws <- function(n_draws = 100, cif_prefix = "spop",
                          time_grid = c(4L, 8L), n_trials = 2, causes = 1:3,
                          cif_value = 0.1) {
  t_idxs <- time_grid + 1L
  cols <- unlist(lapply(causes, function(cause) {
    unlist(lapply(seq_len(n_trials), function(tr) {
      paste0(cif_prefix, "_cif_0", cause, "[", tr, ",", t_idxs, "]")
    }))
  }))
  mat <- matrix(cif_value, nrow = n_draws, ncol = length(cols),
                dimnames = list(NULL, cols))
  as_draws_matrix(mat)
}

# =============================================================================
# flat_cif = TRUE: column naming and shape
# =============================================================================

test_that("flat_cif=TRUE extracts from 1D column names (no trial index)", {
  stan_data <- make_cif_stan_data(n_patients = 4, n_trials = 1)
  draws     <- make_flat_draws(cif_value = 0.1)

  p <- plot_competing_risks_cif(stan_data, draws, flat_cif = TRUE)

  expect_s3_class(p, "gg")
})

test_that("flat_cif=TRUE: model_cif has one row per time point (not per trial × time)", {
  stan_data <- make_cif_stan_data(n_patients = 4, n_trials = 1)
  # time_step=4, max_all_t=8 → two time points
  draws <- make_flat_draws(cif_value = 0.2, time_grid = c(4L, 8L))

  # Intercept model_cif by wrapping the internal call — we verify via the plot
  # data layer that exactly 2 time points × 3 causes appear in model_cif_long.
  p <- plot_competing_risks_cif(stan_data, draws, flat_cif = TRUE, time_step = 4L)
  plot_data <- ggplot2::ggplot_build(p$patches$plots[[1]])$data[[1]]

  # 2 time points × 3 causes = 6 rows in the CIF ribbon layer
  expect_equal(nrow(plot_data), 6L)
})

test_that("flat_cif=TRUE: median CIF value matches draws input", {
  stan_data <- make_cif_stan_data(n_patients = 4, n_trials = 1)
  draws     <- make_flat_draws(cif_value = 0.3)

  p <- plot_competing_risks_cif(stan_data, draws, flat_cif = TRUE, time_step = 4L)

  # Extract ribbon layer from first panel — ymin and ymax should both be ~0.3
  # (constant draws → lo = med = hi = 0.3)
  ribbon_data <- ggplot2::ggplot_build(p$patches$plots[[1]])$data[[1]]
  expect_true(all(abs(ribbon_data$ymin - 0.3) < 1e-6))
  expect_true(all(abs(ribbon_data$ymax - 0.3) < 1e-6))
})

test_that("flat_cif=TRUE: missing column returns NA without error", {
  stan_data <- make_cif_stan_data(n_patients = 4, n_trials = 1)
  # Only provide cause 1; causes 2 and 3 missing → should produce NA, not error
  draws <- make_flat_draws(cif_value = 0.1, causes = 1L)

  expect_no_error(
    plot_competing_risks_cif(stan_data, draws, flat_cif = TRUE, time_step = 4L)
  )
})

# =============================================================================
# flat_cif = FALSE: 2D column naming still works (sclc path unchanged)
# =============================================================================

test_that("flat_cif=FALSE (default) extracts from [trial,t] column names", {
  stan_data <- make_cif_stan_data(n_patients = 4, n_trials = 2)
  draws     <- make_2d_draws(n_trials = 2, cif_value = 0.15)

  p <- plot_competing_risks_cif(stan_data, draws, flat_cif = FALSE, time_step = 4L)

  expect_s3_class(p, "gg")
})

test_that("flat_cif=FALSE: model_cif has n_trials × n_time_points rows", {
  stan_data <- make_cif_stan_data(n_patients = 4, n_trials = 2)
  draws     <- make_2d_draws(n_trials = 2, cif_value = 0.15, time_grid = c(4L, 8L))

  p <- plot_competing_risks_cif(stan_data, draws, flat_cif = FALSE, time_step = 4L)
  # 2 trials × 2 times × 3 causes = 12 rows in ribbon layer
  ribbon_data <- ggplot2::ggplot_build(p$patches$plots[[1]])$data[[1]]
  expect_equal(nrow(ribbon_data), 12L)
})

# =============================================================================
# flat_cif=TRUE rejects 2D column names (wrong path guard)
# =============================================================================

test_that("flat_cif=TRUE with 2D draws returns NA CIF (no crash, wrong columns)", {
  stan_data <- make_cif_stan_data(n_patients = 4, n_trials = 1)
  # 2D column names won't match flat lookup → extract_cif returns NA rows
  draws <- make_2d_draws(n_trials = 1, cif_value = 0.5)

  # Should not error; CIF values will be NA since columns won't be found
  expect_no_error(
    plot_competing_risks_cif(stan_data, draws, flat_cif = TRUE, time_step = 4L)
  )
})
