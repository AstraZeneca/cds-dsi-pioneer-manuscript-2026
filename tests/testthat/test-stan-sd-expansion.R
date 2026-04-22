library(testthat)
library(here)
library(posterior)
library(purrr)
library(dplyr)

# R oracle: expand SD array from mode vector, FE hyperparams, and RE free values.
# Mirrors the Stan transformed_parameters SD expansion pattern.
r_expand_sd <- function(mode, fe_sd, re_sd_values) {
  n_levels <- length(mode)
  sd_out <- numeric(n_levels)
  sd_idx <- 0L
  for (lv in seq_len(n_levels)) {
    if (mode[lv] == 1L) {
      sd_out[lv] <- fe_sd[lv]
    } else if (mode[lv] == 2L) {
      sd_idx <- sd_idx + 1L
      sd_out[lv] <- re_sd_values[sd_idx]
    } else {
      sd_out[lv] <- 0.0
    }
  }
  sd_out
}

run_sd_expansion_test <- function(mode, fe_sd, re_sd_values) {
  n_levels <- length(mode)
  n_re_levels <- sum(mode == 2L)
  stopifnot(length(re_sd_values) == n_re_levels)

  stan_data <- list(
    n_levels    = n_levels,
    mode        = as.array(as.integer(mode)),
    fe_sd       = as.array(fe_sd),
    n_re_levels = n_re_levels,
    re_sd_values = if (n_re_levels > 0) as.array(re_sd_values) else array(numeric(0), dim = 0)
  )

  fit <- test_stan_function(
    here("tests", "testthat", "stan", "test_sd_expansion.stan"),
    data = stan_data
  )
  d <- as_draws_df(fit$draws())

  expected <- r_expand_sd(mode, fe_sd, re_sd_values)

  for (lv in seq_len(n_levels)) {
    expect_equal(
      as.numeric(d[[paste0("sd_expanded[", lv, "]")]]),
      expected[lv],
      tolerance = 1e-6,
      label = sprintf("sd_expanded[%d] (mode=%d)", lv, mode[lv])
    )
  }
}

test_that("SD expansion: all none (mode=0)", {
  run_sd_expansion_test(
    mode         = c(0L, 0L, 0L),
    fe_sd        = c(0.1, 0.2, 0.3),
    re_sd_values = numeric(0)
  )
})

test_that("SD expansion: all FE (mode=1)", {
  run_sd_expansion_test(
    mode         = c(1L, 1L, 1L),
    fe_sd        = c(0.01, 0.02, 0.03),
    re_sd_values = numeric(0)
  )
})

test_that("SD expansion: all RE (mode=2)", {
  run_sd_expansion_test(
    mode         = c(2L, 2L, 2L),
    fe_sd        = c(0.1, 0.2, 0.3),
    re_sd_values = c(0.4, 0.5, 0.6)
  )
})

test_that("SD expansion: mixed none/FE/RE in 3-level hierarchy", {
  run_sd_expansion_test(
    mode         = c(0L, 1L, 2L),
    fe_sd        = c(0.1, 0.05, 0.2),
    re_sd_values = c(0.35)
  )
})

test_that("SD expansion: mixed FE/none/RE in 3-level hierarchy", {
  run_sd_expansion_test(
    mode         = c(1L, 0L, 2L),
    fe_sd        = c(0.01, 0.99, 0.3),
    re_sd_values = c(0.7)
  )
})

test_that("SD expansion: RE/FE/none/RE/FE in 5-level hierarchy", {
  run_sd_expansion_test(
    mode         = c(2L, 1L, 0L, 2L, 1L),
    fe_sd        = c(0.1, 0.02, 0.3, 0.4, 0.005),
    re_sd_values = c(0.6, 0.8)
  )
})

test_that("SD expansion: single FE level", {
  run_sd_expansion_test(
    mode         = c(1L),
    fe_sd        = c(1e-4),
    re_sd_values = numeric(0)
  )
})

test_that("SD expansion: single RE level", {
  run_sd_expansion_test(
    mode         = c(2L),
    fe_sd        = c(0.1),
    re_sd_values = c(0.5)
  )
})

test_that("SD expansion: typical 3-level pioneer config (none/re/re)", {
  # trial=none, arm=re, patient=re
  run_sd_expansion_test(
    mode         = c(0L, 2L, 2L),
    fe_sd        = c(0.1, 0.35, 0.35),
    re_sd_values = c(0.35, 0.40)
  )
})

test_that("SD expansion: FE arm level (fe use case)", {
  # trial=none, arm=fe, patient=re
  run_sd_expansion_test(
    mode         = c(0L, 1L, 2L),
    fe_sd        = c(0.1, 1e-4, 0.35),
    re_sd_values = c(0.40)
  )
})
