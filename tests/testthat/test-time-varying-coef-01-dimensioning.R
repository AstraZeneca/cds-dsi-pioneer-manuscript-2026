# Tests for the time_varying_coef_01 dimensioning rule.
#
# The 0->1 multistate transition has three time-varying-covariate operating
# modes that map onto two different coefficient-vector lengths:
#
#   continuous-time      (visit_gated_01 = 0): length = n_time_varying_covar
#   latent visit-gated   (visit_gated_01 = 1, latent_01 = 1): same — full set
#   observed visit-gated (visit_gated_01 = 1, latent_01 = 0): single coef
#
# Three layers must agree on the length: the Stan parameter declaration in
# stan/modules/multistate/parameters.stan, the prior generator in
# r/priors.R::get_multistate_priors, and the initializer in
# r/sclc/initializers_fixed.R / r/initializers_ms.R. Any mismatch
# triggers a "dims declared=(N); dims found=(M)" error from Stan at sample
# time. These tests pin the R-side contracts so that the cheap failure
# happens here rather than after a recompile + sample.

library(testthat)
library(tibble)        # for lst() used inside get_multistate_priors / initializers
library(rlang)         # for is_null() etc.
library(stringr)       # for str_c() used inside priors.R

source(here::here("r", "priors.R"))

# ---------------------------------------------------------------------------
# get_multistate_priors: prior length matches the rule
# ---------------------------------------------------------------------------

test_that("get_multistate_priors: continuous mode -> length n_time_varying_covar", {
  res <- get_multistate_priors(
    n_levels = 1L, n_time_varying_covar = 3L, n_time_invariant_covar = 2L,
    enable_ms_visit_gated_01 = 0L,
    enable_ms_visit_gated_latent_01 = 0L
  )
  expect_length(res$time_varying_coef_01_mean, 3L)
  expect_length(res$time_varying_coef_01_sd, 3L)
})

test_that("get_multistate_priors: latent visit-gated -> length n_time_varying_covar", {
  res <- get_multistate_priors(
    n_levels = 1L, n_time_varying_covar = 3L, n_time_invariant_covar = 2L,
    enable_ms_visit_gated_01 = 1L,
    enable_ms_visit_gated_latent_01 = 1L
  )
  expect_length(res$time_varying_coef_01_mean, 3L)
  expect_length(res$time_varying_coef_01_sd, 3L)
})

test_that("get_multistate_priors: observed visit-gated -> length 1", {
  res <- get_multistate_priors(
    n_levels = 1L, n_time_varying_covar = 3L, n_time_invariant_covar = 2L,
    enable_ms_visit_gated_01 = 1L,
    enable_ms_visit_gated_latent_01 = 0L
  )
  expect_length(res$time_varying_coef_01_mean, 1L)
  expect_length(res$time_varying_coef_01_sd, 1L)
})

test_that("get_multistate_priors: default args treat both flags as off (continuous)", {
  res <- get_multistate_priors(
    n_levels = 1L, n_time_varying_covar = 4L, n_time_invariant_covar = 2L
  )
  expect_length(res$time_varying_coef_01_mean, 4L)
})

# ---------------------------------------------------------------------------
# Initializers: time_varying_coef_01 init length matches the rule
# ---------------------------------------------------------------------------
#
# The initializers are wrapped in a stan_data env. Build a minimal stan_data
# with just the fields ms_init_values_fixed reads, then inspect the resulting
# init list. Anything not reached by the time_varying_coef_01 branch is left
# as a placeholder.

minimal_ms_stan_data <- function(visit_gated_01, visit_gated_latent_01,
                                  n_time_varying_covar = 3L) {
  list(
    n_time_varying_covar      = n_time_varying_covar,
    n_time_invariant_covar    = 0L,
    n_levels                  = 1L,
    n_groups_per_level        = c(patient = 10L),
    n_patients                = 10L,
    enable_ms_01              = 1L,
    enable_ms_02              = 0L,
    enable_ms_03              = 0L,
    enable_ms_12              = 0L,
    enable_ms_32              = 0L,
    ms_time_scale_12          = 1L,
    enable_ms_level_baseline_hazard = c(patient = 0L),
    enable_ms_level_cov       = c(patient = 0L),
    enable_ms_pop_time_varying_cov = 1L,
    enable_ms_pop_time_invariant_cov = 0L,
    enable_ms_visit_gated_01  = visit_gated_01,
    enable_ms_visit_gated_latent_01 = visit_gated_latent_01,
    enable_ms_02_time_varying_cov = 0L,
    enable_ms_03_time_invariant_cov = 0L,
    enable_ms_03_time_varying_cov = 0L,
    enable_ms_32_time_invariant_cov = 0L,
    enable_ms_12_entry_covar  = 0L,
    enable_ms_32_entry_covar  = 0L,
    max_all_t                 = 50L,
    ms_max_sojourn_t          = 50L,
    ms_max_sojourn_t_32       = 50L,
    ms_gp_grid_step           = 4L,
    enable_states_full_grid   = 0L
  )
}

test_that("ms_init_values_fixed: continuous mode init has length n_time_varying_covar", {
  source(here::here("r", "sclc", "initializers_fixed.R"))
  env <- minimal_ms_stan_data(visit_gated_01 = 0L, visit_gated_latent_01 = 0L,
                              n_time_varying_covar = 3L)
  init <- ms_init_values_fixed(env)
  expect_length(init$time_varying_coef_01, 3L)
})

test_that("ms_init_values_fixed: latent visit-gated init has length n_time_varying_covar", {
  source(here::here("r", "sclc", "initializers_fixed.R"))
  env <- minimal_ms_stan_data(visit_gated_01 = 1L, visit_gated_latent_01 = 1L,
                              n_time_varying_covar = 3L)
  init <- ms_init_values_fixed(env)
  expect_length(init$time_varying_coef_01, 3L)
})

test_that("ms_init_values_fixed: observed visit-gated init has length 1", {
  source(here::here("r", "sclc", "initializers_fixed.R"))
  env <- minimal_ms_stan_data(visit_gated_01 = 1L, visit_gated_latent_01 = 0L,
                              n_time_varying_covar = 3L)
  init <- ms_init_values_fixed(env)
  expect_length(init$time_varying_coef_01, 1L)
})

# ---------------------------------------------------------------------------
# State-3 covariate dimensioning rule
# ---------------------------------------------------------------------------
#
# 0->3 and 3->2 covariates are gated by per-transition flags
# (enable_ms_03_time_invariant_cov, enable_ms_03_time_varying_cov,
# enable_ms_32_time_invariant_cov). When the flag is OFF, both the Stan
# parameter declaration and the priors hyperparam declaration must be
# length 0 so R-side priors that send length 0 round-trip cleanly. When ON,
# they're full length (n_time_varying_covar / n_time_invariant_covar).
#
# Pre-merge bug B1: hyperparams.stan declared TI hyperparam vectors as
# unconditional `vector[n_time_invariant_covar]` while parameters.stan and
# priors.R both followed the gated rule, causing data-load failure.

test_that("get_multistate_priors: 0->3 TV coef length follows enable_ms_03_time_varying_cov", {
  off <- get_multistate_priors(
    n_levels = 1L, n_time_varying_covar = 3L, n_time_invariant_covar = 2L,
    enable_ms_03_time_varying_cov = 0L
  )
  expect_length(off$time_varying_coef_03_mean, 0L)
  expect_length(off$time_varying_coef_03_sd, 0L)

  on <- get_multistate_priors(
    n_levels = 1L, n_time_varying_covar = 3L, n_time_invariant_covar = 2L,
    enable_ms_03_time_varying_cov = 1L
  )
  expect_length(on$time_varying_coef_03_mean, 3L)
  expect_length(on$time_varying_coef_03_sd, 3L)
})

test_that("get_multistate_priors: 0->3 TI coef length follows enable_ms_03_time_invariant_cov", {
  off <- get_multistate_priors(
    n_levels = 1L, n_time_varying_covar = 3L, n_time_invariant_covar = 2L,
    enable_ms_03_time_invariant_cov = 0L
  )
  expect_length(off$time_invariant_coef_03_mean, 0L)
  expect_length(off$time_invariant_coef_03_sd, 0L)

  on <- get_multistate_priors(
    n_levels = 1L, n_time_varying_covar = 3L, n_time_invariant_covar = 2L,
    enable_ms_03_time_invariant_cov = 1L
  )
  expect_length(on$time_invariant_coef_03_mean, 2L)
  expect_length(on$time_invariant_coef_03_sd, 2L)
})

test_that("get_multistate_priors: 3->2 TI coef length follows enable_ms_32_time_invariant_cov", {
  off <- get_multistate_priors(
    n_levels = 1L, n_time_varying_covar = 3L, n_time_invariant_covar = 2L,
    enable_ms_32_time_invariant_cov = 0L
  )
  expect_length(off$time_invariant_coef_32_mean, 0L)
  expect_length(off$time_invariant_coef_32_sd, 0L)

  on <- get_multistate_priors(
    n_levels = 1L, n_time_varying_covar = 3L, n_time_invariant_covar = 2L,
    enable_ms_32_time_invariant_cov = 1L
  )
  expect_length(on$time_invariant_coef_32_mean, 2L)
  expect_length(on$time_invariant_coef_32_sd, 2L)
})

test_that("get_multistate_priors: defaults disable all state-3 covariates", {
  res <- get_multistate_priors(
    n_levels = 1L, n_time_varying_covar = 3L, n_time_invariant_covar = 2L
  )
  expect_length(res$time_varying_coef_03_mean, 0L)
  expect_length(res$time_invariant_coef_03_mean, 0L)
  expect_length(res$time_invariant_coef_32_mean, 0L)
})
