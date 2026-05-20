library(testthat)
library(posterior)
library(dplyr)
library(tidyr)
library(purrr)

source(here::here("r", "diagnostics.R"))

# ── Fixtures ──────────────────────────────────────────────────────────────────

make_draws_array <- function(n_chains = 2, n_iter = 100, params, chain_offset = 0) {
  set.seed(42)
  arr <- array(
    data = rnorm(n_chains * n_iter * length(params)),
    dim = c(n_iter, n_chains, length(params)),
    dimnames = list(
      iteration = NULL,
      chain = paste0("chain:", seq_len(n_chains)),
      variable = params
    )
  )
  if (chain_offset != 0) arr[, 1, ] <- arr[, 1, ] + chain_offset
  posterior::as_draws_array(arr)
}

make_nuts_params <- function(n_chains = 2, n_iter = 100,
                             n_divergent = rep(0, n_chains),
                             treedepths = NULL,
                             energies = NULL) {
  map(seq_len(n_chains), \(chain) {
    e  <- if (!is.null(energies)) energies[[chain]] else cumsum(rnorm(n_iter, sd = 5, mean = 0)) + 100
    td <- if (!is.null(treedepths)) treedepths[[chain]] else rep(5L, n_iter)
    divs <- c(rep(1, n_divergent[chain]), rep(0, n_iter - n_divergent[chain]))
    tibble(Chain = chain, Iteration = seq_len(n_iter)) |>
      expand_grid(Parameter = c("accept_stat__", "divergent__", "energy__",
                                "n_leapfrog__", "stepsize__", "treedepth__")) |>
      mutate(Value = case_when(
        Parameter == "accept_stat__" ~ 0.85,
        Parameter == "divergent__"   ~ divs[Iteration],
        Parameter == "energy__"      ~ e[Iteration],
        Parameter == "n_leapfrog__"  ~ 31,
        Parameter == "stepsize__"    ~ 0.05,
        Parameter == "treedepth__"   ~ as.double(td[Iteration])
      ))
  }) |>
    list_rbind()
}

# ── check_convergence ─────────────────────────────────────────────────────────

test_that("check_convergence returns list with pop and patient elements", {
  draws_pop     <- make_draws_array(params = c("tr_loc_pop", "frac_logit_loc_pop", "measure_sd_sld"))
  draws_patient <- make_draws_array(params = paste0("tr_intercept_patient[", 1:20, "]"))

  result <- check_convergence(draws_pop, draws_patient)

  expect_type(result, "list")
  expect_named(result, c("pop", "patient"))
})

test_that("check_convergence pop summary has one row per parameter, sorted by rhat descending", {
  pop_params    <- c("tr_loc_pop", "frac_logit_loc_pop", "measure_sd_sld")
  draws_pop     <- make_draws_array(params = pop_params)
  draws_patient <- make_draws_array(params = paste0("param[", 1:5, "]"))

  result <- check_convergence(draws_pop, draws_patient)

  expect_equal(nrow(result$pop), length(pop_params))
  expect_true("variable" %in% names(result$pop))
  expect_true("rhat"     %in% names(result$pop))
  expect_true("ess_bulk" %in% names(result$pop))
  expect_true("ess_tail" %in% names(result$pop))
  # sorted descending by rhat
  expect_true(all(diff(result$pop$rhat) <= 0))
})

test_that("check_convergence patient summary has correct structure and counts params", {
  n_patient_params <- 30
  draws_pop     <- make_draws_array(params = c("tr_loc_pop", "measure_sd_sld"))
  draws_patient <- make_draws_array(params = paste0("tr_intercept_patient[", seq_len(n_patient_params), "]"))

  result <- check_convergence(draws_pop, draws_patient)

  expect_equal(nrow(result$patient), 1L)
  expect_equal(result$patient$group, "patient_effects")
  expect_equal(result$patient$n_params, n_patient_params)
  expect_true(all(c("max_rhat", "n_rhat_bad", "min_ess_bulk", "min_ess_tail") %in% names(result$patient)))
})

test_that("check_convergence reports low rhat for well-mixed draws", {
  draws_pop     <- make_draws_array(params = c("tr_loc_pop", "measure_sd_sld"))
  draws_patient <- make_draws_array(params = paste0("param[", 1:10, "]"))

  result <- check_convergence(draws_pop, draws_patient)

  expect_true(max(result$pop$rhat, na.rm = TRUE) < 1.05)
  expect_true(result$patient$max_rhat < 1.05)
})

test_that("check_convergence detects poorly-mixed chains via elevated rhat", {
  draws_pop     <- make_draws_array(params = c("tr_loc_pop"), chain_offset = 20)
  draws_patient <- make_draws_array(params = paste0("param[", 1:5, "]"), chain_offset = 20)

  result <- check_convergence(draws_pop, draws_patient)

  expect_true(max(result$pop$rhat, na.rm = TRUE) > 1.1)
  expect_true(result$patient$n_rhat_bad >= 1L)
})

# ── summarize_nuts ────────────────────────────────────────────────────────────

test_that("summarize_nuts returns one row per chain with expected columns", {
  nuts_df <- make_nuts_params(n_chains = 3, n_iter = 100)

  result <- summarize_nuts(nuts_df)

  expect_equal(nrow(result), 3L)
  expect_true(all(c("Chain", "n_iter", "n_divergent", "n_max_treedepth",
                    "mean_accept", "ebfmi") %in% names(result)))
  expect_equal(result$Chain, 1:3)
})

test_that("summarize_nuts counts divergences correctly per chain", {
  nuts_df <- make_nuts_params(n_chains = 2, n_iter = 100, n_divergent = c(5, 12))

  result <- summarize_nuts(nuts_df)

  expect_equal(result$n_divergent[result$Chain == 1], 5L)
  expect_equal(result$n_divergent[result$Chain == 2], 12L)
})

test_that("summarize_nuts counts treedepth hits at default max_treedepth = 10", {
  # chain 1: all treedepth = 10 (5 hits at default); chain 2: all treedepth = 9 (0 hits)
  nuts_df <- make_nuts_params(
    n_chains = 2, n_iter = 50,
    treedepths = list(c(rep(10, 5), rep(4, 45)), rep(9, 50))
  )

  result <- summarize_nuts(nuts_df)

  expect_equal(result$n_max_treedepth[result$Chain == 1], 5L)
  expect_equal(result$n_max_treedepth[result$Chain == 2], 0L)
})

test_that("summarize_nuts respects custom max_treedepth", {
  # With max_treedepth = 8: treedepth 9 and 10 both count as hits
  nuts_df <- make_nuts_params(
    n_chains = 1, n_iter = 10,
    treedepths = list(c(rep(10, 3), rep(9, 4), rep(7, 3)))
  )

  result_default <- summarize_nuts(nuts_df)              # max = 10
  result_lower   <- summarize_nuts(nuts_df, max_treedepth = 8)  # max = 8

  expect_equal(result_default$n_max_treedepth, 3L)  # only treedepth == 10
  expect_equal(result_lower$n_max_treedepth,   7L)  # treedepth >= 9
})

test_that("summarize_nuts computes mean_accept correctly", {
  nuts_df <- make_nuts_params(n_chains = 1, n_iter = 4)
  # All accept_stat__ = 0.85 from fixture
  result <- summarize_nuts(nuts_df)
  expect_equal(result$mean_accept, 0.85)
})

test_that("summarize_nuts computes ebfmi as var(diff(energy)) / var(energy) per chain", {
  set.seed(7)
  energy <- cumsum(rnorm(100, sd = 3)) + 50
  expected_ebfmi <- stats::var(diff(energy)) / stats::var(energy)

  nuts_df <- make_nuts_params(n_chains = 1, n_iter = 100, energies = list(energy))
  result  <- summarize_nuts(nuts_df)

  expect_equal(result$ebfmi, expected_ebfmi, tolerance = 1e-10)
})

test_that("summarize_nuts n_iter matches number of iterations per chain", {
  nuts_df <- make_nuts_params(n_chains = 2, n_iter = 73)
  result  <- summarize_nuts(nuts_df)
  expect_true(all(result$n_iter == 73L))
})
