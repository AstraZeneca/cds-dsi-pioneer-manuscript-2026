library(testthat)
library(here)
library(stringr)

# Task 11 Step 3 (plan): long-horizon plateau assertion for the warped LFO forecast
# path. The growth arm must approach init_growth + growth_rate/kappa and never exceed it.

test_that("warped growth arm plateaus at init + growth_rate/kappa and never exceeds it", {
  kappa <- 0.05
  growth_rate <- 0.08
  init_growth <- log(2)
  # Forecast out to a very long horizon (weeks): anchor at 0, then far out.
  times <- c(0, seq(4, 4000, by = 4))
  data <- list(
    n_t = length(times),
    times = times,
    init_decrease = log(5),
    init_growth = init_growth,
    decrease_rate = 0.03,
    growth_rate = growth_rate,
    kappa = kappa
  )
  fit <- test_stan_function(
    here::here("tests/testthat/stan/test_gr_decay_plateau.stan"), data
  )
  draws <- posterior::as_draws_df(fit$draws())

  plateau <- get_stan_val(draws, "plateau")
  final_growth <- get_stan_val(draws, "final_growth")
  max_growth <- get_stan_val(draws, "max_growth")

  expected_plateau <- init_growth + growth_rate / kappa
  # 1e-6 absorbs CSV round-trip precision on the passed-in reals (init/growth_rate/kappa).
  expect_equal(plateau, expected_plateau, tolerance = 1e-6)

  # At horizon 4000 wk with kappa=0.05, exp(-kappa*t) ~ 0 => final ~= plateau.
  expect_equal(final_growth, expected_plateau, tolerance = 1e-3)

  # The growth arm is monotone increasing and bounded ABOVE by the plateau:
  # it must never exceed it (the whole point — no explosion).
  expect_lte(max_growth, expected_plateau + 1e-6)
})

test_that("smaller kappa => higher plateau (monotone in 1/kappa)", {
  base <- list(
    n_t = 3L, times = c(0, 100, 1000),
    init_decrease = log(5), init_growth = log(2),
    decrease_rate = 0.03, growth_rate = 0.08
  )
  p_big_k <- {
    fit <- test_stan_function(
      here::here("tests/testthat/stan/test_gr_decay_plateau.stan"),
      c(base, list(kappa = 0.1))
    )
    get_stan_val(posterior::as_draws_df(fit$draws()), "plateau")
  }
  p_small_k <- {
    fit <- test_stan_function(
      here::here("tests/testthat/stan/test_gr_decay_plateau.stan"),
      c(base, list(kappa = 0.02))
    )
    get_stan_val(posterior::as_draws_df(fit$draws()), "plateau")
  }
  # plateau = init + growth_rate/kappa, so smaller kappa => strictly higher plateau.
  expect_gt(p_small_k, p_big_k)
})
