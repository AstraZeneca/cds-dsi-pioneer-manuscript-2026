library(testthat)
library(here)

# Regression guard for BUG 1 (full-model forecast ignored kappa: it called the non-decay
# RNG with growth factor==1). This drives generate_patient_states_with_means_decay_rng —
# the function the fix routes through — and asserts the deterministic forecast actually
# DEPENDS on kappa. A regression to the non-decay RNG / factor==1 makes the forecast
# kappa-independent, so the two-kappa comparison below collapses and the test fails.

stan_file <- here::here("tests/testthat/stan/test_gr_decay_full_forecast.stan")

base_data <- list(
  baseline_week = 0,
  n_forecast = 60L,
  horizon_week = 800,      # long horizon so the plateau (and thus kappa) bites hard
  init_decrease = log(5),
  init_growth = log(2),
  decrease_rate = 0.03,
  growth_rate = 0.08,
  baseline_obs_value = 7,  # baseline SLD, linear units
  measure_sd = 0.1         # > 0 (RNG scale); only noise-free *_mean_* outputs are asserted
)

run_at_kappa <- function(kappa) {
  fit <- test_stan_function(stan_file, c(base_data, list(kappa = kappa)))
  draws <- posterior::as_draws_df(fit$draws())
  list(
    final_log_obs = get_stan_val(draws, "final_forecast_log_obs"),
    final_growth  = get_stan_val(draws, "final_growth_state"),
    plateau       = get_stan_val(draws, "growth_plateau")
  )
}

test_that("the full-model forecast RNG actually uses kappa (smaller kappa => higher long-horizon burden)", {
  small_k <- run_at_kappa(0.02)   # weak attenuation, high plateau
  large_k <- run_at_kappa(0.20)   # strong attenuation, low plateau

  # Sanity: plateau = init_growth + growth_rate/kappa, strictly higher for smaller kappa.
  expect_gt(small_k$plateau, large_k$plateau)

  # The KEY assertion BUG 1 would violate: the forecast must DIFFER between the two kappas.
  # Under the bug (factor==1, kappa ignored) these would be identical.
  expect_gt(small_k$final_growth - large_k$final_growth, 0.1)
  expect_gt(small_k$final_log_obs - large_k$final_log_obs, 0.1)
})

test_that("each forecast growth arm stays bounded by its kappa plateau (no explosion)", {
  for (kappa in c(0.02, 0.05, 0.20)) {
    res <- run_at_kappa(kappa)
    # Bounded above by the plateau (the whole point of the warp). Small tolerance for the
    # finite horizon not having fully reached the asymptote.
    expect_lte(res$final_growth, res$plateau + 1e-6)
  }
})
