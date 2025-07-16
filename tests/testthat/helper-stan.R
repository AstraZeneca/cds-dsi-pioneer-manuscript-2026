# tests/testthat/helper-stan.R
# 
# Helper functions for testing Stan functions, particularly estimate_kaplan_meier
#
# IMPORTANT NOTE about Stan's estimate_kaplan_meier output convention:
# - The survival vector starts with S(0), S(1), S(2), ..., S(max_t)
# - S(0) implicitly accounts for immediate events (pfs=0)
# - If no immediate events: S(0) = 1.0
# - If immediate events occur: S(0) = (n_patients - immediate_events) / n_patients
# - When comparing with R's survfit, account for this indexing difference
#
library(cmdstanr)
library(testthat)

#' Compile and run a Stan function test
#' @param stan_file Path to the test Stan file
#' @param data List of data to pass to Stan
#' @param ... Additional arguments to pass to cmdstan_model$sample()
test_stan_function <- function(stan_file, data, ...) {
  # Convert to absolute path if needed
  if (!file.exists(stan_file)) {
    stan_file <- here::here(stan_file)
  }
  # Always force recompilation to ensure includes are up to date
  mod <- cmdstan_model(
    stan_file,
    include_paths = here::here("stan"),
    quiet = FALSE,
    force_recompile = TRUE
  )
  fit <- mod$sample(
    data = data,
    chains = 1,
    iter_sampling = 1,
    iter_warmup = 0,
    fixed_param = TRUE,
   show_messages = TRUE,
    refresh = 0,
    ...
  )
  return(fit)
}

#' Create mock PFS data for testing estimate_kaplan_meier
create_mock_km_data <- function(n_patients = 20, max_t = 100, pfs_offset = 0) {
  # Ensure non-monotonic PFS values by using sample without replacement first, then sampling
  # Use a mix of times including some early ones to ensure non-monotonic ordering
  pfs_times <- c(
    sample(0:(max_t %/% 4), min(n_patients %/% 3, 5), replace = TRUE),  # Some early times
    sample((max_t %/% 4):(3 * max_t %/% 4), n_patients - min(n_patients %/% 3, 5), replace = TRUE)  # Later times
  )
  # Shuffle to ensure non-monotonic ordering
  pfs_times <- sample(pfs_times)
  
  list(
    n_patients = n_patients,
    event_time = pfs_times,
    right_censored = rbinom(n_patients, 1, 0.3),
    max_t = max_t,
    pfs_offset = pfs_offset
  )
}
