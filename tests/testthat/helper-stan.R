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
  fit
}

#' Create mock PFS data for testing estimate_kaplan_meier
create_mock_km_data <- function(n_patients = 20, max_t = 100, pfs_offset = 0) {
  # Ensure non-monotonic PFS values by using sample without replacement first, then sampling
  # Use a mix of times including some early ones to ensure non-monotonic ordering
  pfs_times <- c(
    sample(0:(max_t %/% 4), min(n_patients %/% 3, 5), replace = TRUE), # Some early times
    sample(
      (max_t %/% 4):(3 * max_t %/% 4),
      n_patients - min(n_patients %/% 3, 5),
      replace = TRUE
    ) # Later times
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

#' Extract a scalar value from a posterior draws data frame by named Stan variable.
#'
#' draws_array is always 3D [iterations, chains, variables]. Multi-dimensional
#' Stan arrays are stored as named scalars, e.g. x[1,2,3]. Use this helper
#' instead of numeric indexing which breaks when variable names sort
#' alphabetically (e.g. x[10,1] sorts before x[2,1]).
#'
#' @param draws_df A draws_df from posterior::as_draws_df(fit$draws())
#' @param var Stan variable name (string)
#' @param ... Integer indices (one per dimension); omit for scalar variables
#' @return Numeric scalar (first iteration value)
get_stan_val <- function(draws_df, var, ...) {
  indices <- c(...)
  vname   <- if (length(indices) == 0) var else sprintf("%s[%s]", var, paste(indices, collapse = ","))
  as.numeric(draws_df[[vname]][1])
}
