# Test for survival_time_rng in Stan
library(testthat)
library(cmdstanr)



# Define test cases for survival_time_rng (all required Stan data variables)
# Use only probabilities in [1e-5, 1-1e-5] to avoid Inf/-Inf in logit
cases <- list(
  # Typical case: log_cond_prob_surv = log(0.8), T = 10, not censored, obs_surv_time = 0, interval_censored = 0
  list(log_cond_prob_surv = rep(log(0.8), 10), T = 10L, obs_surv_time = 0L, right_censored = 0L, interval_censored = 0L),
  # Censored case: log_cond_prob_surv = log(0.5), T = 5, censored, obs_surv_time = 0, interval_censored = 0
  list(log_cond_prob_surv = rep(log(0.5), 5), T = 5L, obs_surv_time = 0L, right_censored = 1L, interval_censored = 0L),
  # Degenerate: log_cond_prob_surv = log(1-1e-5), T = 10, not censored, obs_surv_time = 0, interval_censored = 0
  list(log_cond_prob_surv = rep(log(1-1e-5), 10), T = 10L, obs_surv_time = 0L, right_censored = 0L, interval_censored = 0L),
  # Edge: log_cond_prob_surv = log(1e-5), T = 10, not censored, obs_surv_time = 0, interval_censored = 0
  list(log_cond_prob_surv = rep(log(1e-5), 10), T = 10L, obs_surv_time = 0L, right_censored = 0L, interval_censored = 0L),
  # min_time > 0: log_cond_prob_surv = log(0.7), T = 10, not censored, obs_surv_time = 3, interval_censored = 0
  list(log_cond_prob_surv = rep(log(0.7), 10), T = 10L, obs_surv_time = 3L, right_censored = 0L, interval_censored = 0L)
)

N_CASES <- length(cases)
MAX_T <- max(sapply(cases, function(x) x$T))

# Number of draws per case for Stan RNG loop
N_DRAWS <- 1000L
stan_data <- list(
  N_CASES = N_CASES,
  MAX_T = MAX_T,
  N_DRAWS = N_DRAWS,
  T = vapply(cases, function(x) x$T, integer(1)),
  log_cond_prob_surv = array(
    unlist(lapply(cases, function(x) {
      pad <- rep(log(1-1e-5), MAX_T - x$T)
      c(x$log_cond_prob_surv, pad)
    })),
    dim = c(N_CASES, MAX_T)
  ),
  obs_surv_time = vapply(cases, function(x) x$obs_surv_time, integer(1)),
  right_censored = vapply(cases, function(x) x$right_censored, integer(1)),
  interval_censored = vapply(cases, function(x) x$interval_censored, integer(1))
)

stan_file_path <- here::here("tests", "testthat", "stan", "test_survival_time_rng.stan")
cat("Resolved Stan file path: ", stan_file_path, "\n")
cat("File exists? ", file.exists(stan_file_path), "\n")

# Run Stan model with a fixed seed for reproducibility
mod <- cmdstan_model(stan_file_path, include_paths = here::here("stan"))
fit <- mod$sample(
  data = stan_data,
  seed = 12345,
  chains = 1,
  iter_sampling = 1000, # Number of RNG draws for empirical checks
  iter_warmup = 0,
  fixed_param = TRUE
)


draws <- fit$draws()

# Extract outputs as [case, draw] arrays
get_draw_matrix <- function(varname, n_cases, n_draws) {
  # Extracts [case, draw] matrix from draws_df
  df <- posterior::as_draws_df(draws)
  mat <- matrix(NA_integer_, nrow = n_cases, ncol = n_draws)
  for (i in seq_len(n_cases)) {
    for (j in seq_len(n_draws)) {
      vname <- sprintf("%s[%d,%d]", varname, i, j)
      # Each row in df is a draw, so use the column for this [i,j]
      mat[i, j] <- as.numeric(df[[vname]])[1]  # Only one draw per row in fixed_param
    }
  }
  mat
}

 n_draws <- N_DRAWS
sampled_time <- get_draw_matrix("sampled_time", N_CASES, n_draws)
sampled_censored <- get_draw_matrix("sampled_censored", N_CASES, n_draws)

# Test: empirical mean and quantiles for each case
for (i in seq_len(N_CASES)) {
  times <- sampled_time[i, ]
  cens <- sampled_censored[i, ]
  # Print unique values for debugging
  cat(sprintf("Case %d: unique censored values: %s\n", i, paste(unique(cens), collapse=", ")))
  # Check that all times are within [0, T[i]]
  expect_true(all(times >= 0 & times <= stan_data$T[i]))
  # Check that censored is always 0 or 1
  expect_true(all(cens %in% c(0, 1)))
  # For degenerate case (log(1-1e-5)), time should always be 0
  if (abs(stan_data$log_cond_prob_surv[i,1] - log(1-1e-5)) < 1e-8) {
    expect_true(all(times == 0))
    expect_true(all(cens == 0))
  }
  # For edge case (log(1e-5)), time should always be T[i]
  if (abs(stan_data$log_cond_prob_surv[i,1] - log(1e-5)) < 1e-8) {
    expect_true(all(times == stan_data$T[i]))
  }
}

# Optionally, print empirical means for manual inspection
cat("Empirical means for sampled_time per case:\n")
print(rowMeans(sampled_time))
cat("Empirical means for sampled_censored per case:\n")
print(rowMeans(sampled_censored))
