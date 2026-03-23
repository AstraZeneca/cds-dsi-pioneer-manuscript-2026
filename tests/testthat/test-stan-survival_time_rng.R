# Test for survival_time_rng in Stan
library(testthat)
library(cmdstanr)


# Define test cases for survival_time_rng (all required Stan data variables)
# Use only probabilities in [1e-5, 1-1e-5] to avoid Inf/-Inf in logit
cases <- list(
  # Typical case: log_cond_prob_surv = log(0.8), T = 10, not censored, obs_surv_time = 0, interval_censored = 0
  list(
    log_cond_prob_surv = rep(log(0.8), 10),
    T = 10L,
    obs_surv_time = 0L,
    right_censored = 0L,
    interval_censored = 0L
  ),
  # Censored case: log_cond_prob_surv = log(0.5), T = 5, censored, obs_surv_time = 0, interval_censored = 0
  list(
    log_cond_prob_surv = rep(log(0.5), 5),
    T = 5L,
    obs_surv_time = 0L,
    right_censored = 1L,
    interval_censored = 0L
  ),
  # Degenerate: log_cond_prob_surv = log(1-1e-5), T = 10, not censored, obs_surv_time = 0, interval_censored = 0
  list(
    log_cond_prob_surv = rep(log(1 - 1e-5), 10),
    T = 10L,
    obs_surv_time = 0L,
    right_censored = 0L,
    interval_censored = 0L
  ),
  # Edge: log_cond_prob_surv = log(1e-5), T = 10, not censored, obs_surv_time = 0, interval_censored = 0
  list(
    log_cond_prob_surv = rep(log(1e-5), 10),
    T = 10L,
    obs_surv_time = 0L,
    right_censored = 0L,
    interval_censored = 0L
  ),
  # min_time > 0: log_cond_prob_surv = log(0.7), T = 10, not censored, obs_surv_time = 3, interval_censored = 0
  list(
    log_cond_prob_surv = rep(log(0.7), 10),
    T = 10L,
    obs_surv_time = 3L,
    right_censored = 0L,
    interval_censored = 0L
  )
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
      pad <- rep(log(1 - 1e-5), MAX_T - x$T)
      c(x$log_cond_prob_surv, pad)
    })),
    dim = c(N_CASES, MAX_T)
  ),
  obs_surv_time = vapply(cases, function(x) x$obs_surv_time, integer(1)),
  right_censored = vapply(cases, function(x) x$right_censored, integer(1)),
  interval_censored = vapply(cases, function(x) x$interval_censored, integer(1))
)

test_that("survival_time_rng: empirical distribution properties are correct", {
  stan_file_path <- here::here(
    "tests",
    "testthat",
    "stan",
    "test_survival_time_rng.stan"
  )

  # Run with 1000 fixed_param draws to accumulate empirical RNG distribution
  mod <- cmdstan_model(stan_file_path, include_paths = here::here("stan"))
  fit <- mod$sample(
    data = stan_data,
    seed = 12345,
    chains = 1,
    iter_sampling = 1000,
    iter_warmup = 0,
    fixed_param = TRUE
  )

  draws <- fit$draws()

  # Extract [case, draw] matrix for a given output variable
  get_draw_matrix <- function(varname, n_cases, n_draws) {
    df <- posterior::as_draws_df(draws)
    mat <- matrix(NA_integer_, nrow = n_cases, ncol = n_draws)
    for (i in seq_len(n_cases)) {
      for (j in seq_len(n_draws)) {
        vname <- sprintf("%s[%d,%d]", varname, i, j)
        mat[i, j] <- as.numeric(df[[vname]])[1]
      }
    }
    mat
  }

  sampled_time     <- get_draw_matrix("sampled_time",     N_CASES, N_DRAWS)
  sampled_censored <- get_draw_matrix("sampled_censored", N_CASES, N_DRAWS)

  for (i in seq_len(N_CASES)) {
    times <- sampled_time[i, ]
    cens  <- sampled_censored[i, ]

    expect_true(
      all(times >= 0 & times <= stan_data$T[i]),
      label = sprintf("Case %d: sampled times within [0, T]", i)
    )
    expect_true(
      all(cens %in% c(0L, 1L)),
      label = sprintf("Case %d: censored flag is binary", i)
    )

    # Degenerate case: survival prob ≈ 1 → event always at time 0
    if (abs(stan_data$log_cond_prob_surv[i, 1] - log(1 - 1e-5)) < 1e-8) {
      expect_true(all(times == 0),   label = sprintf("Case %d: degenerate high-surv → time==0", i))
      expect_true(all(cens  == 0L),  label = sprintf("Case %d: degenerate high-surv → not censored", i))
    }

    # Edge case: survival prob ≈ 0 → event always at max time T
    if (abs(stan_data$log_cond_prob_surv[i, 1] - log(1e-5)) < 1e-8) {
      expect_true(all(times == stan_data$T[i]),
                  label = sprintf("Case %d: near-zero surv → time==T", i))
    }
  }
})
