# Phase 2: Simple 1D Laplace experiment
# ============================================================================
# Validates manual Laplace approximation by comparing three fitting modes:
#   mode 0 = Full HMC (z_i sampled as parameters)
#   mode 1 = Laplace marginalization (z_i integrated out)
#   mode 2 = Exact marginalization (known analytical form)
#
# For this Gaussian-Gaussian model, modes 1 and 2 should give identical
# results, and mode 0 should match within MCSE.
# ============================================================================

library(cmdstanr)
library(posterior)
library(dplyr)
library(tibble)
library(purrr)
library(tidyr)

# Use CmdStan develop (has laplace_marginal_tol support)
set_cmdstan_path("~/.cmdstan/cmdstan-develop")

# --- Simulate data ---
set.seed(42)
N <- 30
true_mu <- 2.0
true_sigma <- 1.5
obs_sd <- 1.0

z_true <- rnorm(N)
y <- rnorm(N, mean = true_mu + true_sigma * z_true, sd = obs_sd)

stan_data <- list(N = N, y = y, obs_sd = obs_sd)

# --- Compile model ---
model <- cmdstan_model(
  "stan/experiments/laplace_simple_test.stan",
  quiet = FALSE
)

# --- Fit all three modes ---
fit_modes <- list()
mode_names <- c("full_hmc", "laplace", "exact")

for (mode in 0:2) {
  cat(sprintf("\n=== Fitting mode %d (%s) ===\n", mode, mode_names[mode + 1]))

  fit <- model$sample(
    data = c(stan_data, list(laplace_mode = mode)),
    seed = 123,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 1000,
    iter_sampling = 2000,
    refresh = 500
  )

  fit_modes[[mode_names[mode + 1]]] <- fit
}

# --- Compare posteriors ---
cat("\n\n=== Posterior comparison ===\n\n")

summaries <- map(mode_names, \(name) {
  fit_modes[[name]]$summary(variables = c("mu", "sigma")) |>
    mutate(mode = name, .before = 1)
}) |>
  list_rbind()

print(summaries |> select(mode, variable, mean, sd, q5, q95, rhat, ess_bulk))

# --- Check agreement between Laplace and Exact ---
cat("\n\n=== Laplace vs Exact (should be near-identical) ===\n")

laplace_draws <- fit_modes$laplace$draws(variables = c("mu", "sigma"), format = "df")
exact_draws <- fit_modes$exact$draws(variables = c("mu", "sigma"), format = "df")
hmc_draws <- fit_modes$full_hmc$draws(variables = c("mu", "sigma"), format = "df")

compare_means <- tibble(
  parameter = c("mu", "sigma"),
  hmc_mean = c(mean(hmc_draws$mu), mean(hmc_draws$sigma)),
  laplace_mean = c(mean(laplace_draws$mu), mean(laplace_draws$sigma)),
  exact_mean = c(mean(exact_draws$mu), mean(exact_draws$sigma)),
  hmc_sd = c(sd(hmc_draws$mu), sd(hmc_draws$sigma)),
  laplace_sd = c(sd(laplace_draws$mu), sd(laplace_draws$sigma)),
  exact_sd = c(sd(exact_draws$mu), sd(exact_draws$sigma))
)

print(compare_means)

# Check that Laplace and Exact are within Monte Carlo SE
cat("\n=== Difference (Laplace - Exact) / MCSE ===\n")
n_eff <- min(
  fit_modes$laplace$summary("mu")$ess_bulk,
  fit_modes$exact$summary("mu")$ess_bulk
)
mcse_mu <- compare_means$exact_sd[1] / sqrt(n_eff)
mcse_sigma <- compare_means$exact_sd[2] / sqrt(n_eff)

cat(sprintf("mu:    diff = %.4f, MCSE = %.4f, ratio = %.2f\n",
            abs(compare_means$laplace_mean[1] - compare_means$exact_mean[1]),
            mcse_mu,
            abs(compare_means$laplace_mean[1] - compare_means$exact_mean[1]) / mcse_mu))
cat(sprintf("sigma: diff = %.4f, MCSE = %.4f, ratio = %.2f\n",
            abs(compare_means$laplace_mean[2] - compare_means$exact_mean[2]),
            mcse_sigma,
            abs(compare_means$laplace_mean[2] - compare_means$exact_mean[2]) / mcse_sigma))

cat("\nPhase 2 complete. If Laplace/Exact ratio < 2, the implementation is correct.\n")
