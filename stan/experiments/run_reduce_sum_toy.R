library(cmdstanr)

set.seed(42)
N <- 100
y <- rnorm(N, mean = 0.2, sd = 0.8)

mod <- cmdstan_model(exe_file = "stan/experiments/reduce_sum_toy")
mod$.__enclos_env__$private$cpp_options_$stan_threads <- TRUE

fit <- mod$sample(
  data = list(N = N, y = y, grainsize = 1),
  chains = 4,
  parallel_chains = 4,
  threads_per_chain = 2,
  iter_warmup = 500,
  iter_sampling = 500,
  seed = 42,
  refresh = 250
)

fit$summary(c("mu", "sigma"))
