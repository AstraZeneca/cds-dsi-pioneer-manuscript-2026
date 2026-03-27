library(targets)
library(tidybayes)
library(dplyr)
library(stringr)
library(purrr)
library(posterior)
library(rlang)
library(qs2)
source("r/util.R")

store <- file.path(
  Sys.getenv("DOMINO_DATASETS_DIR"),
  "analysis-results",
  Sys.getenv("DOMINO_STARTING_USERNAME"),
  "pioneer/propensity-full/_targets"
)

cat("Loading fit...\n")
fit <- tar_read(pioneer_fit_res_posterior_propensity, store = store)

cat("Selecting likelihood_weight draws...\n")
w_draws <- select_draws(fit, starts_with("likelihood_weight"))

cat("Spreading rvars...\n")
w <- w_draws |>
  spread_rvars(likelihood_weight[patient]) |>
  median_qi(likelihood_weight, .width = c(0.5, 0.9))

out_path <- "/mnt/artifacts/propensity_weights_median_qi.qs"
cat("Saving to", out_path, "\n")
qs_save(w, out_path)

cat("Done. Rows:", nrow(w), "\n")
print(quantile(w$likelihood_weight, c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1)))
