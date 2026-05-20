if (file.exists("~/.Rprofile")) {
  source("~/.Rprofile")
}

# Find and activate renv - try relative path first (when R starts from project root),
# then fall back to here::here() for other cases (like Quarto subdirectories)
if (file.exists("renv/activate.R")) {
  source("renv/activate.R")
} else if (requireNamespace("here", quietly = TRUE)) {
  source(here::here("renv", "activate.R"))
}

# if (is_domino) {
data_path <- Sys.getenv("DOMINO_DATASETS_DIR")
output_path <- file.path(data_path, "analysis-results", Sys.getenv("DOMINO_STARTING_USERNAME"))
artifacts_path <- file.path(Sys.getenv("DOMINO_ARTIFACTS_DIR"), Sys.getenv("DOMINO_STARTING_USERNAME"))
fit_output_timestamp <- FALSE

library(conflicted)

conflicts_prefer(
  dplyr::filter,
  dplyr::lag,
  posterior::sd,
  posterior::mad,
  rlang::set_names,
  purrr::flatten_dbl,
)

options(yaml.eval.expr = TRUE)

# Set default TAR_RUN if not provided by the shell environment
if (Sys.getenv("TAR_RUN") == "") {
  Sys.setenv(TAR_RUN = "main")
}

# Set AZ colour scheme
AZ_plum <- "#830051"
AZ_gold <- "#F0AB00"
AZ_turquoise <- "#68D2DF"
AZ_darkpurple <- "#3C1053"
AZ_green <- "#C4D600"
AZ_navy <- "#003865"
AZ_platinum <- "#9DB0AC"
AZ_darkgrey <- "#3F4444"
AZ_pink <- "#D0006F"
AZ_lightpurple <- "#AE53DE"
AZ_palette <- c(
  AZ_plum,
  AZ_gold,
  AZ_navy,
  AZ_turquoise,
  AZ_darkpurple,
  AZ_green,
  AZ_platinum,
  AZ_darkgrey,
  AZ_pink,
  AZ_lightpurple
)

# RECIST category colors (standardized across all plots)
RECIST_COLORS <- c(
  "CR" = AZ_green,
  "PR" = AZ_turquoise,
  "SD" = AZ_gold,
  "PD" = AZ_plum
)

# Helper function to create RECIST color scales
scale_fill_recist <- function(...) {
  ggplot2::scale_fill_manual(values = RECIST_COLORS, name = "RECIST Category", ...)
}

scale_color_recist <- function(...) {
  ggplot2::scale_color_manual(values = RECIST_COLORS, name = "RECIST Category", ...)
}

init_project <- function(output_path = output_path, artifacts_path = artifacts_path) {
  library(magrittr)
  library(tidyverse)
  library(rlang)
  library(targets)
  library(tarchetypes)
  library(stantargets)
  library(crew)
  library(autometric)
  library(here)
  library(cmdstanr)
  set_cmdstan_path("~/.cmdstan/cmdstan-2.38.0")
  library(posterior)
  library(tidybayes)
  library(qs2)
  library(recipes)

  source(here("r", "util.R"))
  source(here("r", "multi_level_hierarchy.R"))  # Required by priors.R
  source(here("r", "priors.R"))
  source(here("r", "posterior.R"))
  source(here("r", "prepare_analysis_data.R"))
  source(here("r", "initializers_ms.R"))
  source(here("r", "initializers.R"))
  source(here("r", "accuracy.R"))
  source(here("r", "state_space.R"))
  source(here("r", "plot_functions.R"))
  source(here("r", "parquet_draws.R"))
  source(here("r", "diagnostics.R"))

  tar_project <- Sys.getenv("TAR_PROJECT", "sclc")

  if (tar_project == "sclc") {
    source(here("r", "sclc", "priors.R"))
    source(here("r", "multistate.R"))
    source(here("r", "sclc", "prepare_analysis_data.R"))
    source(here("r", "sclc", "accuracy.R"))
    source(here("r", "sclc", "initializers.R"))
    source(here("r", "sclc", "plot_functions.R"))
  } else if (tar_project == "pioneer") {
    source(here("r", "pioneer", "prepare_analysis_data.R"))
    source(here("r", "pioneer", "prepare_laplace_data.R"))
    source(here("r", "pioneer", "priors.R"))
    source(here("r", "pioneer", "initializers.R"))
    source(here("r", "pioneer", "plot_functions.R"))
    source(here("r", "pioneer", "prepare_ms_standalone_lfo_data.R"))
  }

  source(here("r", "targets_tidyselect.R"))
}
