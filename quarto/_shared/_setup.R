# Shared setup for PIONEER Quarto documents
# This file should be sourced by trial-specific _setup.R files.
#
# Usage: source(here::here("quarto", "_shared", "_setup.R"))
#
# Provides: renv init, common packages, shared utility sources, ggplot2 theme.
# Does NOT load any trial-specific data or set targets stores.

# here::here() auto-detects project root from .git directory
repo_root <- here::here()

# Only source .Rprofile if init_project isn't already defined
# (which means .Rprofile was already sourced by R startup)
if (!exists("init_project", mode = "function")) {
  # Activate renv first (using absolute path)
  renv_activate <- file.path(repo_root, "renv", "activate.R")
  if (file.exists(renv_activate)) {
    source(renv_activate)
  }

  # Then source .Rprofile with chdir to handle relative paths
  rprofile_path <- file.path(repo_root, ".Rprofile")
  if (file.exists(rprofile_path)) {
    source(rprofile_path, chdir = TRUE)
  }
}

# Initialize the project (loads renv, sets paths, etc.)
init_project()

# Load required libraries
library(loo)
library(ggdist)
library(bayesplot)
library(colorspace)
library(ggrepel)
library(tidycmprsk)
library(gt)
library(patchwork)
library(ggalluvial)

# Source shared utility functions
source(here::here("r", "util.R"))
source(here::here("r", "plot_functions.R"))
source(here::here("r", "table_functions.R"))
source(here::here("r", "targets_tidyselect.R"))

# Set default ggplot theme with proper margins to prevent caption cutoff
ggplot2::theme_set(
  theme_minimal(base_family = "") +
  ggplot2::theme(
    plot.margin = ggplot2::margin(5, 10, 20, 20, "pt"),
    plot.caption = ggplot2::element_text(hjust = 0, margin = ggplot2::margin(15, 0, 0, 0, "pt"))
  )
)
