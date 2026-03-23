# Shared setup for SCLC-01 Analysis Website
# This file is sourced by all analysis pages

# Get the repository root (2 levels up from quarto/website/)
# Working directory is quarto/website/ due to execute-dir: project in _quarto.yml
repo_root <- normalizePath(file.path(getwd(), "../.."))

# Set up here package to use the repo root
here::i_am("quarto/website/_setup.R")

# Check if init_project exists, if not source .Rprofile
if (!exists("init_project", mode = "function")) {
  rprofile_path <- file.path(repo_root, ".Rprofile")
  if (file.exists(rprofile_path)) {
    source(rprofile_path)
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

# Source utility functions using here() which now works after init_project()
source(here("r", "util.R"))
source(here("r", "plot_functions.R"))
source(here("r", "table_functions.R"))
source(here("r", "targets_tidyselect.R"))
source(here("r", "sclc", "plot_functions.R"))
source(here("r", "sclc", "table_functions.R"))
source(here("r", "sclc", "accuracy.R"))

# Define target stores (shared across website and other scripts)
source(here("r", "sclc", "target_stores.R"))

# Set default targets store (note: per CLAUDE.md, avoid tar_config_set for global state)
# Instead, use explicit store arguments in tar_read/tar_load calls
# This line is kept for backward compatibility but should be migrated to explicit stores
tar_config_set(store = analysis_store)

# Set default ggplot theme with proper margins to prevent caption cutoff
ggplot2::theme_set(
  theme_minimal(base_family = "Arial") +
  ggplot2::theme(
    plot.margin = ggplot2::margin(5, 10, 20, 20, "pt"),
    plot.caption = ggplot2::element_text(hjust = 0, margin = ggplot2::margin(15, 0, 0, 0, "pt"))
  )
)

