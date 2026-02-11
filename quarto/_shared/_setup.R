# Shared setup for SCLC-01 Quarto documents
# This file can be sourced by website pages and presentations
#
# Usage: source(here::here("quarto", "_shared", "_setup.R"))
# Or:    source("../../_shared/_setup.R") from quarto/website/
# Or:    source("../../_shared/_setup.R") from quarto/presentations/*/

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

# Source utility functions using here() which now works after init_project()
source(here::here("r", "util.R"))
source(here::here("r", "plot_functions.R"))
source(here::here("r", "table_functions.R"))
source(here::here("r", "targets_tidyselect.R"))
source(here::here("r", "sclc", "plot_functions.R"))
source(here::here("r", "sclc", "table_functions.R"))
source(here::here("r", "sclc", "accuracy.R"))

# Define target stores (shared across website and other scripts)
source(here::here("r", "sclc", "target_stores.R"))

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

# Load common data used across multiple pages
all_analysis_data <- tar_read(all_analysis_data_ctdna_jan26)
all_stan_data <- tar_read(all_stan_data_ctdna_jan26)

# === Jan26 DCO Data (current DCO - used across pages) ===

# KM data
jan26_km_data <- tar_read(all_tumor_ssls_km_rvar_ctdna_jan26)
jan26_cond_km_data <- tar_read(all_tumor_ssls_cond_km_rvar_ctdna_jan26)
jan26_obs_km_all <- tar_read(km_trial_pfs_ctdna_jan26)
jan26_obs_km_first_line <- tar_read(km_trial_pfs_naive_ctdna_jan26)

# ORR data
jan26_orr_data <- tar_read(all_tumor_ssls_orr_rvar_ctdna_jan26)
jan26_cond_orr_data <- tar_read(all_tumor_ssls_cond_orr_rvar_ctdna_jan26)

# Median PFS data
jan26_quant_pfs_data <- tar_read(all_tumor_ssls_trial_quant_pfs_ctdna_jan26)
jan26_cond_quant_pfs_data <- tar_read(all_tumor_ssls_cond_quant_pfs_ctdna_jan26)

# PFS-n data
jan26_pfs_n_data <- tar_read(all_tumor_ssls_forecast_target_pfs_n_rvar_ctdna_jan26)
jan26_cond_pfs_n_data <- tar_read(all_tumor_ssls_cond_forecast_target_pfs_n_rvar_ctdna_jan26)
pfs_timepoints <- tar_read(pfs_timepoints)

# PDL1 observed KM data
jan26_obs_km_pdl1 <- tar_read(km_trial_pfs_pdl1_ctdna_jan26)
jan26_obs_km_pdl1_naive <- tar_read(km_trial_pfs_pdl1_naive_ctdna_jan26)
jan26_obs_km_part_e_pdl1 <- tar_read(km_trial_pfs_part_e_pdl1_ctdna_jan26)
