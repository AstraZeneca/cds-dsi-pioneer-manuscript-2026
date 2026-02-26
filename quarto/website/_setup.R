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

# === OS Data ===

# OS KM data
jan26_os_km_data <- tar_read(all_tumor_ssls_km_os_rvar_ctdna_jan26)
jan26_cond_os_km_data <- tar_read(all_tumor_ssls_cond_km_os_rvar_ctdna_jan26)
jan26_obs_os_km_all <- tar_read(km_trial_os_ctdna_jan26)
jan26_obs_os_km_first_line <- tar_read(km_trial_os_naive_ctdna_jan26)

# Median OS data
jan26_quant_os_data <- tar_read(all_tumor_ssls_trial_quant_os_ctdna_jan26)
jan26_cond_quant_os_data <- tar_read(all_tumor_ssls_cond_quant_os_ctdna_jan26)

# OS-n data
jan26_os_n_data <- tar_read(all_tumor_ssls_forecast_os_n_rvar_ctdna_jan26)
jan26_cond_os_n_data <- tar_read(all_tumor_ssls_cond_forecast_os_n_rvar_ctdna_jan26)

# PDL1 observed OS KM data
jan26_obs_os_km_pdl1 <- tar_read(km_trial_os_pdl1_ctdna_jan26)
jan26_obs_os_km_pdl1_naive <- tar_read(km_trial_os_pdl1_naive_ctdna_jan26)
jan26_obs_os_km_part_e_pdl1 <- tar_read(km_trial_os_part_e_pdl1_ctdna_jan26)
