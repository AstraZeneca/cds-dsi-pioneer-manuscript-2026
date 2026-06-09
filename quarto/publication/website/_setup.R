# Setup for Publication Website
# This file is sourced by all publication analysis pages.
#
# Sources shared setup (common libs + theme) then loads publication-specific
# functions and data from the targets store.

here::i_am("quarto/publication/website/_setup.R")

# Load shared setup: renv, packages, utility sources, ggplot2 theme
source(here::here("quarto", "_shared", "_setup.R"))

# Source sclc plot, table, and accuracy functions for use in publication
source(here::here("r", "sclc", "plot_functions.R"))
source(here::here("r", "sclc", "table_functions.R"))
source(here::here("r", "sclc", "accuracy.R"))
source(here::here("r", "accuracy.R"))

# Define target stores (shared across website and other scripts)
source(here::here("r", "publication", "target_stores.R"))
