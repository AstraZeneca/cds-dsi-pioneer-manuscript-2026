if (file.exists("~/.Rprofile")) {
  source("~/.Rprofile")
}

options(
  renv.config.synchronized.check = FALSE
)

source("renv/activate.R")

# Lazy install: installs from lockfile if missing, then loads
lazy_lib <- function(pkg) {
  pkg <- deparse(substitute(pkg))
  if (!requireNamespace(pkg, quietly = TRUE)) {
    lockfile <- renv::lockfile_read()
    if (!pkg %in% names(lockfile$Packages)) {
      stop(sprintf("Package '%s' is not in renv.lock", pkg))
    }
    renv::install(pkg, prompt = FALSE)
  }
  library(pkg, character.only = TRUE)
}

# Hook into loadNamespace to auto-install when using pkg::func()
local({
  original <- base::loadNamespace
  wrapper <- function(package, ...) {
    pkg_name <- as.character(package)
    if (!isNamespaceLoaded(pkg_name) && !nzchar(system.file(package = pkg_name))) {
      lockfile <- renv::lockfile_read()
      if (pkg_name %in% names(lockfile$Packages)) {
        renv::install(pkg_name, prompt = FALSE)
      }
    }
    original(package, ...)
  }
  unlockBinding("loadNamespace", baseenv())
  assign("loadNamespace", wrapper, baseenv())
  lockBinding("loadNamespace", baseenv())
})

# if (is_domino) {
data_path <- Sys.getenv("DOMINO_DATASETS_DIR")
output_path <- file.path(data_path, "analysis-results", Sys.getenv("DOMINO_STARTING_USERNAME"))
artifacts_path <- file.path(Sys.getenv("DOMINO_ARTIFACTS_DIR"), Sys.getenv("DOMINO_STARTING_USERNAME"))
fit_output_timestamp <- FALSE

lazy_lib(conflicted)

conflicts_prefer(
  dplyr::filter,
  dplyr::lag,
  posterior::sd,
  posterior::mad,
  rlang::set_names,
  purrr::flatten_dbl,
)

options(yaml.eval.expr = TRUE)

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

init_project <- function(output_path = output_path, artifacts_path = artifacts_path) {
  lazy_lib(magrittr)
  lazy_lib(tidyverse)
  lazy_lib(rlang)
  lazy_lib(targets)
  lazy_lib(tarchetypes)
  lazy_lib(stantargets)
  lazy_lib(crew)
  lazy_lib(autometric)
  lazy_lib(here)
  lazy_lib(cmdstanr)
  lazy_lib(posterior)
  lazy_lib(tidybayes)
  lazy_lib(qs2)
  lazy_lib(recipes)

  source(here("r", "util.R"))
  source(here("r", "priors.R"))
  source(here("r", "posterior.R"))
  source(here("r", "prepare_analysis_data.R"))
  source(here("r", "initializers.R"))
  source(here("r", "accuracy.R"))
  source(here("r", "state_space.R"))
  source(here("r", "plot_functions.R"))
  source(here("r", "parquet_draws.R"))

  source(here("r", "sclc", "priors.R"))
  source(here("r", "sclc", "prepare_analysis_data.R"))
  source(here("r", "sclc", "accuracy.R"))
  source(here("r", "sclc", "plot_functions.R"))

  source(here("r", "targets_tidyselect.R"))
}
