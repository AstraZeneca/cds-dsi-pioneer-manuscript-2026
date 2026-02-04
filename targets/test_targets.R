# library(rlang)
# library(targets)
# library(tarchetypes)
# library(crew)
# library(here)
# library(qs2)
#
# tm <- tar_map(
#   values = list(x = 1:5),
#   tar_target(t1, data.frame(a = x)),
#   tar_target(t2, data.frame(b = x))
# )
#
# all_t <- list(
#   tar_target(t3, 100)
# )
#
# list(
#   all_t,
#   tm,
#   tar_combine(ct, tm[["t1"]], command = dplyr::bind_rows(!!!.x))
# )

library(magrittr)
library(tidyverse)
library(rlang)
library(targets)
library(tarchetypes)
library(crew)
library(autometric)
library(here)
library(cmdstanr)
library(posterior)
library(tidybayes)
library(qs2)

source(here("r", "util.R"))
source(here("r", "priors.R"))
source(here("r", "posterior.R"))
source(here("r", "prepare_analysis_data.R"))
source(here("r", "initializers.R"))
source(here("r", "accuracy.R"))
source(here("r", "state_space.R"))

source(here("r", "sclc", "priors.R"))
source(here("r", "sclc", "prepare_analysis_data.R"))

sclc_output_path <- file.path(output_path, "sclc")
sclc_artifacts_path <- file.path(artifacts_path, "sclc")

# fs::dir_create(output_path) # Data folder
# fs::dir_create(artifacts_path) # Artifacts folder
fs::dir_create(file.path(sclc_output_path, "fit"))
fs::dir_create(file.path(sclc_artifacts_path, "models"))
fs::dir_create(file.path(sclc_artifacts_path, "crew_logs"))

# We're using different controllers for different types of targets depending on the resources they end up using.

# crew_options <- crew_options_local(log_directory = file.path(artifacts_path, "crew_logs"))
controller_default <- crew_controller_local(
  name = "default",
  workers = 20,
)
controller_intense <- crew_controller_local(
  name = "intense",
  workers = 10,
)
controller_sampling <- crew_controller_local(
  name = "sampling",
  workers = 22,
)

tar_option_set(
  packages = c(
    "magrittr",
    "tidyverse",
    "rlang",
    "here",
    "targets",
    "cmdstanr",
    "tidybayes",
    "posterior",
    "priorsense",
    "loo"
  ),
  controller = crew_controller_group(
    controller_default,
    controller_intense,
    controller_sampling
  ),
  resources = tar_resources(crew = tar_resources_crew(controller = "default")),
  memory = "transient",
  garbage_collection = TRUE,
  format = "qs",
  error = "continue",
  seed = 26091468
)

# Some global settings. Some of the Boolean ones are used with tar_skip() to skip some targets depending on the project. Still reports errors at runtime though.

cat("TAR_PROJECT =", Sys.getenv("TAR_PROJECT"), "\n")
cat("User =", Sys.getenv("DOMINO_STARTING_USERNAME"), "\n")
cat("Targets store =", tar_config_get("store"), "\n")
cat("RStudio = ", rstudioapi::isAvailable() || rstudioapi::isJob(), "\n")
cat("data_path =", data_path, "\n")
cat("output_path =", sclc_output_path, "\n")
cat("artifacts_path =", sclc_artifacts_path, "\n")

tumor_ssls <- tar_map(
  tibble(
    fit_data = c(FALSE, TRUE),
    base_name = c("prior_tumor_ssls", "tumor_ssls"),
    iter_warmup = c(300, 500),
    iter_sampling = 500,
    fit_type = c("prior", "posterior"),
  ),
  names = "fit_type",

  tar_target(t1, tibble(fit_type))
)

sclc_targets <- lst(
  tumor_ssls,

  tar_combine(
    all_tumor_ssls_rates_rvar,
    tumor_ssls[["t1"]],
    command = dplyr::bind_rows(!!!.x)
  )
  # # tar_combine(all_tumor_ssls_patient_rates_rvar, tar_select_targets(tumor_ssls, starts_with("tumor_ssls_patient_rates")), command = bind_rows(!!!.x)),
  # tar_combine(all_tumor_ssls_patient_rates_rvar, tumor_ssls[["tumor_ssls_rates"]], command = bind_rows(!!!.x)),
  # tar_combine(all_tumor_ssls_noise_sd_rvar, tar_select_targets(tumor_ssls, starts_with("tumor_ssls_noise_sd_rvar"))),
  # tar_combine(all_tumor_ssls_log_growth_lag_rvar, tar_select_targets(tumor_ssls, starts_with("tumor_ssls_log_growth_lag"))),
  # # tar_combine(all_tumor_ssls_patient_log_growth_lag, tar_select_targets(tumor_ssls, starts_with("tumor_ssls_patient_log_growth_lag"))),
  # tar_combine(all_tumor_ssls_patient_decrease_prop_rvar, tar_select_targets(tumor_ssls, starts_with("tumor_ssls_patient_decrease_prop"))),
  # tar_combine(all_tumor_ssls_patient_states_rvar, tar_select_targets(tumor_ssls, starts_with("tumor_ssls_patient_states_rvar")), storage = "worker"),
)

# sclc_targets <- list(
#   # sclc_targets,
#   tumor_ssls,
#
#   tar_combine(all_tumor_ssls_rates_rvar, tumor_ssls[["tumor_ssls_rates_rvar"]], dplyr::bind_rows(!!!.x))
# )
