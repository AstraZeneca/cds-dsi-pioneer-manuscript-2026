
if (file.exists("~/.Rprofile")) source("~/.Rprofile")

source("renv/activate.R")

# azcore::azcore_module_load("GLPK")
# #azcore::azcore_module_load("CMake")
# azcore::azcore_bundle_dynload("libglpk.so.40")
# azcore::azcore_module_load("git/2.38.1-GCCcore-10.3.0-nodocs")

library(conflicted)

conflicts_prefer(
  dplyr::filter, dplyr::lag,
  posterior::sd, posterior::mad,
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
