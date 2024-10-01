
if (file.exists("~/.Rprofile")) source("~/.Rprofile")

source("renv/activate.R")

# azcore::azcore_module_load("GLPK")
# #azcore::azcore_module_load("CMake")
# azcore::azcore_bundle_dynload("libglpk.so.40")
azcore::azcore_module_load("git/2.38.1-GCCcore-10.3.0-nodocs")

library(conflicted)

conflicts_prefer(
  dplyr::filter, dplyr::lag,
  posterior::sd, posterior::mad,
  rlang::set_names,
)

