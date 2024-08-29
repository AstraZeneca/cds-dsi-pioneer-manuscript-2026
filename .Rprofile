
if (file.exists("~/.Rprofile")) source("~/.Rprofile")

source("renv/activate.R")

library(conflicted)

conflicts_prefer(
  dplyr::filter, 
  posterior::sd, posterior::mad
)

azcore::azcore_module_load("GLPK")
#azcore::azcore_module_load("CMake")
azcore::azcore_bundle_dynload("libglpk.so.40")
azcore::azcore_module_load("git/2.43.3-GCCcore-10.3.0")
