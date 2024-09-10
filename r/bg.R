library(targets)

tar_make(
  c(
    "bootstrap_confirmed_resp_pfs_res", 
    "pfs_conf_resp_bootstrap_cr_median_pfs", 
    "pfs_conf_resp_bootstrap_pfs_median_pfs", 
    "fixed_bootstrap_cr_median_pfs",
    "fixed_bootstrap_pfs_median_pfs",
    "bootstrap_cr_median_pfs_draws",
    "pfs_conf_resp_bootstrap_param"
  ),
  use_crew = TRUE
)
