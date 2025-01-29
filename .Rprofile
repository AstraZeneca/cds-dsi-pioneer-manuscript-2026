
if (file.exists("~/.Rprofile")) source("~/.Rprofile")

source("renv/activate.R")

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

is_domino <- !is.na(Sys.getenv("IS_DOMINO")) && Sys.getenv("IS_DOMINO") == "true"

if (is_domino) {
  output_path <- "/mnt/artifacts"
  data_path <- "/mnt/data/endometrial-to-lung"
} else {
  user <- Sys.info()["user"]
  output_path <- file.path("/scratch", user, "pioneer")
  data_path <- "/scratch/ewfteams/dpo0160"
  
  options(
    renv.config.external.library = "/opt/scp/services/azcore/coreutils/0.1.0/libraries/R/4.3.1",
    renv.config.ignored.packages = c("azcore", "rseed")
  )
}
