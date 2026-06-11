pub_disease_suffix <- ""  # blank: SCLC store predates tar_map disease suffix

here::i_am("quarto/publication/sclc/_setup.R")

source(here::here("quarto", "_shared", "_setup.R"))

source(here::here("r", "sclc", "plot_functions.R"))
source(here::here("r", "sclc", "table_functions.R"))
source(here::here("r", "sclc", "accuracy.R"))

source(here::here("r", "publication", "target_stores.R"))
pub_store <- pub_store_sclc
source(here::here("quarto", "publication", "_setup_common.R"))
