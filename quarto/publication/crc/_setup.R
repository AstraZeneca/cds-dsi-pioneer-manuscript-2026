pub_disease_suffix <- "crc"

here::i_am("quarto/publication/crc/_setup.R")

source(here::here("quarto", "_shared", "_setup.R"))

source(here::here("r", "sclc", "plot_functions.R"))
source(here::here("r", "sclc", "table_functions.R"))
source(here::here("r", "sclc", "accuracy.R"))

source(here::here("r", "publication", "target_stores.R"))
pub_store <- pub_store_crc
source(here::here("quarto", "publication", "_setup_common.R"))
