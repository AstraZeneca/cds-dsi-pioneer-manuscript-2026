pub_disease_suffix <- "sclc"
lfo_disease_suffix <- "sclc"  # LFO targets live inside the disease map → always suffixed

here::i_am("quarto/publication/sclc/_setup.R")

source(here::here("quarto", "_shared", "_setup.R"))

source(here::here("r", "accuracy.R"))  # get_oos_confusion_marix + OOS plot/table helpers

source(here::here("r", "publication", "target_stores.R"))
pub_store <- pub_store_sclc
source(here::here("quarto", "publication", "_setup_common.R"))
