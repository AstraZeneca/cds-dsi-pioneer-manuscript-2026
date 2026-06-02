# quarto/publication/_setup_common.R
# pub_tar_read() reads a target from pub_store, appending pub_disease_suffix
# when non-empty. Set pub_disease_suffix before sourcing this file.
#
# pub_disease_suffix <- ""     → reads unsuffixed name (legacy SCLC store)
# pub_disease_suffix <- "crc"  → reads <name>_crc

pub_tar_read <- function(name, store = pub_store) {
  full_name <- if (nchar(pub_disease_suffix) > 0) {
    paste0(name, "_", pub_disease_suffix)
  } else {
    name
  }
  targets::tar_read_raw(full_name, store = store)
}
