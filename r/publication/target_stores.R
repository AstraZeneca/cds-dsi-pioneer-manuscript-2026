pub_store_sclc <- file.path(output_path, "publication", "main", "_targets")
pub_store_crc  <- file.path(output_path, "publication", "crc",  "_targets")

# pub_store is set by each subsite's _setup.R via pub_store_sclc or pub_store_crc
if (!exists("pub_store")) {
  pub_store <- pub_store_sclc
}
