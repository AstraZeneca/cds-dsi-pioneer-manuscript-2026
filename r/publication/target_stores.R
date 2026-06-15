pub_store_sclc <- file.path(output_path, "publication", Sys.getenv("TAR_RUN", "main"), "_targets")
pub_store_crc  <- file.path(output_path, "publication", "crc",  "_targets")

# pub_store is set by each subsite's _setup.R via pub_store_sclc or pub_store_crc
if (!exists("pub_store")) {
  pub_store <- pub_store_sclc
}

# LFO cross-validation store. Resolved per run from TAR_RUN, mirroring how the
# pipeline store itself is selected — the dual-disease LFO targets are built into
# whichever store the run targets (publication/lfo for SCLC, a crc-suffixed run
# for CRC), and read back by name with the disease suffix (pub_disease_suffix).
if (!exists("lfo_store")) {
  lfo_store <- file.path(output_path, "publication", Sys.getenv("TAR_RUN", "lfo"), "_targets")
}
