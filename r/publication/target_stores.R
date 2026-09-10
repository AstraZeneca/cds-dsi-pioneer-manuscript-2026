pub_store_sclc <- file.path(output_path, "publication", Sys.getenv("TAR_RUN", "main"), "_targets")
pub_store_crc  <- file.path(output_path, "publication", "crc",  "_targets")

# pub_store is set by each subsite's _setup.R via pub_store_sclc or pub_store_crc
if (!exists("pub_store")) {
  pub_store <- pub_store_sclc
}

# LFO cross-validation store. The LFO targets often live in a SEPARATE run from
# the main posterior results (e.g. main results in `gompertz`, LFO in
# `gompertz-lfo`), because LFO is fit as its own pipeline. So it gets its own env
# var: LFO_TAR_RUN takes precedence, then TAR_RUN, then the literal "lfo". This
# lets a single render point pub_store at one run and lfo_store at another:
#   TAR_RUN=gompertz LFO_TAR_RUN=gompertz-lfo quarto render ...
# The dual-disease LFO targets are read back by name with the disease suffix
# (pub_disease_suffix).
if (!exists("lfo_store")) {
  lfo_run <- Sys.getenv("LFO_TAR_RUN", Sys.getenv("TAR_RUN", "lfo"))
  lfo_store <- file.path(output_path, "publication", lfo_run, "_targets")
}
