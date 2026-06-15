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
  obj <- targets::tar_read_raw(full_name, store = store)
  .repair_rvar_dims(obj)
}

# lfo_tar_read() reads an LFO target from lfo_store, always appending the disease
# name (lfo_disease_suffix). Unlike pub_disease_suffix (blank for the legacy SCLC
# store), the LFO targets live inside the disease tar_map, so SCLC names are
# suffixed _sclc just like CRC is _crc. Set lfo_disease_suffix in each subsite's
# _setup.R before sourcing this file.
lfo_tar_read <- function(name, store = lfo_store) {
  full_name <- paste0(name, "_", lfo_disease_suffix)
  obj <- targets::tar_read_raw(full_name, store = store)
  .repair_rvar_dims(obj)
}

# lfo_tar_read_pattern() reads all branches of an LFO pattern target and binds
# them with bind_rows(). tar_read_raw() on the parent name calls targets'
# tar_vec_c() which uses vctrs::vec_c() — this fails on rvar columns when branch
# dimensions differ. Reading branches individually and binding avoids the issue.
lfo_tar_read_pattern <- function(name, store = lfo_store) {
  full_name <- paste0(name, "_", lfo_disease_suffix)
  branches <- targets::tar_meta(store = store, fields = c("name", "children")) |>
    dplyr::filter(name == full_name) |>
    dplyr::pull(children) |>
    unlist()
  purrr::map(branches, \(b) {
    obj <- targets::tar_read_raw(b, store = store)
    .repair_rvar_dims(obj)
  }) |>
    dplyr::bind_rows()
}

# qs2 drops the dim attribute on rvar wrappers during deserialisation (length=0).
# Repair by reconstructing from the draws matrix, which is always correct.
.repair_rvar_dims <- function(obj) {
  if (!is.data.frame(obj)) return(obj)
  rvar_cols <- which(vapply(obj, inherits, logical(1), "rvar"))
  for (j in rvar_cols) {
    rv <- obj[[j]]
    if (length(rv) == 0L) {
      dm <- posterior::draws_of(rv)
      n  <- ncol(dm)
      if (n > 0L) {
        ndraws <- nrow(dm)
        nch    <- posterior::nchains(rv)
        obj[[j]] <- posterior::rvar(
          array(dm, dim = c(ndraws, n)),
          nchains = nch
        )
      }
    }
  }
  obj
}
