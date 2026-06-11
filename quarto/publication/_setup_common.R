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
