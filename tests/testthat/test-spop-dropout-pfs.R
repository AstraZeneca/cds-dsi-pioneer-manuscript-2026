# Invariants for the spop dropout PFS-from-OS convention (spec 2026-05-30).
# Reads the regenerated endpoint draws from the publication store.
test_that("spop dropout PFS-from-OS invariants hold", {
  store <- file.path("/mnt/data/analysis-results",
                     Sys.getenv("DOMINO_STARTING_USERNAME"),
                     "publication/main/_targets")
  skip_if_not(dir.exists(store), "publication store not present")

  fit <- targets::tar_read(tumor_ssls_res_posterior, store = store)
  ad  <- targets::tar_read(all_analysis_data, store = store)
  np  <- nrow(ad)

  pull <- function(v) {
    m <- posterior::as_draws_matrix(fit$draws(variables = v))
    m[, match(paste0(v, "[", seq_len(np), "]"), colnames(m)), drop = FALSE]
  }
  s_pfs  <- pull("spop_pfs");  s_rc  <- pull("spop_right_censored")
  s_os   <- pull("spop_os");   s_osc <- pull("spop_os_censored")
  drop   <- pull("spop_is_dropout")

  # (1) PFS <= OS pointwise (every draw, every patient)
  expect_true(all(s_pfs <= s_os),
              info = "PFS must never exceed OS")

  # (2) For dropouts: PFS == OS and PFS-censoring == OS-censoring (graft identity)
  dmask <- drop == 1L
  expect_true(all(s_pfs[dmask] == s_os[dmask]),
              info = "dropout PFS must equal OS by construction")
  expect_true(all(s_rc[dmask] == s_osc[dmask]),
              info = "dropout PFS censoring must equal OS censoring")
})
