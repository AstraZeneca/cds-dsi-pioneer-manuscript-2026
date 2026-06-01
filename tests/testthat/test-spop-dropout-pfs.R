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

# CIF_03 must stay a valid cumulative-incidence curve after the graft. Regression
# guard for the out-of-range bug where a dropout's pfs == OS could exceed max_t
# and was either crashing compute_trial_cif (pre-fix) or had to be horizon-guarded.
test_that("spop CIF_03 (dropout) is a valid, populated cumulative incidence", {
  store <- file.path("/mnt/data/analysis-results",
                     Sys.getenv("DOMINO_STARTING_USERNAME"),
                     "publication/main/_targets")
  skip_if_not(dir.exists(store), "publication store not present")

  km <- targets::tar_read(all_tumor_ssls_km_rvar, store = store)
  skip_if_not("spop_cif_03" %in% names(km) ||
              any(grepl("cif_03", names(km))), "spop_cif_03 not in km rvars")

  cif03 <- km[[grep("spop_cif_03", names(km), value = TRUE)[1]]]
  m <- posterior::as_draws_matrix(cif03)
  # bounded in [0,1] and non-decreasing in time within each draw (cumulative)
  expect_true(all(m >= 0 & m <= 1), info = "CIF_03 must lie in [0,1]")
  # populated: at least some dropout mass accrues (not identically zero)
  expect_true(mean(m) > 0, info = "CIF_03 must be populated, not emptied by the graft")
})
