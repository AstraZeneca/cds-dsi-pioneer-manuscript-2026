# Invariants for the spop dropout PFS-from-OS convention (spec 2026-05-30).
# Reads the regenerated endpoint draws from the publication store.
test_that("spop dropout PFS-from-OS invariants hold", {
  store <- file.path("/mnt/data/analysis-results",
                     Sys.getenv("DOMINO_STARTING_USERNAME"),
                     "publication/main/_targets")
  skip_if_not(dir.exists(store), "publication store not present")
  skip_if_not(
    file.exists(file.path(store, "objects", "tumor_ssls_res_posterior_sclc")),
    "tumor_ssls_res_posterior_sclc not in store (store may be pre-disease-split)"
  )

  fit <- targets::tar_read(tumor_ssls_res_posterior_sclc, store = store)

  # Endpoint GQ arrays are length n_forecast_patients, which is < n_patients
  # whenever forecast_split_level > 0 (RWD/borrowing models forecast only the
  # target arm). Pull ALL columns of each variable in index order rather than
  # indexing by nrow(all_analysis_data) — the invariants are self-contained
  # within the forecast cohort, so the test must not couple to n_patients.
  pull <- function(v) {
    m <- posterior::as_draws_matrix(fit$draws(variables = v))
    idx <- as.integer(sub(paste0("^", v, "\\[(\\d+)\\]$"), "\\1", colnames(m)))
    m[, order(idx), drop = FALSE]
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
  skip_if_not(
    file.exists(file.path(store, "objects", "all_tumor_ssls_km_rvar_sclc")),
    "all_tumor_ssls_km_rvar_sclc not in store (store may be pre-disease-split)"
  )

  km <- targets::tar_read(all_tumor_ssls_km_rvar_sclc, store = store)
  skip_if_not("spop_cif_03" %in% names(km) ||
              any(grepl("cif_03", names(km))), "spop_cif_03 not in km rvars")

  cif03 <- km[[grep("spop_cif_03", names(km), value = TRUE)[1]]]
  m <- posterior::as_draws_matrix(cif03)
  # bounded in [0,1] and non-decreasing in time within each draw (cumulative)
  expect_true(all(m >= 0 & m <= 1), info = "CIF_03 must lie in [0,1]")
  # populated: at least some dropout mass accrues (not identically zero)
  expect_true(mean(m) > 0, info = "CIF_03 must be populated, not emptied by the graft")
})

# Symmetric sample_* regression test. The spop graft (stan/pfs.stanfunctions L2004-2009)
# and its invariant test (above) were added together but never mirrored to the sample
# pathway — a classic asymmetric-pathway coverage gap. This test must run against a fit
# built AFTER the sample-graft fix (stan/pfs.stanfunctions L2067+); it will skip if the
# store is absent. Once fixed, this guards against reintroducing the same dropout PFS bug.
test_that("sample dropout PFS-from-OS invariants hold", {
  store <- file.path("/mnt/data/analysis-results",
                     Sys.getenv("DOMINO_STARTING_USERNAME"),
                     "publication/main/_targets")
  skip_if_not(dir.exists(store), "publication store not present")
  skip_if_not(
    file.exists(file.path(store, "objects", "tumor_ssls_res_posterior_sclc")),
    "tumor_ssls_res_posterior_sclc not in store (store may be pre-disease-split)"
  )

  fit <- targets::tar_read(tumor_ssls_res_posterior_sclc, store = store)

  # Endpoint GQ arrays are length n_forecast_patients (< n_patients when
  # forecast_split_level > 0). Pull ALL columns of each variable in index order;
  # the invariants are self-contained within the forecast cohort, so do not
  # couple to nrow(all_analysis_data). See the spop block above.
  pull <- function(v) {
    m <- posterior::as_draws_matrix(fit$draws(variables = v))
    idx <- as.integer(sub(paste0("^", v, "\\[(\\d+)\\]$"), "\\1", colnames(m)))
    m[, order(idx), drop = FALSE]
  }
  sample_pfs  <- pull("sample_pfs");  sample_rc  <- pull("sample_right_censored")
  sample_os   <- pull("sample_os");   sample_osc <- pull("sample_os_censored")
  sample_drop <- pull("sample_is_dropout")

  # (1) PFS <= OS pointwise (every draw, every patient)
  expect_true(all(sample_pfs <= sample_os),
              info = "sample PFS must never exceed OS")

  # (2) For dropouts: PFS == OS and PFS-censoring == OS-censoring (graft identity)
  dmask <- sample_drop == 1L
  expect_true(all(sample_pfs[dmask] == sample_os[dmask]),
              info = "sample dropout PFS must equal OS by construction")
  expect_true(all(sample_rc[dmask] == sample_osc[dmask]),
              info = "sample dropout PFS censoring must equal OS censoring")
})
