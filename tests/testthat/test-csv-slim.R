# Tests for project_cmdstan_csv() in r/publication/csv_slim.R

source(here::here("r", "publication", "csv_slim.R"))

# Helper: build a tiny synthetic CmdStan CSV string
make_fake_csv <- function(n_rows = 3, n_cols = 6) {
  header  <- paste(c("lp__", "accept_stat__", paste0("param_", seq_len(n_cols - 2))), collapse = ",")
  rows    <- replicate(n_rows, paste(round(runif(n_cols), 4), collapse = ","))
  comment <- c("# stan_version_major = 2", "# num_samples = 3", "# num_warmup = 0")
  footer  <- c("# Elapsed Time: 1 seconds (Warm-up)", "#               2 seconds (Sampling)")
  paste(c(comment, header, rows, footer), collapse = "\n")
}

write_fake_csv <- function(content, path) {
  writeLines(content, path)
  path
}

test_that("project_cmdstan_csv keeps comment lines verbatim", {
  tmp_in  <- withr::local_tempfile(fileext = ".csv")
  tmp_out <- withr::local_tempfile(fileext = ".csv")
  write_fake_csv(make_fake_csv(), tmp_in)

  project_cmdstan_csv(tmp_in, tmp_out, keep_col_names = c("lp__", "param_1"))

  out_lines <- readLines(tmp_out)
  comment_lines <- out_lines[startsWith(out_lines, "#")]
  expect_length(comment_lines, 5L)   # 3 header + 2 footer comments
  expect_true(all(startsWith(comment_lines, "#")))
})

test_that("project_cmdstan_csv projects to exactly the kept columns", {
  tmp_in  <- withr::local_tempfile(fileext = ".csv")
  tmp_out <- withr::local_tempfile(fileext = ".csv")
  write_fake_csv(make_fake_csv(n_rows = 3, n_cols = 6), tmp_in)

  project_cmdstan_csv(tmp_in, tmp_out, keep_col_names = c("lp__", "param_2", "param_4"))

  out_lines <- readLines(tmp_out)
  data_lines <- out_lines[!startsWith(out_lines, "#")]
  # first data line is header
  header_fields <- strsplit(data_lines[[1]], ",")[[1]]
  expect_equal(header_fields, c("lp__", "param_2", "param_4"))
  # remaining lines have same field count
  row_counts <- lengths(strsplit(data_lines[-1], ","))
  expect_true(all(row_counts == 3L))
})

test_that("project_cmdstan_csv errors when keep_col_names has unknown columns", {
  tmp_in  <- withr::local_tempfile(fileext = ".csv")
  tmp_out <- withr::local_tempfile(fileext = ".csv")
  write_fake_csv(make_fake_csv(), tmp_in)

  expect_error(
    project_cmdstan_csv(tmp_in, tmp_out, keep_col_names = c("lp__", "nonexistent_col")),
    "not found in CSV header"
  )
  expect_false(file.exists(tmp_out))
})

test_that("project_cmdstan_csv preserves all data values for kept columns", {
  set.seed(42)
  tmp_in  <- withr::local_tempfile(fileext = ".csv")
  tmp_out <- withr::local_tempfile(fileext = ".csv")
  # 4 cols: lp__, accept_stat__, param_1, param_2
  rows <- list(c(1.1, 0.9, 2.2, 3.3), c(4.4, 0.8, 5.5, 6.6))
  header  <- "lp__,accept_stat__,param_1,param_2"
  content <- paste(c("# num_samples = 2",
                     header,
                     paste(rows[[1]], collapse = ","),
                     paste(rows[[2]], collapse = ",")),
                   collapse = "\n")
  writeLines(content, tmp_in)

  project_cmdstan_csv(tmp_in, tmp_out, keep_col_names = c("lp__", "param_2"))

  out_lines <- readLines(tmp_out)
  data_lines <- out_lines[!startsWith(out_lines, "#")]
  row1 <- as.numeric(strsplit(data_lines[[2]], ",")[[1]])  # data_lines[[1]] is header
  row2 <- as.numeric(strsplit(data_lines[[3]], ",")[[1]])
  expect_equal(row1, c(1.1, 3.3))
  expect_equal(row2, c(4.4, 6.6))
})

test_that("slim_cmdstan_csv_files writes to output_dir with same filename", {
  tmp_in_dir  <- withr::local_tempdir()
  tmp_out_dir <- withr::local_tempdir()

  # Write a fake CSV with a CmdStan-style filename
  fake_csv <- file.path(tmp_in_dir, "model-202601010000-1-abc123.csv")
  content <- paste(c(
    "# num_samples = 1",
    "lp__,accept_stat__,param_1,param_2",
    "1.1,0.9,2.2,3.3"
  ), collapse = "\n")
  writeLines(content, fake_csv)

  result <- slim_cmdstan_csv_files(
    input_csv    = fake_csv,
    output_dir   = tmp_out_dir,
    keep_col_names = c("lp__", "accept_stat__", "param_1")
  )

  expect_equal(basename(result), "model-202601010000-1-abc123.csv")
  expect_true(file.exists(result))
  out_lines <- readLines(result)
  data_lines <- out_lines[!startsWith(out_lines, "#")]
  expect_equal(strsplit(data_lines[[1]], ",")[[1]], c("lp__", "accept_stat__", "param_1"))
})

test_that("build_crc_keep_col_names matches expected publication column families", {
  tmp_csv <- withr::local_tempfile(fileext = ".csv")

  # Synthetic header covering the key pattern families
  # Include index suffixes like real CmdStan output
  header_cols <- c(
    # sampler (exact match via CRC_SLIM_KEEP_COLS)
    "lp__", "accept_stat__", "stepsize__", "treedepth__",
    "n_leapfrog__", "divergent__", "energy__",
    # exact named cols
    "measure_sd_sld", "rep_patient_log_sld.1", "forecast_patient_log_sld.2",
    "rep_recist.3", "forecast_obs_recist.4", "recist_confusion_matrix",
    # _pop$ family
    "tr_loc_pop", "frac_logit_loc_pop",
    # pop_ family
    "pop_log_decrease_rate.1",
    # _sd_level_ family
    "tr_sd_level_intercept.1.2",
    # time_invariant/time_varying_coef
    "time_invariant_coef_01.1", "time_varying_coef_12.2",
    # patient family
    "tr_loc_patient.1", "frac_log_growth_patient.2",
    # patient_log rate
    "patient_log_decrease_rate.1.3", "patient_log_growth_rate.2.4",
    # spop/sample KM
    "sample_target_km_est.1.200", "spop_pfs_km_est.2.100",
    # os endpoints
    "sample_os.1", "spop_os_censored.2",
    # should NOT be kept
    "log_cond_surv_01.1.200", "states_full_grid.1.1.218",
    "ms_time_varying_covar_01.1.200"
  )
  content <- paste(c("# num_samples = 1",
                     paste(header_cols, collapse = ","),
                     paste(rep("0.1", length(header_cols)), collapse = ",")),
                   collapse = "\n")
  writeLines(content, tmp_csv)

  kept <- build_crc_keep_col_names(tmp_csv)

  # All expected columns kept
  expect_true("lp__" %in% kept)
  expect_true("measure_sd_sld" %in% kept)
  expect_true("rep_patient_log_sld.1" %in% kept)
  expect_true("tr_loc_pop" %in% kept)
  expect_true("pop_log_decrease_rate.1" %in% kept)
  expect_true("tr_sd_level_intercept.1.2" %in% kept)
  expect_true("time_invariant_coef_01.1" %in% kept)
  expect_true("tr_loc_patient.1" %in% kept)
  expect_true("patient_log_decrease_rate.1.3" %in% kept)
  expect_true("sample_target_km_est.1.200" %in% kept)
  expect_true("sample_os.1" %in% kept)
  expect_true("spop_os_censored.2" %in% kept)

  # Columns that should NOT be kept
  expect_false("log_cond_surv_01.1.200" %in% kept)
  expect_false("states_full_grid.1.1.218" %in% kept)
  expect_false("ms_time_varying_covar_01.1.200" %in% kept)
})
