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
