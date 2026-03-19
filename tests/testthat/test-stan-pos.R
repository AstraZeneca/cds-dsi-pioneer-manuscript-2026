library(testthat)

# Build padded input matrices for the batched test_pos_all.stan harness.
# group_sizes: list of integer vectors (one per case)
# flat_data:   list of integer vectors (one per case)
build_pos_data <- function(group_sizes_list, flat_data_list) {
  n_cases    <- length(group_sizes_list)
  max_groups <- max(sapply(group_sizes_list, length))
  max_size   <- max(sapply(flat_data_list, length))

  gs_mat <- matrix(0L, nrow = n_cases, ncol = max_groups)
  fd_mat <- matrix(0L, nrow = n_cases, ncol = max_size)
  n_grp  <- integer(n_cases)
  n_flat <- integer(n_cases)

  for (i in seq_len(n_cases)) {
    ng <- length(group_sizes_list[[i]])
    nf <- length(flat_data_list[[i]])
    gs_mat[i, seq_len(ng)] <- group_sizes_list[[i]]
    fd_mat[i, seq_len(nf)] <- flat_data_list[[i]]
    n_grp[i]  <- ng
    n_flat[i] <- nf
  }

  list(
    N_CASES    = n_cases,
    MAX_GROUPS = max_groups,
    MAX_SIZE   = max_size,
    group_sizes = gs_mat,
    flat_data   = fd_mat,
    n_groups    = n_grp,
    n_flat      = n_flat
  )
}

test_that("create_pos and get_pos_size match R oracle for multiple cases", {
  # Case 1: sizes c(2, 3, 1) → pos = c(1, 3, 6, 7)
  # Case 2: sizes c(0, 4)    → pos = c(1, 1, 5)
  # Case 3: sizes c(5)       → pos = c(1, 6)
  group_sizes <- list(c(2L, 3L, 1L), c(0L, 4L), c(5L))
  flat_data   <- list(c(10L, 20L, 30L, 40L, 50L, 60L),
                      c(1L, 2L, 3L, 4L),
                      c(7L, 5L, 9L, 3L, 6L))

  data <- build_pos_data(group_sizes, flat_data)
  fit  <- test_stan_function(
    here::here("tests", "testthat", "stan", "test_pos_all.stan"), data
  )
  d <- posterior::as_draws_df(fit$draws())

  # Case 1 pos: c(1, 3, 6, 7)
  expected_pos1 <- r_create_pos(c(2L, 3L, 1L))
  for (i in seq_along(expected_pos1)) {
    expect_equal(get_stan_val(d, "pos_out", 1, i), expected_pos1[i],
      label = paste0("case1 pos_out[", i, "]"))
  }

  # Case 2 pos: c(1, 1, 5)
  expected_pos2 <- r_create_pos(c(0L, 4L))
  for (i in seq_along(expected_pos2)) {
    expect_equal(get_stan_val(d, "pos_out", 2, i), expected_pos2[i],
      label = paste0("case2 pos_out[", i, "]"))
  }

  # pos_size for case 1: c(2, 3, 1)
  expect_equal(get_stan_val(d, "pos_size_out", 1, 1), 2L)
  expect_equal(get_stan_val(d, "pos_size_out", 1, 2), 3L)
  expect_equal(get_stan_val(d, "pos_size_out", 1, 3), 1L)

  # pos_size for case 2: group 1 empty (0), group 2 has 4
  expect_equal(get_stan_val(d, "pos_size_out", 2, 1), 0L)
  expect_equal(get_stan_val(d, "pos_size_out", 2, 2), 4L)
})

test_that("get_max_pos and get_min_pos per group", {
  # Case: sizes c(2, 3), flat = c(5, 2, 8, 1, 9)
  # group 1: {5, 2} → min=2, max=5
  # group 2: {8, 1, 9} → min=1, max=9
  group_sizes <- list(c(2L, 3L))
  flat_data   <- list(c(5L, 2L, 8L, 1L, 9L))

  data <- build_pos_data(group_sizes, flat_data)
  fit  <- test_stan_function(
    here::here("tests", "testthat", "stan", "test_pos_all.stan"), data
  )
  d <- posterior::as_draws_df(fit$draws())

  expect_equal(get_stan_val(d, "min_out", 1, 1), 2L)
  expect_equal(get_stan_val(d, "max_out", 1, 1), 5L)
  expect_equal(get_stan_val(d, "min_out", 1, 2), 1L)
  expect_equal(get_stan_val(d, "max_out", 1, 2), 9L)
})

test_that("get_int (last element of each group via last_int_out)", {
  # get_int(x, pos, i, get_pos_size(pos, i)) = last element of group i
  # sizes c(3, 2), flat = c(10, 20, 30, 40, 50)
  # group 1 last: 30, group 2 last: 50
  group_sizes <- list(c(3L, 2L))
  flat_data   <- list(c(10L, 20L, 30L, 40L, 50L))

  data <- build_pos_data(group_sizes, flat_data)
  fit  <- test_stan_function(
    here::here("tests", "testthat", "stan", "test_pos_all.stan"), data
  )
  d <- posterior::as_draws_df(fit$draws())

  expect_equal(get_stan_val(d, "last_int_out", 1, 1), 30L)
  expect_equal(get_stan_val(d, "last_int_out", 1, 2), 50L)
})

test_that("all-zero sizes: pos is all 1s, group sizes are all 0", {
  group_sizes <- list(c(0L, 0L, 0L))
  flat_data   <- list(c(999L))  # unused, but MAX_SIZE must be >= 1

  data <- build_pos_data(group_sizes, flat_data)
  fit  <- test_stan_function(
    here::here("tests", "testthat", "stan", "test_pos_all.stan"), data
  )
  d <- posterior::as_draws_df(fit$draws())

  # pos = c(1, 1, 1, 1)
  for (i in 1:4) {
    expect_equal(get_stan_val(d, "pos_out", 1, i), 1L,
      label = paste0("all-zero pos_out[", i, "]"))
  }
  for (i in 1:3) {
    expect_equal(get_stan_val(d, "pos_size_out", 1, i), 0L,
      label = paste0("all-zero pos_size_out[", i, "]"))
  }
})
