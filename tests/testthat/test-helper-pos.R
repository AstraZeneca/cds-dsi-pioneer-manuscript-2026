library(testthat)

test_that("r_create_pos: basic cumulative positions", {
  pos <- r_create_pos(c(2L, 3L, 1L))
  # groups of size 2, 3, 1 → pos = c(1, 3, 6, 7)
  expect_equal(pos, c(1L, 3L, 6L, 7L))
  expect_length(pos, 4L)
})

test_that("r_create_pos: zero-size group handled correctly", {
  pos <- r_create_pos(c(2L, 0L, 3L))
  # c(1, 3, 3, 6)
  expect_equal(pos, c(1L, 3L, 3L, 6L))
})

test_that("r_create_pos: all-zero sizes yields all-ones pos", {
  pos <- r_create_pos(c(0L, 0L, 0L))
  expect_equal(pos, c(1L, 1L, 1L, 1L))
})

test_that("r_get_pos: first group", {
  pos    <- c(1L, 3L, 3L, 6L)
  result <- r_get_pos(pos, 1L)
  expect_equal(result, c(1L, 2L))
})

test_that("r_get_pos: empty group (start > end)", {
  pos    <- c(1L, 3L, 3L, 6L)
  result <- r_get_pos(pos, 2L)
  # group 2 is empty: start=3, end=2
  expect_equal(result[1], 3L)
  expect_equal(result[2], 2L)
})

test_that("r_get_pos_size_single: empty group → 0", {
  pos <- c(1L, 3L, 3L, 6L)
  expect_equal(r_get_pos_size_single(pos, 2L), 0L)
  expect_equal(r_get_pos_size_single(pos, 1L), 2L)
  expect_equal(r_get_pos_size_single(pos, 3L), 3L)
})

test_that("r_get_pos_size_all: matches diff(pos)", {
  pos <- c(1L, 3L, 3L, 6L)
  expect_equal(r_get_pos_size_all(pos), c(2L, 0L, 3L))
})

test_that("r_get_pos_total_size: sum of all group sizes", {
  pos <- c(1L, 3L, 3L, 6L)
  expect_equal(r_get_pos_total_size(pos), 5L)
  # all-zero case
  expect_equal(r_get_pos_total_size(c(1L, 1L, 1L, 1L)), 0L)
})

test_that("r_get_int_sub_array: correct elements extracted", {
  full <- c(10L, 20L, 30L, 40L, 50L)
  pos  <- c(1L, 3L, 3L, 6L)
  # group 1: elements 1:2 = c(10, 20)
  expect_equal(r_get_int_sub_array(full, pos, 1L), c(10L, 20L))
  # group 3: elements 3:5 = c(30, 40, 50)
  expect_equal(r_get_int_sub_array(full, pos, 3L), c(30L, 40L, 50L))
})

test_that("r_get_int: correct single element", {
  full <- c(10L, 20L, 30L, 40L, 50L)
  pos  <- c(1L, 3L, 3L, 6L)
  # group 1, element 2 → full[pos[1] + 2 - 1] = full[2] = 20
  expect_equal(r_get_int(full, pos, 1L, 2L), 20L)
  # group 3, element 2 → full[pos[3] + 2 - 1] = full[4] = 40
  expect_equal(r_get_int(full, pos, 3L, 2L), 40L)
})

test_that("r_get_min_pos / r_get_max_pos", {
  full <- c(5L, 2L, 8L, 1L, 9L)
  pos  <- c(1L, 3L, 6L)
  # group 1: c(5, 2) → min=2, max=5
  expect_equal(r_get_min_pos(full, pos, 1L), 2L)
  expect_equal(r_get_max_pos(full, pos, 1L), 5L)
  # group 2: c(8, 1, 9) → min=1, max=9
  expect_equal(r_get_min_pos(full, pos, 2L), 1L)
  expect_equal(r_get_max_pos(full, pos, 2L), 9L)
})

test_that("r_create_enabled_pos: disabled group contributes 0", {
  result <- r_create_enabled_pos(c(3L, 5L, 2L), c(1L, 0L, 1L))
  # group 1 enabled (3), group 2 disabled (0), group 3 enabled (2) → c(1, 4, 4, 6)
  expect_equal(result, c(1L, 4L, 4L, 6L))
})

test_that("r_compute_n_enabled_groups: sums enabled group sizes only", {
  expect_equal(r_compute_n_enabled_groups(c(3L, 5L, 2L), c(1L, 0L, 1L)), 5L)
  expect_equal(r_compute_n_enabled_groups(c(3L, 5L), c(0L, 0L)), 0L)
  expect_equal(r_compute_n_enabled_groups(c(3L, 5L), c(1L, 1L)), 8L)
})

test_that("r_get_global_group_idx: level offset + group_id - 1", {
  level_pos <- c(1L, 4L, 9L)
  # level=2, group=3 → 4 + 3 - 1 = 6
  expect_equal(r_get_global_group_idx(level_pos, 2L, 3L), 6L)
  # level=1, group=1 → 1 + 1 - 1 = 1
  expect_equal(r_get_global_group_idx(level_pos, 1L, 1L), 1L)
})

test_that("r_validate_pos: non-decreasing input passes through unchanged", {
  pos <- c(1L, 3L, 3L, 6L)
  expect_equal(r_validate_pos(pos), pos)
})
