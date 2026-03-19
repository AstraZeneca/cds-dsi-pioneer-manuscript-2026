# tests/testthat/helper-pos.R
#
# Pure R oracles for pos.stanfunctions.
# Position arrays encode ragged arrays: for group sizes [s1, s2, ..., sn],
# create_pos returns [1, 1+s1, 1+s1+s2, ..., 1+sum(si)] (1-based).

# create_pos: cumulative sum, 1-based, length = n+1
r_create_pos <- function(n_x) {
  c(1L, 1L + cumsum(as.integer(n_x)))
}

# get_pos: (start, end) indices for group i (1-based)
# Empty group: start = end + 1 (no elements)
r_get_pos <- function(pos, i) c(pos[i], pos[i + 1L] - 1L)

# get_pos_size: number of elements in group i
r_get_pos_size_single <- function(pos, i) pos[i + 1L] - pos[i]

# get_pos_size (all groups): diff of consecutive pos entries
r_get_pos_size_all <- function(pos) diff(as.integer(pos))

# get_pos_total_size: total number of elements across all groups
r_get_pos_total_size <- function(pos) pos[length(pos)] - 1L

# get_int_sub_array: extract elements for group n from full integer array
r_get_int_sub_array <- function(full, pos, n) {
  s <- r_get_pos(pos, n)
  full[s[1]:s[2]]
}

# get_int: element j within group n (1-based j)
r_get_int <- function(full, pos, n, j) full[pos[n] + j - 1L]

# get_min_pos / get_max_pos for group i
r_get_min_pos <- function(full, pos, i) {
  s <- r_get_pos(pos, i)
  min(full[s[1]:s[2]])
}
r_get_max_pos <- function(full, pos, i) {
  s <- r_get_pos(pos, i)
  max(full[s[1]:s[2]])
}

# create_enabled_pos: like create_pos but disabled groups contribute 0
r_create_enabled_pos <- function(n_groups, enabled) {
  n   <- length(n_groups)
  pos <- integer(n + 1L)
  pos[1L] <- 1L
  for (i in seq_len(n)) {
    pos[i + 1L] <- pos[i] + if (enabled[i] != 0L) n_groups[i] else 0L
  }
  pos
}

# compute_n_enabled_groups: total elements in enabled groups
r_compute_n_enabled_groups <- function(n_groups, enabled) {
  sum(n_groups[enabled != 0L])
}

# get_global_group_idx: level_pos[level] + group_id - 1 (1-based)
r_get_global_group_idx <- function(level_pos, level, group_id) {
  level_pos[level] + group_id - 1L
}

# validate_pos: returns pos unchanged for non-decreasing input; errors otherwise
r_validate_pos <- function(pos) {
  stopifnot(all(diff(pos) >= 0L))
  pos
}
