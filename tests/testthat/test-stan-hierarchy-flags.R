library(testthat)
library(here)
library(posterior)
library(stringr)
library(dplyr)
library(purrr)

# Build the Stan data list for one (n_levels, n_patients, n_forecast_patients)
# configuration. All flag permutations for enable_intercept and enable_slope
# are expanded into a batch of combos.
make_hierarchy_flag_combos <- function(
    n_levels,
    n_patients,
    n_forecast_patients,
    n_groups_per_level
) {
  # All 2^n_levels combos for intercept x 2^n_levels for slope
  ei_list <- replicate(n_levels, 0:1, simplify = FALSE) |>
    setNames(str_c("ei_", seq_len(n_levels)))
  es_list <- replicate(n_levels, 0:1, simplify = FALSE) |>
    setNames(str_c("es_", seq_len(n_levels)))
  flag_grid <- expand.grid(c(ei_list, es_list))

  n_combos <- nrow(flag_grid)
  n_levels_max <- n_levels  # only one n_levels in this batch

  # Build padded arrays (n_combos x n_levels_max)
  ei_mat <- matrix(0L, n_combos, n_levels_max)
  es_mat <- matrix(0L, n_combos, n_levels_max)
  for (k in seq_len(n_levels)) {
    ei_mat[, k] <- as.integer(flag_grid[[str_c("ei_", k)]])
    es_mat[, k] <- as.integer(flag_grid[[str_c("es_", k)]])
  }

  # n_forecast_groups: use n_groups_per_level, substitute n_forecast_patients at last level
  fg <- n_groups_per_level
  fg[n_levels] <- n_forecast_patients

  fg_mat <- matrix(rep(fg, each = n_combos), nrow = n_combos)

  # Patient group memberships: assign each patient to group 1 at every level,
  # except at the patient level (n_levels) use patient identity 1..n_patients.
  # Make some patients "background" by assigning gid > n_forecast_patients at the
  # patient level when n_patients > n_forecast_patients.
  plg <- matrix(1L, n_patients, n_levels_max)
  for (i in seq_len(n_patients)) {
    plg[i, n_levels] <- i  # patient identity at last level
  }

  list(
    n_combos               = n_combos,
    n_patients             = n_patients,
    n_levels_max           = n_levels_max,
    n_forecast_patients    = n_forecast_patients,
    enable_intercept_combos = ei_mat,
    enable_slope_combos    = es_mat,
    n_levels_combos        = rep(n_levels, n_combos),
    n_forecast_groups_combos = fg_mat,
    patient_level_groups   = plg,
    # Keep R-side oracle data for checking:
    .flag_grid     = flag_grid,
    .fg            = fg,
    .n_levels      = n_levels,
    .n_combos      = n_combos
  )
}

# R oracle: compute flat_idx for one patient, one level, one flag vector
r_flat_idx_for_patient <- function(i, lv, n_levels, n_forecast_patients,
                                    fg, enable, plg) {
  if (enable[lv] == 0L) return(1L)
  gid <- plg[i, lv]
  pos <- r_create_enabled_pos(fg, enable)
  if (lv == n_levels && gid > n_forecast_patients) return(1L)
  r_get_global_group_idx(pos, lv, gid)
}

run_hierarchy_flag_test <- function(n_levels, n_patients, n_forecast_patients,
                                     n_groups_per_level) {
  dat <- make_hierarchy_flag_combos(
    n_levels, n_patients, n_forecast_patients, n_groups_per_level
  )

  stan_data <- dat[!startsWith(names(dat), ".")]
  fit <- test_stan_function(
    here("tests", "testthat", "stan", "test_hierarchy_flags_all.stan"),
    data = stan_data
  )
  d <- as_draws_df(fit$draws())

  flag_grid    <- dat$.flag_grid
  fg           <- dat$.fg
  n_lvs        <- dat$.n_levels
  n_combos     <- dat$.n_combos
  plg          <- dat$patient_level_groups

  for (c in seq_len(n_combos)) {
    ei <- as.integer(flag_grid[c, str_c("ei_", seq_len(n_lvs))])
    es <- as.integer(flag_grid[c, str_c("es_", seq_len(n_lvs))])

    # --- n_enabled_groups ---
    expect_equal(
      get_stan_val(d, "out_n_groups_intercept", c),
      r_compute_n_enabled_groups(fg, ei),
      label = sprintf("combo %d n_groups_intercept (ei=%s)", c, paste(ei, collapse = ""))
    )
    expect_equal(
      get_stan_val(d, "out_n_groups_slope", c),
      r_compute_n_enabled_groups(fg, es),
      label = sprintf("combo %d n_groups_slope (es=%s)", c, paste(es, collapse = ""))
    )

    # --- pos arrays (1..n_levels+1) ---
    expected_pos_int <- r_create_enabled_pos(fg, ei)
    expected_pos_slp <- r_create_enabled_pos(fg, es)
    for (k in seq_len(n_lvs + 1L)) {
      expect_equal(
        get_stan_val(d, "out_pos_intercept", c, k),
        expected_pos_int[k],
        label = sprintf("combo %d pos_intercept[%d] (ei=%s)", c, k, paste(ei, collapse = ""))
      )
      expect_equal(
        get_stan_val(d, "out_pos_slope", c, k),
        expected_pos_slp[k],
        label = sprintf("combo %d pos_slope[%d] (es=%s)", c, k, paste(es, collapse = ""))
      )
    }

    # --- flat_idx matrices ---
    for (i in seq_len(n_patients)) {
      for (lv in seq_len(n_lvs)) {
        expect_equal(
          get_stan_val(d, "out_flat_idx_intercept", c, i, lv),
          r_flat_idx_for_patient(i, lv, n_lvs, n_forecast_patients, fg, ei, plg),
          label = sprintf("combo %d flat_idx_int[%d,%d] (ei=%s)", c, i, lv, paste(ei, collapse = ""))
        )
        expect_equal(
          get_stan_val(d, "out_flat_idx_slope", c, i, lv),
          r_flat_idx_for_patient(i, lv, n_lvs, n_forecast_patients, fg, es, plg),
          label = sprintf("combo %d flat_idx_slp[%d,%d] (es=%s)", c, i, lv, paste(es, collapse = ""))
        )
      }
    }
  }
}

test_that(
  "compute_level_module_flags: correct results for 2-level hierarchy, all flag permutations",
  {
    # 2 levels: trial (2 groups) and patient (identity).
    # n_patients=5, n_forecast_patients=3 (patients 4,5 are background).
    run_hierarchy_flag_test(
      n_levels           = 2L,
      n_patients         = 5L,
      n_forecast_patients = 3L,
      n_groups_per_level  = c(2L, 5L)  # 2 trials, 5 patients total
    )
  }
)

test_that(
  "compute_level_module_flags: correct results for 3-level hierarchy, all flag permutations",
  {
    # 3 levels: population (1 group), trial (3 groups), patient (identity).
    # n_patients=4, n_forecast_patients=4 (no background patients).
    run_hierarchy_flag_test(
      n_levels           = 3L,
      n_patients         = 4L,
      n_forecast_patients = 4L,
      n_groups_per_level  = c(1L, 3L, 4L)
    )
  }
)
