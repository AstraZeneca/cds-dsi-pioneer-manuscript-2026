// Test harness for compute_level_module_flags in stan/hierarchy.stanfunctions.
// Runs all flag permutations in a single fixed_param pass.
functions {
  #include "util.stanfunctions"
  #include "pos.stanfunctions"
  #include "hierarchy.stanfunctions"
}

data {
  int<lower=1> n_combos;
  int<lower=1> n_patients;
  int<lower=1> n_levels_max;
  int<lower=1> n_forecast_patients;
  // Per-combo flag arrays (padded to n_levels_max):
  array[n_combos, n_levels_max] int enable_intercept_combos;
  array[n_combos, n_levels_max] int enable_slope_combos;
  array[n_combos] int n_levels_combos;
  // n_forecast_groups per level, padded with 1s (padded levels unused)
  array[n_combos, n_levels_max] int n_forecast_groups_combos;
  // Patient group memberships — shared across all combos (padded to n_levels_max)
  array[n_patients, n_levels_max] int patient_level_groups;
}

generated quantities {
  array[n_combos] int out_n_groups_intercept;
  array[n_combos] int out_n_groups_slope;
  // Position arrays padded to n_levels_max+1 per combo
  array[n_combos, n_levels_max + 1] int out_pos_intercept;
  array[n_combos, n_levels_max + 1] int out_pos_slope;
  // Flat index matrices: n_combos x n_patients x n_levels_max
  array[n_combos, n_patients, n_levels_max] int out_flat_idx_intercept;
  array[n_combos, n_patients, n_levels_max] int out_flat_idx_slope;

  for (c in 1:n_combos) {
    int nl = n_levels_combos[c];

    // Slice flag vectors to actual n_levels for this combo
    array[nl] int ei;
    array[nl] int es;
    array[nl] int fg;
    for (k in 1:nl) {
      ei[k] = enable_intercept_combos[c, k];
      es[k] = enable_slope_combos[c, k];
      fg[k] = n_forecast_groups_combos[c, k];
    }

    // Slice patient_level_groups columns to nl
    array[n_patients, nl] int plg;
    for (i in 1:n_patients) {
      for (k in 1:nl) {
        plg[i, k] = patient_level_groups[i, k];
      }
    }

    int n_int;
    array[nl + 1] int pos_int;
    array[n_patients, nl] int flat_int;
    int n_slp;
    array[nl + 1] int pos_slp;
    array[n_patients, nl] int flat_slp;

    (n_int, pos_int, flat_int, n_slp, pos_slp, flat_slp) =
      compute_level_module_flags(n_patients, nl, n_forecast_patients, fg, plg, ei, es);

    out_n_groups_intercept[c] = n_int;
    out_n_groups_slope[c]     = n_slp;

    // Copy pos arrays, pad remainder with 0
    for (k in 1:(nl + 1)) out_pos_intercept[c, k] = pos_int[k];
    // Padding loops: when nl == n_levels_max (typical single-n_levels test),
    // these ranges are zero-length and the loops are skipped.
    // They are here to support future multi-n_levels batch tests.
    for (k in (nl + 2):(n_levels_max + 1)) out_pos_intercept[c, k] = 0;
    for (k in 1:(nl + 1)) out_pos_slope[c, k] = pos_slp[k];
    for (k in (nl + 2):(n_levels_max + 1)) out_pos_slope[c, k] = 0;

    // Copy flat_idx, pad remainder with 0
    for (i in 1:n_patients) {
      for (k in 1:nl) {
        out_flat_idx_intercept[c, i, k] = flat_int[i, k];
        out_flat_idx_slope[c, i, k]     = flat_slp[i, k];
      }
      for (k in (nl + 1):n_levels_max) {
        out_flat_idx_intercept[c, i, k] = 0;
        out_flat_idx_slope[c, i, k]     = 0;
      }
    }
  }
}
