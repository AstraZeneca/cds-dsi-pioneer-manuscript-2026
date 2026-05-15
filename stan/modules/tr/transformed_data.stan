// tr/transformed_data.stan
// Compute enabled group counts, position arrays, and pre-computed indices for
// efficient vectorized operations in transformed_parameters.

int n_enabled_groups_tr_intercept;
array[n_levels + 1] int enabled_level_pos_tr_intercept;
array[n_patients, n_levels] int patient_tr_intercept_flat_idx;
int n_enabled_groups_tr_slope;
array[n_levels + 1] int enabled_level_pos_tr_slope;
array[n_patients, n_levels] int patient_tr_slope_flat_idx;
(n_enabled_groups_tr_intercept, enabled_level_pos_tr_intercept, patient_tr_intercept_flat_idx,
 n_enabled_groups_tr_slope, enabled_level_pos_tr_slope, patient_tr_slope_flat_idx) =
  compute_level_module_flags(n_patients, n_levels, n_forecast_patients,
    n_forecast_groups_per_level, patient_level_groups,
    enable_level_intercept_tr, enable_level_cov_tr);

// Split intercept groups into NCP (_raw_) and CP (_cp_) buckets by level mode.
// Groups at NONE levels are skipped; FE/RE/RE_GP go to raw, RE_CP goes to cp.
int n_raw_groups_tr_intercept;
int n_cp_groups_tr_intercept;
array[n_levels + 1] int raw_level_pos_tr_intercept;
array[n_levels + 1] int cp_level_pos_tr_intercept;
(n_raw_groups_tr_intercept, raw_level_pos_tr_intercept,
 n_cp_groups_tr_intercept,  cp_level_pos_tr_intercept) =
  split_cp_ncp_pos(n_levels, n_forecast_groups_per_level, enable_level_intercept_tr);

// Count RE-only levels for right-sizing the SD parameter array — RE_CP also needs
// a free SD hyperparameter (same as RE), so both count here.
int n_re_levels_tr_intercept = 0;
for (lv in 1:n_levels)
  if (enable_level_intercept_tr[lv] == LEVEL_MODE_RE ||
      enable_level_intercept_tr[lv] == LEVEL_MODE_RE_CP) n_re_levels_tr_intercept += 1;

// Split slope groups the same way, but routed by a combined mask: slope is only
// enabled when enable_level_cov_tr[lv] == 1 AND the level has an intercept.
// When routed, the slope parameterization follows the level's intercept mode.
array[n_levels] int tr_slope_mode;
for (lv in 1:n_levels) {
  tr_slope_mode[lv] = enable_level_cov_tr[lv] ? enable_level_intercept_tr[lv] : 0;
}
int n_raw_groups_tr_slope;
int n_cp_groups_tr_slope;
array[n_levels + 1] int raw_level_pos_tr_slope;
array[n_levels + 1] int cp_level_pos_tr_slope;
(n_raw_groups_tr_slope, raw_level_pos_tr_slope,
 n_cp_groups_tr_slope,  cp_level_pos_tr_slope) =
  split_cp_ncp_pos(n_levels, n_forecast_groups_per_level, tr_slope_mode);

// ===== SD sub-hierarchy sizing (issue #110) =====
// All derived from enable_sd_level_intercept_mode_tr (declared in flags.stan).

// Per-level active flag: sum over row L > 0.
array[n_levels] int has_sd_subhierarchy_tr_intercept;
int n_subhier_active_tr_intercept = 0;
for (lv in 1:n_levels) {
  int row_sum = 0;
  for (sub_lv in 1:n_levels) {
    row_sum += enable_sd_level_intercept_mode_tr[lv, sub_lv];
  }
  has_sd_subhierarchy_tr_intercept[lv] = (row_sum > 0) ? 1 : 0;
  n_subhier_active_tr_intercept += has_sd_subhierarchy_tr_intercept[lv];
}

// Position arrays via split_sd_cp_ncp_pos. Total counts size the raw/cp buckets.
int n_raw_groups_tr_log_sd_intercept;
array[n_levels, n_levels + 1] int raw_pos_tr_log_sd_intercept;
int n_cp_groups_tr_log_sd_intercept;
array[n_levels, n_levels + 1] int cp_pos_tr_log_sd_intercept;
(
  n_raw_groups_tr_log_sd_intercept,
  raw_pos_tr_log_sd_intercept,
  n_cp_groups_tr_log_sd_intercept,
  cp_pos_tr_log_sd_intercept
) = split_sd_cp_ncp_pos(
  n_levels,
  n_forecast_groups_per_level,
  enable_sd_level_intercept_mode_tr
);

// Hyperscale count: one per (L, ℓ) with mode ∈ {RE, RE_CP} (not FE — FE uses data value).
int n_sd_hyperscales_tr_intercept = 0;
for (lv in 1:n_levels) {
  for (sub_lv in 1:(lv - 1)) {
    int m = enable_sd_level_intercept_mode_tr[lv, sub_lv];
    if (m == LEVEL_MODE_RE || m == LEVEL_MODE_RE_CP) {
      n_sd_hyperscales_tr_intercept += 1;
    }
  }
}

// Flat hyperscale index: hyperscale_idx_tr_intercept[L, ℓ] = 1-based position in
// tr_sd_hyperscale_level_intercept_raw, or 0 if not a RE/RE_CP entry.
array[n_levels, n_levels] int hyperscale_idx_tr_intercept;
{
  int idx = 0;
  for (lv in 1:n_levels) {
    for (sub_lv in 1:n_levels) {
      int m = enable_sd_level_intercept_mode_tr[lv, sub_lv];
      if (m == LEVEL_MODE_RE || m == LEVEL_MODE_RE_CP) {
        idx += 1;
        hyperscale_idx_tr_intercept[lv, sub_lv] = idx;
      } else {
        hyperscale_idx_tr_intercept[lv, sub_lv] = 0;
      }
    }
  }
}

// Pop-vector index: pop_idx_tr_intercept[L] = 1-based position in
// tr_log_sd_level_intercept_pop, or 0 if L has no active sub-hierarchy.
array[n_levels] int pop_idx_tr_intercept;
{
  int idx = 0;
  for (lv in 1:n_levels) {
    if (has_sd_subhierarchy_tr_intercept[lv]) {
      idx += 1;
      pop_idx_tr_intercept[lv] = idx;
    } else {
      pop_idx_tr_intercept[lv] = 0;
    }
  }
}

// Subgroup flat index: for each location-level-lv group g, which sub-level-ℓ group
// does it map to? Derived from patient_level_groups by picking a representative
// patient per lv-group (well-defined because the implementation is nested).
// Flat layout: subgroup_idx_tr_intercept_flat[lv][sub_lv][g_lv] but stored as a
// jagged array via enabled_level_pos-style indexing.
//
// Simplification for Phase 1: store only the lookup we need — for each location
// level lv with active sub-hierarchy and each active sub-level ℓ, a vector of
// length n_forecast_groups_per_level[lv] giving the ℓ-group id for each lv-group.
array[n_levels, n_levels, max(n_forecast_groups_per_level)] int subgroup_idx_tr_intercept;
for (lv in 1:n_levels) {
  if (!has_sd_subhierarchy_tr_intercept[lv]) continue;
  int ngroups_lv = n_forecast_groups_per_level[lv];
  for (sub_lv in 1:(lv - 1)) {
    if (enable_sd_level_intercept_mode_tr[lv, sub_lv] == LEVEL_MODE_NONE) continue;
    // Build a representative-patient map from lv-group → ℓ-group.
    // patient_level_groups[i, lv] = group id at level lv for patient i.
    // For each lv-group g, find the first patient with that group id and read its ℓ-group.
    array[ngroups_lv] int seen = zeros_int_array(ngroups_lv);
    for (i in 1:n_patients) {
      int g_lv = patient_level_groups[i, lv];
      if (g_lv >= 1 && g_lv <= ngroups_lv && seen[g_lv] == 0) {
        subgroup_idx_tr_intercept[lv, sub_lv, g_lv] = patient_level_groups[i, sub_lv];
        seen[g_lv] = 1;
      }
    }
  }
}
