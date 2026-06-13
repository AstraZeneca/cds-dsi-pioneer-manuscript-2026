// gr_decay/transformed_data.stan
// Flat indices + bucket positions for the level-indexed log(kappa) hierarchy.
// Mirrors frac/transformed_data.stan. Computed unconditionally (cheap; all counts
// collapse to 0 when every level mode is NONE, so the pop-only path is bit-exact).

int n_enabled_groups_gr_decay_intercept;
array[n_levels + 1] int enabled_level_pos_gr_decay_intercept;
array[n_patients, n_levels] int patient_gr_decay_intercept_flat_idx;
int n_enabled_groups_gr_decay_slope;
array[n_levels + 1] int enabled_level_pos_gr_decay_slope;
array[n_patients, n_levels] int patient_gr_decay_slope_flat_idx;
(n_enabled_groups_gr_decay_intercept, enabled_level_pos_gr_decay_intercept, patient_gr_decay_intercept_flat_idx,
 n_enabled_groups_gr_decay_slope, enabled_level_pos_gr_decay_slope, patient_gr_decay_slope_flat_idx) =
  compute_level_module_flags(n_patients, n_levels, n_forecast_patients,
    n_forecast_groups_per_level, patient_level_groups,
    enable_level_intercept_gr_decay, enable_level_cov_gr_decay);

int n_re_levels_gr_decay_intercept = 0;
for (lv in 1:n_levels)
  if (enable_level_intercept_gr_decay[lv] == LEVEL_MODE_RE ||
      enable_level_intercept_gr_decay[lv] == LEVEL_MODE_RE_CP) n_re_levels_gr_decay_intercept += 1;

int n_raw_groups_gr_decay_intercept;
int n_cp_groups_gr_decay_intercept;
array[n_levels + 1] int raw_level_pos_gr_decay_intercept;
array[n_levels + 1] int cp_level_pos_gr_decay_intercept;
(n_raw_groups_gr_decay_intercept, raw_level_pos_gr_decay_intercept,
 n_cp_groups_gr_decay_intercept,  cp_level_pos_gr_decay_intercept) =
  split_cp_ncp_pos(n_levels, n_forecast_groups_per_level, enable_level_intercept_gr_decay);

array[n_levels] int gr_decay_slope_mode;
for (lv in 1:n_levels) {
  gr_decay_slope_mode[lv] = enable_level_cov_gr_decay[lv] ? enable_level_intercept_gr_decay[lv] : 0;
}
int n_raw_groups_gr_decay_slope;
int n_cp_groups_gr_decay_slope;
array[n_levels + 1] int raw_level_pos_gr_decay_slope;
array[n_levels + 1] int cp_level_pos_gr_decay_slope;
(n_raw_groups_gr_decay_slope, raw_level_pos_gr_decay_slope,
 n_cp_groups_gr_decay_slope,  cp_level_pos_gr_decay_slope) =
  split_cp_ncp_pos(n_levels, n_forecast_groups_per_level, gr_decay_slope_mode);
