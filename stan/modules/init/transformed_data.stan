// init/transformed_data.stan
// Compute enabled group counts, position arrays, and pre-computed indices for
// efficient vectorized operations in transformed_parameters.

int n_enabled_groups_init_intercept;
array[n_levels + 1] int enabled_level_pos_init_intercept;
array[n_patients, n_levels] int patient_init_intercept_flat_idx;
int n_enabled_groups_init_slope;
array[n_levels + 1] int enabled_level_pos_init_slope;
array[n_patients, n_levels] int patient_init_slope_flat_idx;
(n_enabled_groups_init_intercept, enabled_level_pos_init_intercept, patient_init_intercept_flat_idx,
 n_enabled_groups_init_slope, enabled_level_pos_init_slope, patient_init_slope_flat_idx) =
  compute_level_module_flags(n_patients, n_levels, n_forecast_patients,
    n_forecast_groups_per_level, patient_level_groups,
    enable_level_intercept_init, enable_level_cov_init);

// Count RE and RE_CP levels for right-sizing the SD parameter array
int n_re_levels_init_intercept = 0;
for (lv in 1:n_levels)
  if (enable_level_intercept_init[lv] == LEVEL_MODE_RE ||
      enable_level_intercept_init[lv] == LEVEL_MODE_RE_CP) n_re_levels_init_intercept += 1;

// Split intercept groups into NCP (_raw_) and CP (_cp_) buckets by level mode.
int n_raw_groups_init_intercept;
int n_cp_groups_init_intercept;
array[n_levels + 1] int raw_level_pos_init_intercept;
array[n_levels + 1] int cp_level_pos_init_intercept;
(n_raw_groups_init_intercept, raw_level_pos_init_intercept,
 n_cp_groups_init_intercept,  cp_level_pos_init_intercept) =
  split_cp_ncp_pos(n_levels, n_forecast_groups_per_level, enable_level_intercept_init);

// Split slope groups the same way: slope is only enabled when
// enable_level_cov_init[lv] == 1. When routed, slope parameterization follows the
// level's intercept mode.
array[n_levels] int init_slope_mode;
for (lv in 1:n_levels) {
  init_slope_mode[lv] = enable_level_cov_init[lv] ? enable_level_intercept_init[lv] : 0;
}
int n_raw_groups_init_slope;
int n_cp_groups_init_slope;
array[n_levels + 1] int raw_level_pos_init_slope;
array[n_levels + 1] int cp_level_pos_init_slope;
(n_raw_groups_init_slope, raw_level_pos_init_slope,
 n_cp_groups_init_slope,  cp_level_pos_init_slope) =
  split_cp_ncp_pos(n_levels, n_forecast_groups_per_level, init_slope_mode);
