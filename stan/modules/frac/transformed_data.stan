// frac/transformed_data.stan
// Compute enabled group counts, position arrays, and pre-computed indices for
// efficient vectorized operations in transformed_parameters.

int n_enabled_groups_frac_intercept;
array[n_levels + 1] int enabled_level_pos_frac_intercept;
array[n_patients, n_levels] int patient_frac_intercept_flat_idx;
int n_enabled_groups_frac_slope;
array[n_levels + 1] int enabled_level_pos_frac_slope;
array[n_patients, n_levels] int patient_frac_slope_flat_idx;
(n_enabled_groups_frac_intercept, enabled_level_pos_frac_intercept, patient_frac_intercept_flat_idx,
 n_enabled_groups_frac_slope, enabled_level_pos_frac_slope, patient_frac_slope_flat_idx) =
  compute_level_module_flags(n_patients, n_levels, n_forecast_patients,
    n_forecast_groups_per_level, patient_level_groups,
    enable_level_intercept_frac, enable_level_cov_frac);

// Count RE and RE_CP levels for right-sizing the SD parameter array
int n_re_levels_frac_intercept = 0;
for (lv in 1:n_levels)
  if (enable_level_intercept_frac[lv] == LEVEL_MODE_RE ||
      enable_level_intercept_frac[lv] == LEVEL_MODE_RE_CP) n_re_levels_frac_intercept += 1;

// Intercept bucket counts + positions (raw NCP vs CP)
int n_raw_groups_frac_intercept;
int n_cp_groups_frac_intercept;
array[n_levels + 1] int raw_level_pos_frac_intercept;
array[n_levels + 1] int cp_level_pos_frac_intercept;
(n_raw_groups_frac_intercept, raw_level_pos_frac_intercept,
 n_cp_groups_frac_intercept,  cp_level_pos_frac_intercept) =
  split_cp_ncp_pos(n_levels, n_forecast_groups_per_level, enable_level_intercept_frac);

// Slope bucket counts + positions — slope mode mirrors intercept mode when slope enabled
array[n_levels] int frac_slope_mode;
for (lv in 1:n_levels) {
  frac_slope_mode[lv] = enable_level_cov_frac[lv] ? enable_level_intercept_frac[lv] : 0;
}
int n_raw_groups_frac_slope;
int n_cp_groups_frac_slope;
array[n_levels + 1] int raw_level_pos_frac_slope;
array[n_levels + 1] int cp_level_pos_frac_slope;
(n_raw_groups_frac_slope, raw_level_pos_frac_slope,
 n_cp_groups_frac_slope,  cp_level_pos_frac_slope) =
  split_cp_ncp_pos(n_levels, n_forecast_groups_per_level, frac_slope_mode);
