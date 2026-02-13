// ============================================================================
// Multistate Hazard Model Transformed Data
// ============================================================================

// --- Validate Final State Upper Bound ---
int ms_max_state = (enable_ms_02 || enable_ms_12) ? 2 : 1;
for (i in 1:n_patients) {
  if (ms_final_state[i] > ms_max_state) {
    fatal_error("Patient ", i, " has ms_final_state=", ms_final_state[i],
                " but max reachable state is ", ms_max_state,
                " given transition flags (enable_ms_02=", enable_ms_02,
                ", enable_ms_12=", enable_ms_12, ")");
  }
}

// --- Derived Flags for Time Scale ---
// Which GPs are needed for the 1→2 transition
int need_12_s_gp = enable_ms_12 && (ms_time_scale_12 == 1 || ms_time_scale_12 == 2);
int need_12_t_gp = enable_ms_12 && (ms_time_scale_12 == 0 || ms_time_scale_12 == 2);

// --- Enabled Group Counts for Baseline Hazard N-level Hierarchy ---
int n_enabled_groups_ms_baseline_01 = enable_ms_01 ? compute_n_enabled_groups(
  n_groups_per_level, enable_ms_level_baseline_hazard
) : 0;

int n_enabled_groups_ms_baseline_02 = enable_ms_02 ? compute_n_enabled_groups(
  n_groups_per_level, enable_ms_level_baseline_hazard
) : 0;

int n_enabled_groups_ms_baseline_12_s = need_12_s_gp ? compute_n_enabled_groups(
  n_groups_per_level, enable_ms_level_baseline_hazard
) : 0;

int n_enabled_groups_ms_baseline_12_t = need_12_t_gp ? compute_n_enabled_groups(
  n_groups_per_level, enable_ms_level_baseline_hazard
) : 0;

// --- Position Arrays for Baseline Hazard Level Hierarchy ---
array[n_levels + 1] int enabled_level_pos_ms_baseline = create_enabled_pos(
  n_groups_per_level, enable_ms_level_baseline_hazard
);

// --- Enabled Group Counts for Covariate Slopes ---
int n_enabled_groups_ms_slope = compute_n_enabled_groups(
  n_groups_per_level, enable_ms_level_cov
);

// Position array for enabled slope levels
array[n_levels + 1] int enabled_level_pos_ms_slope = create_enabled_pos(
  n_groups_per_level, enable_ms_level_cov
);

// --- Pre-computed Flat Indices for Patient Lookups ---
// Baseline hazard level indices
array[n_patients, n_levels] int patient_ms_baseline_flat_idx;
{
  for (i in 1:n_patients) {
    for (lv in 1:n_levels) {
      if (enable_ms_level_baseline_hazard[lv]) {
        patient_ms_baseline_flat_idx[i, lv] =
          get_global_group_idx(enabled_level_pos_ms_baseline, lv, patient_level_groups[i, lv]);
      } else {
        patient_ms_baseline_flat_idx[i, lv] = 1;
      }
    }
  }
}

// Covariate slope level indices
array[n_patients, n_levels] int patient_ms_slope_flat_idx;
{
  for (i in 1:n_patients) {
    for (lv in 1:n_levels) {
      if (enable_ms_level_cov[lv]) {
        patient_ms_slope_flat_idx[i, lv] =
          get_global_group_idx(enabled_level_pos_ms_slope, lv, patient_level_groups[i, lv]);
      } else {
        patient_ms_slope_flat_idx[i, lv] = 1;
      }
    }
  }
}

