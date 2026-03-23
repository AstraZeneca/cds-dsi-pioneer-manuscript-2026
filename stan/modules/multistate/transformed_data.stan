// ============================================================================
// Multistate Hazard Model Transformed Data
// ============================================================================

// --- Validate Final State Upper Bound ---
int ms_max_state = (enable_ms_02 || enable_ms_12) ? 2 : 1;
if (enable_ms_03) ms_max_state = 3;
for (i in 1:n_patients) {
  if (ms_final_state[i] > ms_max_state) {
    fatal_error("Patient ", i, " has ms_final_state=", ms_final_state[i],
                " but max reachable state is ", ms_max_state,
                " given transition flags (enable_ms_02=", enable_ms_02,
                ", enable_ms_12=", enable_ms_12,
                ", enable_ms_03=", enable_ms_03, ")");
  }
}

if (!enable_ms_03 && enable_ms_32) {
  fatal_error("enable_ms_32=1 requires enable_ms_03=1");
}

// --- Derived Flags for Time Scale ---
// Which GPs are needed for the 1→2 transition
int need_12_s_gp = enable_ms_12 && (ms_time_scale_12 == 1 || ms_time_scale_12 == 2);
int need_12_t_gp = enable_ms_12 && (ms_time_scale_12 == 0 || ms_time_scale_12 == 2);

// In extended mode (2), both GPs are active and additive. To avoid
// non-identifiability of two intercepts, the clock-forward GP is zero-mean
// and only the sojourn GP carries the intercept.
int ms_12_t_has_intercept = need_12_t_gp && !need_12_s_gp;  // Only in pure Markov mode

// --- Enabled Group Counts for Baseline Hazard N-level Hierarchy ---
int n_enabled_groups_ms_baseline_01 = enable_ms_01 ? compute_n_enabled_groups(
  n_hmc_groups_per_level, enable_ms_level_baseline_hazard
) : 0;

int n_enabled_groups_ms_baseline_02 = enable_ms_02 ? compute_n_enabled_groups(
  n_hmc_groups_per_level, enable_ms_level_baseline_hazard
) : 0;

int n_enabled_groups_ms_baseline_12_s = need_12_s_gp ? compute_n_enabled_groups(
  n_hmc_groups_per_level, enable_ms_level_baseline_hazard
) : 0;

int n_enabled_groups_ms_baseline_12_t = need_12_t_gp ? compute_n_enabled_groups(
  n_hmc_groups_per_level, enable_ms_level_baseline_hazard
) : 0;

int n_enabled_groups_ms_baseline_03 = enable_ms_03 ? compute_n_enabled_groups(
  n_hmc_groups_per_level, enable_ms_level_baseline_hazard
) : 0;

int n_enabled_groups_ms_baseline_32 = enable_ms_32 ? compute_n_enabled_groups(
  n_hmc_groups_per_level, enable_ms_level_baseline_hazard
) : 0;

// --- GP-Only Boolean Mask and Group Counts ---
// For intercept-only mode (mode==1), we don't need eta vectors.
// GP-only arrays track which levels use full GP (mode==2).
array[n_levels] int ms_level_baseline_is_gp;
for (lv in 1:n_levels) {
  ms_level_baseline_is_gp[lv] = (enable_ms_level_baseline_hazard[lv] == 2) ? 1 : 0;
}

int n_gp_groups_ms_baseline = compute_n_enabled_groups(
  n_hmc_groups_per_level, ms_level_baseline_is_gp
);

int n_gp_groups_ms_baseline_01 = enable_ms_01 ? n_gp_groups_ms_baseline : 0;
int n_gp_groups_ms_baseline_02 = enable_ms_02 ? n_gp_groups_ms_baseline : 0;
int n_gp_groups_ms_baseline_12_s = need_12_s_gp ? n_gp_groups_ms_baseline : 0;
int n_gp_groups_ms_baseline_12_t = need_12_t_gp ? n_gp_groups_ms_baseline : 0;
int n_gp_groups_ms_baseline_03 = enable_ms_03 ? n_gp_groups_ms_baseline : 0;
int n_gp_groups_ms_baseline_32 = enable_ms_32 ? n_gp_groups_ms_baseline : 0;

// GP-only position array (for indexing into eta matrices)
array[n_levels + 1] int gp_level_pos_ms_baseline = create_enabled_pos(
  n_hmc_groups_per_level, ms_level_baseline_is_gp
);

// --- Position Arrays for Baseline Hazard Level Hierarchy ---
// (includes both intercept-only and GP levels — any truthy flag)
array[n_levels + 1] int enabled_level_pos_ms_baseline = create_enabled_pos(
  n_hmc_groups_per_level, enable_ms_level_baseline_hazard
);

// --- Enabled Group Counts for Covariate Slopes ---
int n_enabled_groups_ms_slope = compute_n_enabled_groups(
  n_hmc_groups_per_level, enable_ms_level_cov
);

// Position array for enabled slope levels
array[n_levels + 1] int enabled_level_pos_ms_slope = create_enabled_pos(
  n_hmc_groups_per_level, enable_ms_level_cov
);

// --- GP Coarse Knot Grids ---
// Knot counts per time domain
int n_ms_gp_cal_knots      = (max_all_t          + ms_gp_grid_step - 1) %/% ms_gp_grid_step;
int n_ms_gp_sojourn_knots  = (ms_max_sojourn_t   + ms_gp_grid_step - 1) %/% ms_gp_grid_step;
int n_ms_gp_sojourn_32_knots = (ms_max_sojourn_t_32 + ms_gp_grid_step - 1) %/% ms_gp_grid_step;

// Knot positions (real-valued, at j * ms_gp_grid_step for j=1..n_knots)
array[n_ms_gp_cal_knots] real ms_gp_cal_t;
for (j in 1:n_ms_gp_cal_knots)
  ms_gp_cal_t[j] = j * ms_gp_grid_step * 1.0;

array[n_ms_gp_sojourn_knots] real ms_gp_sojourn_t;
for (j in 1:n_ms_gp_sojourn_knots)
  ms_gp_sojourn_t[j] = j * ms_gp_grid_step * 1.0;

array[n_ms_gp_sojourn_32_knots] real ms_gp_sojourn_32_t;
for (j in 1:n_ms_gp_sojourn_32_knots)
  ms_gp_sojourn_32_t[j] = j * ms_gp_grid_step * 1.0;

// Week-to-nearest-knot mappings
array[max_all_t] int knot_of_cal;
for (t in 1:max_all_t) {
  int idx = (t + ms_gp_grid_step %/% 2) %/% ms_gp_grid_step;
  if (idx < 1) idx = 1;
  if (idx > n_ms_gp_cal_knots) idx = n_ms_gp_cal_knots;
  knot_of_cal[t] = idx;
}

array[ms_max_sojourn_t] int knot_of_sojourn;
for (t in 1:ms_max_sojourn_t) {
  int idx = (t + ms_gp_grid_step %/% 2) %/% ms_gp_grid_step;
  if (idx < 1) idx = 1;
  if (idx > n_ms_gp_sojourn_knots) idx = n_ms_gp_sojourn_knots;
  knot_of_sojourn[t] = idx;
}

array[ms_max_sojourn_t_32] int knot_of_sojourn_32;
for (t in 1:ms_max_sojourn_t_32) {
  int idx = (t + ms_gp_grid_step %/% 2) %/% ms_gp_grid_step;
  if (idx < 1) idx = 1;
  if (idx > n_ms_gp_sojourn_32_knots) idx = n_ms_gp_sojourn_32_knots;
  knot_of_sojourn_32[t] = idx;
}

// --- Pre-computed Flat Indices for Patient Lookups ---
// Baseline hazard level indices
array[n_patients, n_levels] int patient_ms_baseline_flat_idx;
{
  for (i in 1:n_patients) {
    for (lv in 1:n_levels) {
      if (enable_ms_level_baseline_hazard[lv]) {
        patient_ms_baseline_flat_idx[i, lv] =
          (lv == n_levels && patient_level_groups[i, lv] > n_hmc_patients) ? 1
          : get_global_group_idx(enabled_level_pos_ms_baseline, lv, patient_level_groups[i, lv]);
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
          (lv == n_levels && patient_level_groups[i, lv] > n_hmc_patients) ? 1
          : get_global_group_idx(enabled_level_pos_ms_slope, lv, patient_level_groups[i, lv]);
      } else {
        patient_ms_slope_flat_idx[i, lv] = 1;
      }
    }
  }
}

// --- IC Gap for 0→1 Marginalization ---
// ms_ic_gap_01[i] = number of candidate weeks in (T_c, T_d] for patient i.
// = interval_censored[i] + 1 for stochastic-progression patients (gap > 0).
// = 0 otherwise → IC code is bypassed, reducing to the no-IC likelihood.
array[n_patients] int ms_ic_gap_01;
for (i in 1:n_patients) {
  if (ms_censored_01[i] || ms_prog_deterministic[i]) {
    ms_ic_gap_01[i] = 0;
  } else {
    ms_ic_gap_01[i] = interval_censored[i] + 1;
  }
}

