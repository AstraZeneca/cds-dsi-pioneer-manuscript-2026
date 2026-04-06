// ============================================================================
// Multistate Hazard Model Transformed Data
// ============================================================================

// --- Validate Final State Upper Bound ---
int ms_max_state = (enable_ms_02 || enable_ms_12) ? 2 : 1;
if (enable_ms_03) ms_max_state = 3;
for (i in 1:n_patients) {
  // State 3 is always tolerated in data: when enable_ms_03=0 the likelihood
  // treats these patients as right-censored at patient_max_t (time_03).
  if (ms_final_state[i] != 3 && ms_final_state[i] > ms_max_state) {
    fatal_error("Patient ", i, " has ms_final_state=", ms_final_state[i],
                " but max reachable state is ", ms_max_state,
                " given transition flags (enable_ms_02=", enable_ms_02,
                ", enable_ms_12=", enable_ms_12,
                ", enable_ms_03=", enable_ms_03, ")");
  }
}

// --- Derived Flags for Time Scale (B1) ---
// Which GPs are needed for the 1→2 transition
// strict=1: fatal_error on invalid input; is_valid sentinel discarded
int ms_flag_b1_valid;
int need_12_s_gp;
int need_12_t_gp;
int ms_12_t_has_intercept;
(ms_flag_b1_valid, need_12_s_gp, need_12_t_gp, ms_12_t_has_intercept) =
  compute_ms_time_scale_flags(enable_ms_12, ms_time_scale_12, 1);

// --- Level Baseline Hazard Flags (B2) ---
// GP-Only Boolean Mask, group counts, and position arrays
// strict=1: fatal_error on invalid input; is_valid sentinel discarded
int ms_flag_b2_valid;
int any_re_level;
array[n_levels] int ms_level_baseline_is_gp;
int n_gp_groups_ms_baseline;
array[n_levels + 1] int gp_level_pos_ms_baseline;
array[n_levels + 1] int enabled_level_pos_ms_baseline;
(ms_flag_b2_valid, any_re_level, ms_level_baseline_is_gp, n_gp_groups_ms_baseline,
 gp_level_pos_ms_baseline, enabled_level_pos_ms_baseline) =
  compute_ms_level_baseline_flags(
    n_levels, n_forecast_groups_per_level, enable_ms_level_baseline_hazard, 1
  );

// --- Shared enabled group count (needed by B3) ---
int n_enabled_groups_ms_baseline = compute_n_enabled_groups(
  n_forecast_groups_per_level, enable_ms_level_baseline_hazard
);

// --- Per-Transition Group Counts (B3) ---
// Validates enable_ms_32 requires enable_ms_03 (strict=1: fatal_error on violation)
// is_valid sentinel discarded
int ms_flag_b3_valid;
int n_enabled_groups_ms_baseline_01; int n_gp_groups_ms_baseline_01;
int n_enabled_groups_ms_baseline_02; int n_gp_groups_ms_baseline_02;
int n_enabled_groups_ms_baseline_12_s; int n_gp_groups_ms_baseline_12_s;
int n_enabled_groups_ms_baseline_12_t; int n_gp_groups_ms_baseline_12_t;
int n_enabled_groups_ms_baseline_03; int n_gp_groups_ms_baseline_03;
int n_enabled_groups_ms_baseline_32; int n_gp_groups_ms_baseline_32;
(ms_flag_b3_valid,
 n_enabled_groups_ms_baseline_01, n_gp_groups_ms_baseline_01,
 n_enabled_groups_ms_baseline_02, n_gp_groups_ms_baseline_02,
 n_enabled_groups_ms_baseline_12_s, n_gp_groups_ms_baseline_12_s,
 n_enabled_groups_ms_baseline_12_t, n_gp_groups_ms_baseline_12_t,
 n_enabled_groups_ms_baseline_03, n_gp_groups_ms_baseline_03,
 n_enabled_groups_ms_baseline_32, n_gp_groups_ms_baseline_32) =
  compute_ms_transition_group_counts(
    enable_ms_01, enable_ms_02,
    need_12_s_gp, need_12_t_gp,
    enable_ms_03, enable_ms_32,
    n_enabled_groups_ms_baseline, n_gp_groups_ms_baseline, 1
  );

// --- Enabled Group Counts for Covariate Slopes ---
int n_enabled_groups_ms_slope = compute_n_enabled_groups(
  n_forecast_groups_per_level, enable_ms_level_cov
);

// Position array for enabled slope levels
array[n_levels + 1] int enabled_level_pos_ms_slope = create_enabled_pos(
  n_forecast_groups_per_level, enable_ms_level_cov
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
          (lv == n_levels && patient_level_groups[i, lv] > n_forecast_patients) ? 1
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
          (lv == n_levels && patient_level_groups[i, lv] > n_forecast_patients) ? 1
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

// Flat visit-indexed observed-PSA covariate for visit-gated 0->1 mode.
// Declared here so ALL models have the symbol in scope.
// Zero-sized when not in visit-gated observed mode — never accessed.
// Populated by _psa_observed_covar_transformed_data.stan (PSA models only).
// Indexed identically to log_psa_values[v]: visit v in [1, sum(n_patient_visits)].
// Only psa_measured[v]==1 entries are meaningful; others are 0 (never used).
// Not allocated when enable_ms_visit_gated_latent_01=1 (latent PSA used instead).
vector[compute_ms_obs_psa_covar_size(
  enable_ms_visit_gated_01, enable_ms_visit_gated_latent_01, size(t_patient_visits)
)] ms_obs_psa_covar_flat;

