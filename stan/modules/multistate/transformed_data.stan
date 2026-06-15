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

// --- Derived Censoring Indicators ---
// Logical consequences of ms_final_state, ms_time_01, and ms_time_32.
// Factored into derive_ms_censoring_indicators() so the same code path is
// exercised by Stan-level tests.
array[n_patients] int ms_censored_02;
array[n_patients] int ms_censored_12;
array[n_patients] int ms_censored_32;
(ms_censored_02, ms_censored_12, ms_censored_32) = derive_ms_censoring_indicators(
  ms_final_state, ms_time_01, ms_time_32
);

// Guard: the (level, velocity) basis is a pure 2-vs-3 feature mode switch.
// enable_ms_velocity_basis must be paired with the matching feature count, or
// the coefficient arrays mis-size relative to what the builders emit.
if (enable_ms_pop_time_varying_cov) {
  if (enable_ms_velocity_basis && n_time_varying_covar != 2)
    fatal_error("enable_ms_velocity_basis=1 requires n_time_varying_covar=2 (level, velocity); got ",
                n_time_varying_covar);
  if (!enable_ms_velocity_basis && n_time_varying_covar == 2)
    fatal_error("n_time_varying_covar=2 with enable_ms_velocity_basis=0 is ambiguous; ",
                "set enable_ms_velocity_basis=1 for the (level, velocity) basis or use the 3-feature legacy basis.");
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

// share_dead_gp_shape constraints
if (share_dead_gp_shape) {
  if (ms_time_scale_12 != 0)
    fatal_error("share_dead_gp_shape requires ms_time_scale_12=0 (Markov); extended mode not supported");
  if (!need_12_t_gp)
    fatal_error("share_dead_gp_shape requires need_12_t_gp=1");
  if (!enable_ms_02)
    fatal_error("share_dead_gp_shape requires enable_ms_02=1");
}

// --- Intercept Slot Indexing (decomposed schema) ---
// N_TRANS = 6 intercept slots over the additive hazard-intercept channels.
// Order matches flags.stan (enable_ms_level_gp / ms_level_intercept_mode /
// ms_level_intercept_corr_group rows).
int N_MS_INTERCEPT_SLOTS = 6;
int MS_SLOT_01   = 1;
int MS_SLOT_02   = 2;
int MS_SLOT_03   = 3;
int MS_SLOT_12_S = 4;
int MS_SLOT_12_T = 5;
int MS_SLOT_32   = 6;

// --- Per-slot "active" flags (transition gate) ---
// A slot contributes only when its transition channel is live. The 12_s / 12_t
// slots are gated by the time-scale-derived need flags (B1), not enable_ms_12.
array[N_MS_INTERCEPT_SLOTS] int ms_slot_active;
ms_slot_active[MS_SLOT_01]   = enable_ms_01;
ms_slot_active[MS_SLOT_02]   = enable_ms_02;
ms_slot_active[MS_SLOT_03]   = enable_ms_03;
ms_slot_active[MS_SLOT_12_S] = need_12_s_gp;
ms_slot_active[MS_SLOT_12_T] = need_12_t_gp;
ms_slot_active[MS_SLOT_32]   = enable_ms_32;

// --- Reconstruct legacy 0-4 baseline-mode vector per slot ---
// Each slot reconstructs the legacy mode vector the downstream sizing and
// residual/prior dispatch machinery still consumes. When corr_group is all-zero
// and every slot shares the same decomposed config (the bit-identical legacy
// path), these vectors all equal the old enable_ms_level_baseline_hazard.
array[N_MS_INTERCEPT_SLOTS, n_levels] int ms_legacy_mode;
for (k in 1:N_MS_INTERCEPT_SLOTS) {
  ms_legacy_mode[k] = ms_reconstruct_legacy_mode(
    n_levels, ms_slot_active[k],
    enable_ms_level_gp[k], ms_level_intercept_mode[k], 1
  );
}

// --- Validation: enable_ms_32 requires enable_ms_03 (was B3 strict check) ---
if (enable_ms_32 && !enable_ms_03)
  fatal_error("enable_ms_32 requires enable_ms_03");

// --- Validation: corr_group members must be RE, and share a group partition ---
for (lv in 1:n_levels) {
  // collect distinct positive corr groups at this level and validate members
  for (k in 1:N_MS_INTERCEPT_SLOTS) {
    int g = ms_level_intercept_corr_group[k, lv];
    if (g < 0) fatal_error("ms_level_intercept_corr_group[", k, ",", lv, "] must be >= 0");
    // Codes are bounded by the slot count: there can be at most N_MS_INTERCEPT_SLOTS
    // distinct groups at a level, so the block-discovery scan over 1..N covers all
    // realizable codes. Reject larger codes so a typo can't silently drop a block.
    if (g > N_MS_INTERCEPT_SLOTS)
      fatal_error("ms_level_intercept_corr_group[", k, ",", lv, "] = ", g,
                  " exceeds N_MS_INTERCEPT_SLOTS (", N_MS_INTERCEPT_SLOTS, ")");
    if (g > 0) {
      // member must be RE non-centered at this level when active. The MVN
      // assembly u = diag_pre_multiply(sigma, L) * z is itself the NCP form, so
      // only decomposed intercept_mode == 2 (RE-NCP) is admissible; RE-CP (3) and
      // GP (enable_ms_level_gp) would double-parameterize the same channel.
      if (ms_slot_active[k]) {
        if (ms_level_intercept_mode[k, lv] != 2) {
          fatal_error("corr_group member slot ", k, " at level ", lv,
                      " must have ms_level_intercept_mode == 2 (RE-NCP); got ",
                      ms_level_intercept_mode[k, lv]);
        }
        if (enable_ms_level_gp[k, lv]) {
          fatal_error("corr_group member slot ", k, " at level ", lv,
                      " cannot also enable_ms_level_gp (GP residual conflicts ",
                      "with the correlated NCP intercept)");
        }
      }
      // Student-t marginals are incompatible with the MVN/Gaussian-copula path.
      if (enable_student_t_hierarchy && ms_slot_active[k]) {
        fatal_error("corr_group member slot ", k, " at level ", lv,
                    " cannot use enable_student_t_hierarchy (correlated blocks ",
                    "require Gaussian marginals)");
      }
    }
  }
}

// ============================================================================
// CORRELATED INTERCEPT BLOCKS (Phase 2 — cross-transition frailty)
// ============================================================================
// A correlation "block" is a distinct positive corr_group code at a given level
// whose ACTIVE members number >= 2 (a singleton collapses to the scalar path).
// Each block is one MVN row per member transition, sampled via a Cholesky LKJ
// factor and an NCP std-normal matrix:
//     u = diag_pre_multiply(sigma_members, L) * z
// scattered back to each member transition's level-lv intercept vector.
//
// Stan requires uniform inner dimensions across an array of cholesky_factor_corr
// / matrix parameters, so every configured block must share the same dimension
// `ms_corr_dim` (number of active members) and the same group count
// `ms_corr_n_groups` (= n_forecast_groups_per_level at the block's level). Both
// are enforced with fatal_error below. The realistic configuration is a single
// d=2 patient-level block (0->1 correlated with 0->3); the machinery is fully
// general for any homogeneous-dimension set of blocks at any levels.

// First pass: count blocks and determine the (uniform) member dimension. A block
// is keyed by (level, positive group code); we discover distinct codes per level
// by treating the first active member slot that carries a code as the canonical
// representative and counting how many active slots share it.
int n_ms_corr_blocks = 0;
int ms_corr_dim = 0;        // uniform member count across all blocks (0 if none)
int ms_corr_n_groups = 0;   // uniform group count across all blocks (0 if none)
for (lv in 1:n_levels) {
  for (g in 1:N_MS_INTERCEPT_SLOTS) {
    // candidate group code `g` (corr_group codes are small positive ints; we
    // scan the code space using the slot range as an upper bound, since a code
    // can only be shared by <= N_MS_INTERCEPT_SLOTS members).
    int d_g = 0;
    for (k in 1:N_MS_INTERCEPT_SLOTS) {
      if (ms_slot_active[k] && ms_level_intercept_corr_group[k, lv] == g) {
        d_g += 1;
      }
    }
    if (d_g >= 2) {
      n_ms_corr_blocks += 1;
      if (ms_corr_dim == 0) {
        ms_corr_dim = d_g;
        ms_corr_n_groups = n_forecast_groups_per_level[lv];
      } else {
        if (d_g != ms_corr_dim)
          fatal_error("all correlated intercept blocks must share the same ",
                      "member dimension; found ", ms_corr_dim, " and ", d_g);
        if (n_forecast_groups_per_level[lv] != ms_corr_n_groups)
          fatal_error("all correlated intercept blocks must share the same ",
                      "group count; found ", ms_corr_n_groups, " and ",
                      n_forecast_groups_per_level[lv]);
      }
    }
  }
}

// Second pass: record each block's level and member slots (in slot order),
// and build a reverse lookup so each transition's intercept-assembly branch can
// ask "am I a correlated member at this level, and if so which (block, row)?"
//   ms_corr_member_block[k, lv]  = block index b (0 if slot k @ lv not correlated)
//   ms_corr_member_row[k, lv]    = MVN row m within that block (0 otherwise)
array[n_ms_corr_blocks] int ms_corr_block_level;
array[n_ms_corr_blocks, ms_corr_dim] int ms_corr_block_member_slot;
array[N_MS_INTERCEPT_SLOTS, n_levels] int ms_corr_member_block =
  rep_array(0, N_MS_INTERCEPT_SLOTS, n_levels);
array[N_MS_INTERCEPT_SLOTS, n_levels] int ms_corr_member_row =
  rep_array(0, N_MS_INTERCEPT_SLOTS, n_levels);
{
  int b = 0;
  for (lv in 1:n_levels) {
    for (g in 1:N_MS_INTERCEPT_SLOTS) {
      int d_g = 0;
      for (k in 1:N_MS_INTERCEPT_SLOTS) {
        if (ms_slot_active[k] && ms_level_intercept_corr_group[k, lv] == g) d_g += 1;
      }
      if (d_g >= 2) {
        b += 1;
        ms_corr_block_level[b] = lv;
        int m = 0;
        for (k in 1:N_MS_INTERCEPT_SLOTS) {
          if (ms_slot_active[k] && ms_level_intercept_corr_group[k, lv] == g) {
            m += 1;
            ms_corr_block_member_slot[b, m] = k;
            ms_corr_member_block[k, lv] = b;
            ms_corr_member_row[k, lv]   = m;
          }
        }
      }
    }
  }
}

// --- Per-Slot Level Baseline Hazard Flags (B2, generalized) ---
// GP mask / group counts / position arrays, one set per intercept slot.
// strict=1: fatal_error on invalid input; is_valid sentinel discarded.
int ms_flag_b2_valid;
array[N_MS_INTERCEPT_SLOTS] int any_re_level_slot;
array[N_MS_INTERCEPT_SLOTS, n_levels] int ms_level_baseline_is_gp_slot;
array[N_MS_INTERCEPT_SLOTS] int n_gp_groups_ms_baseline_slot;
array[N_MS_INTERCEPT_SLOTS, n_levels + 1] int gp_level_pos_ms_baseline_slot;
array[N_MS_INTERCEPT_SLOTS, n_levels + 1] int enabled_level_pos_ms_baseline_slot;
array[N_MS_INTERCEPT_SLOTS] int n_enabled_groups_ms_baseline_slot;
// raw/cp buckets per slot
array[N_MS_INTERCEPT_SLOTS] int n_raw_groups_ms_baseline_slot;
array[N_MS_INTERCEPT_SLOTS] int n_cp_groups_ms_baseline_slot;
array[N_MS_INTERCEPT_SLOTS, n_levels + 1] int raw_level_pos_ms_baseline_slot;
array[N_MS_INTERCEPT_SLOTS, n_levels + 1] int cp_level_pos_ms_baseline_slot;
for (k in 1:N_MS_INTERCEPT_SLOTS) {
  int valid_k;
  array[n_levels] int is_gp_k;
  int n_gp_k;
  array[n_levels + 1] int gp_pos_k;
  array[n_levels + 1] int en_pos_k;
  (valid_k, any_re_level_slot[k], is_gp_k, n_gp_k, gp_pos_k, en_pos_k) =
    compute_ms_level_baseline_flags(
      n_levels, n_forecast_groups_per_level, ms_legacy_mode[k], 1
    );
  ms_level_baseline_is_gp_slot[k]      = is_gp_k;
  n_gp_groups_ms_baseline_slot[k]      = n_gp_k;
  gp_level_pos_ms_baseline_slot[k]     = gp_pos_k;
  enabled_level_pos_ms_baseline_slot[k] = en_pos_k;
  n_enabled_groups_ms_baseline_slot[k] = compute_n_enabled_groups(
    n_forecast_groups_per_level, ms_legacy_mode[k]
  );
  int n_raw_k;
  array[n_levels + 1] int raw_pos_k;
  int n_cp_k;
  array[n_levels + 1] int cp_pos_k;
  (n_raw_k, raw_pos_k, n_cp_k, cp_pos_k) =
    split_cp_ncp_pos(n_levels, n_forecast_groups_per_level, ms_legacy_mode[k]);
  n_raw_groups_ms_baseline_slot[k]  = n_raw_k;
  n_cp_groups_ms_baseline_slot[k]   = n_cp_k;
  raw_level_pos_ms_baseline_slot[k] = raw_pos_k;
  cp_level_pos_ms_baseline_slot[k]  = cp_pos_k;
}
ms_flag_b2_valid = 1;

// --- Per-transition scalar aliases ---
// Named per-transition views into the per-slot arrays, used by parameters.stan,
// transformed_parameters.stan, and priors.stan. (Each slot now has its own
// enabled / GP / raw / cp counts; nothing is shared across transitions.)
int n_enabled_groups_ms_baseline_01   = n_enabled_groups_ms_baseline_slot[MS_SLOT_01];
int n_gp_groups_ms_baseline_01        = n_gp_groups_ms_baseline_slot[MS_SLOT_01];
int n_enabled_groups_ms_baseline_02   = n_enabled_groups_ms_baseline_slot[MS_SLOT_02];
int n_gp_groups_ms_baseline_02        = n_gp_groups_ms_baseline_slot[MS_SLOT_02];
int n_enabled_groups_ms_baseline_12_s = n_enabled_groups_ms_baseline_slot[MS_SLOT_12_S];
int n_gp_groups_ms_baseline_12_s      = n_gp_groups_ms_baseline_slot[MS_SLOT_12_S];
int n_enabled_groups_ms_baseline_12_t = n_enabled_groups_ms_baseline_slot[MS_SLOT_12_T];
int n_gp_groups_ms_baseline_12_t      = n_gp_groups_ms_baseline_slot[MS_SLOT_12_T];
int n_enabled_groups_ms_baseline_03   = n_enabled_groups_ms_baseline_slot[MS_SLOT_03];
int n_gp_groups_ms_baseline_03        = n_gp_groups_ms_baseline_slot[MS_SLOT_03];
int n_enabled_groups_ms_baseline_32   = n_enabled_groups_ms_baseline_slot[MS_SLOT_32];
int n_gp_groups_ms_baseline_32        = n_gp_groups_ms_baseline_slot[MS_SLOT_32];

int n_raw_groups_ms_baseline_01   = n_raw_groups_ms_baseline_slot[MS_SLOT_01];
int n_cp_groups_ms_baseline_01    = n_cp_groups_ms_baseline_slot[MS_SLOT_01];
int n_raw_groups_ms_baseline_02   = n_raw_groups_ms_baseline_slot[MS_SLOT_02];
int n_cp_groups_ms_baseline_02    = n_cp_groups_ms_baseline_slot[MS_SLOT_02];
int n_raw_groups_ms_baseline_12_s = n_raw_groups_ms_baseline_slot[MS_SLOT_12_S];
int n_cp_groups_ms_baseline_12_s  = n_cp_groups_ms_baseline_slot[MS_SLOT_12_S];
int n_raw_groups_ms_baseline_12_t = n_raw_groups_ms_baseline_slot[MS_SLOT_12_T];
int n_cp_groups_ms_baseline_12_t  = n_cp_groups_ms_baseline_slot[MS_SLOT_12_T];
int n_raw_groups_ms_baseline_03   = n_raw_groups_ms_baseline_slot[MS_SLOT_03];
int n_cp_groups_ms_baseline_03    = n_cp_groups_ms_baseline_slot[MS_SLOT_03];
int n_raw_groups_ms_baseline_32   = n_raw_groups_ms_baseline_slot[MS_SLOT_32];
int n_cp_groups_ms_baseline_32    = n_cp_groups_ms_baseline_slot[MS_SLOT_32];

// --- Enabled Group Counts for Covariate Slopes ---
int n_enabled_groups_ms_slope = compute_n_enabled_groups(
  n_forecast_groups_per_level, enable_ms_level_cov
);

// Position array for enabled slope levels
array[n_levels + 1] int enabled_level_pos_ms_slope = create_enabled_pos(
  n_forecast_groups_per_level, enable_ms_level_cov
);

// --- Slope bucket routing (same pattern as intercepts) ---
// Slope parameterization follows the intercept mode at that level;
// slope is only included when enable_ms_level_cov[lv] == 1.
// Slopes follow the intercept mode at that level. With the legacy single flag
// removed, the slope mode is sourced from the 0->1 slot's reconstructed legacy
// mode vector (slots share one config on the bit-identical legacy path, so this
// reproduces the old enable_ms_level_baseline_hazard[lv] exactly). Slopes are
// not themselves per-transition in the raw/cp routing, so a single
// representative slot drives ms_slope_mode.
array[n_levels] int ms_slope_mode;
for (lv in 1:n_levels) {
  ms_slope_mode[lv] = enable_ms_level_cov[lv] ? ms_legacy_mode[MS_SLOT_01, lv] : 0;
}
int n_raw_groups_ms_slope_shared;
int n_cp_groups_ms_slope_shared;
array[n_levels + 1] int raw_level_pos_ms_slope_shared;
array[n_levels + 1] int cp_level_pos_ms_slope_shared;
(n_raw_groups_ms_slope_shared, raw_level_pos_ms_slope_shared,
 n_cp_groups_ms_slope_shared,  cp_level_pos_ms_slope_shared) =
  split_cp_ncp_pos(n_levels, n_forecast_groups_per_level, ms_slope_mode);

// --- 0->1 Baseline Log-Time Trend Centering Constant ---
// g(t) = log(t) - ms_log_t_centering, where ms_log_t_centering = mean(log(1..max_all_t)).
// Centering keeps the slope mean-zero and preserves the GP intercept's meaning.
// Always computed (cheap); only consumed when enable_ms_baseline_trend_01 == 1.
real ms_log_t_centering;
{
  real log_t_sum = 0;
  for (t in 1:max_all_t) log_t_sum += log(t);
  ms_log_t_centering = log_t_sum / max_all_t;
}

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
// Baseline hazard level indices, one gather table per intercept slot (which
// levels are enabled now varies by transition). Indexed [slot, patient, level].
array[N_MS_INTERCEPT_SLOTS, n_patients, n_levels] int patient_ms_baseline_flat_idx_slot;
{
  for (k in 1:N_MS_INTERCEPT_SLOTS) {
    for (i in 1:n_patients) {
      for (lv in 1:n_levels) {
        if (ms_legacy_mode[k, lv]) {
          patient_ms_baseline_flat_idx_slot[k, i, lv] =
            (lv == n_levels && patient_level_groups[i, lv] > n_forecast_patients) ? 1
            : get_global_group_idx(enabled_level_pos_ms_baseline_slot[k], lv, patient_level_groups[i, lv]);
        } else {
          patient_ms_baseline_flat_idx_slot[k, i, lv] = 1;
        }
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

// Compute ms_patient_idx: compacted index of patients contributing to MS likelihood
array[n_patients] int<lower=0, upper=1> ms_is_target = rep_array(0, n_patients);
int n_ms_patients = 0;
if (ms_split_level == 0) {
  for (j in 1:n_forecast_patients) {
    ms_is_target[forecast_patient_idx[j]] = 1;
  }
  n_ms_patients = n_forecast_patients;
} else {
  for (j in 1:n_forecast_patients) {
    int p = forecast_patient_idx[j];
    int grp = patient_level_groups[p, ms_split_level];
    for (k in 1:n_ms_target_groups) {
      if (grp == ms_target_groups[k]) {
        ms_is_target[p] = 1;
        n_ms_patients += 1;
        break;
      }
    }
  }
}

if (fit_multistate_data == 1 && n_ms_patients == 0) {
  fatal_error("MS likelihood enabled but n_ms_patients == 0 (check ms_split_level/ms_target_groups)");
}

array[n_ms_patients] int ms_patient_idx;
{
  int write_idx = 1;
  for (j in 1:n_forecast_patients) {
    int p = forecast_patient_idx[j];
    if (ms_is_target[p] == 1) {
      ms_patient_idx[write_idx] = p;
      write_idx += 1;
    }
  }
}

// Flat visit-indexed observed-PSA covariate for visit-gated 0->1 mode.
// Declared here so ALL models have the symbol in scope.
// Zero-sized when not in visit-gated observed mode — never accessed.
// Populated by _psa_observed_covar_transformed_data.stan (PSA models only).
// Indexed identically to log_psa_values[v]: visit v in [1, sum(n_patient_visits)].
// Only psa_measured[v]==1 entries are meaningful; others are 0 (never used).
// Not allocated when enable_ms_visit_gated_latent_01=1 (latent PSA used instead).
vector[compute_ms_obs_visit_covar_size(
  enable_ms_visit_gated_01, enable_ms_visit_gated_latent_01, size(t_patient_visits)
)] ms_obs_visit_covar_flat;

