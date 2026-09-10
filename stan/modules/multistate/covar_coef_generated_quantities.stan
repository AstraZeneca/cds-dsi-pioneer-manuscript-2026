// ============================================================================
// BACK-TRANSFORMED TIME-INVARIANT COVARIATE COEFFICIENTS
// ============================================================================
// Recovers original-space coefficients from QR-space parameters, for both the
// population-level coefficients and (when enabled) per-arm (level-2) random
// slopes.  Applies to all enabled transitions: 0→1, 0→2, 0→3, 1→2.
//
// Population back-transform (all transitions):
//   beta = R^{-1} * theta_qr
//   Using: R^{-1} * v = (mdivide_right_tri_low(v', R'))' ,
//   where R' is lower triangular.
//
// Arm-level back-transform (when enable_ms_level_cov[2] == 1):
//   scaled_qr[arm] = raw_level_slope[arm] .* sd_level_slope[2]  (QR space)
//   arm_coef[arm]  = R^{-1} * scaled_qr[arm]                    (original space)
//   arm_level_coef_sd = sd of arm_coef across arms, per covariate
// ============================================================================

// ── Population-level coefficients ────────────────────────────────────────────

vector[enable_ms_01 && enable_ms_pop_time_invariant_cov && n_time_invariant_covar > 0
         ? n_time_invariant_covar : 0] time_invariant_coef_01;
if (enable_ms_01 && enable_ms_pop_time_invariant_cov && n_time_invariant_covar > 0) {
  time_invariant_coef_01 = (mdivide_right_tri_low(
    time_invariant_coef_qr_01', R_covar_design_matrix'
  ))';
}

vector[enable_ms_02 && enable_ms_pop_time_invariant_cov && n_time_invariant_covar > 0
         ? n_time_invariant_covar : 0] time_invariant_coef_02;
if (enable_ms_02 && enable_ms_pop_time_invariant_cov && n_time_invariant_covar > 0) {
  time_invariant_coef_02 = (mdivide_right_tri_low(
    time_invariant_coef_qr_02', R_covar_design_matrix'
  ))';
}

vector[enable_ms_12 && enable_ms_pop_time_invariant_cov && n_time_invariant_covar > 0
         ? n_time_invariant_covar : 0] time_invariant_coef_12;
if (enable_ms_12 && enable_ms_pop_time_invariant_cov && n_time_invariant_covar > 0) {
  time_invariant_coef_12 = (mdivide_right_tri_low(
    time_invariant_coef_qr_12', R_covar_design_matrix'
  ))';
}

vector[enable_ms_03 && enable_ms_pop_time_invariant_cov && n_time_invariant_covar > 0
         ? n_time_invariant_covar : 0] time_invariant_coef_03;
if (enable_ms_03 && enable_ms_pop_time_invariant_cov && n_time_invariant_covar > 0) {
  time_invariant_coef_03 = (mdivide_right_tri_low(
    time_invariant_coef_qr_03', R_covar_design_matrix'
  ))';
}

// ── Arm-level (level-2) random slopes ────────────────────────────────────────
// Only active when enable_ms_level_cov[2] == 1. Level 2 = arm level in the
// three-level hierarchy (population → arm → patient).
// enabled_level_pos_ms_slope[2] gives the 1-based start index of level-2 arms
// in raw_level_slope_0X / raw_level_slope_12.
//
// Declaration sizes use data-block variables directly (Stan requires data
// variables in top-level GQ size expressions, not locally-computed ints).

// n_groups_per_level[2] is the arm count; 0 when level-cov disabled or
// n_time_invariant_covar == 0 (size is data-based ternary, which Stan allows).
matrix[
  enable_ms_pop_time_invariant_cov && n_time_invariant_covar > 0
    && n_levels >= 2 && enable_ms_level_cov[2] ? n_groups_per_level[2] : 0,
  enable_ms_pop_time_invariant_cov && n_time_invariant_covar > 0
    && n_levels >= 2 && enable_ms_level_cov[2] ? n_time_invariant_covar : 0
] arm_coef_01, arm_coef_02, arm_coef_12;

vector[
  enable_ms_pop_time_invariant_cov && n_time_invariant_covar > 0
    && n_levels >= 2 && enable_ms_level_cov[2] ? n_time_invariant_covar : 0
] arm_level_coef_01_sd, arm_level_coef_02_sd, arm_level_coef_12_sd;

if (enable_ms_pop_time_invariant_cov && n_time_invariant_covar > 0
    && n_levels >= 2 && enable_ms_level_cov[2]) {
  int n_slope_arms = n_groups_per_level[2];
  int arm_start    = enabled_level_pos_ms_slope[2];

  for (arm in 1:n_slope_arms) {
    int row = arm_start + arm - 1;

    if (enable_ms_01) {
      // raw_level_slope_01[row] is row_vector; sd_level_slope_01[2] is vector.
      // Scale in row_vector space then back-transform via R^{-1}.
      row_vector[n_time_invariant_covar] sq01 =
        raw_level_slope_01[row] .* sd_level_slope_01[2]';
      arm_coef_01[arm] = mdivide_right_tri_low(sq01, R_covar_design_matrix');
    }

    if (enable_ms_02) {
      row_vector[n_time_invariant_covar] sq02 =
        raw_level_slope_02[row] .* sd_level_slope_02[2]';
      arm_coef_02[arm] = mdivide_right_tri_low(sq02, R_covar_design_matrix');
    }

    if (enable_ms_12) {
      row_vector[n_time_invariant_covar] sq12 =
        raw_level_slope_12[row] .* sd_level_slope_12[2]';
      arm_coef_12[arm] = mdivide_right_tri_low(sq12, R_covar_design_matrix');
    }
  }

  for (c in 1:n_time_invariant_covar) {
    if (enable_ms_01) arm_level_coef_01_sd[c] = sd(arm_coef_01[:, c]);
    if (enable_ms_02) arm_level_coef_02_sd[c] = sd(arm_coef_02[:, c]);
    if (enable_ms_12) arm_level_coef_12_sd[c] = sd(arm_coef_12[:, c]);
  }
}
