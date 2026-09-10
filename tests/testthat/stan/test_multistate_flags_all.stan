// Test harness for B1-B4 multistate flag functions.
// All four function families are tested in a single batch-loop design.
functions {
  #include "util.stanfunctions"
  #include "pos.stanfunctions"
  #include "pfs.stanfunctions"
  #include "multistate.stanfunctions"
}
data {
  // ── B1: compute_ms_time_scale_flags ───────────────────────────────────────
  int<lower=0> n_combos_b1;
  array[n_combos_b1] int enable_ms_12_combos;
  array[n_combos_b1] int ms_time_scale_12_combos;

  // ── B2: compute_ms_level_baseline_flags ───────────────────────────────────
  int<lower=1> n_levels_b2;
  int<lower=0> n_combos_b2;
  array[n_levels_b2] int n_forecast_groups_b2;
  array[n_combos_b2, n_levels_b2] int enable_ms_level_baseline_combos;

  // ── B3: compute_ms_transition_group_counts ────────────────────────────────
  int<lower=0> n_enabled_b3;
  int<lower=0> n_gp_b3;
  int<lower=0> n_combos_b3;
  array[n_combos_b3] int enable_ms_01_combos_b3;
  array[n_combos_b3] int enable_ms_02_combos_b3;
  array[n_combos_b3] int need_12_s_gp_combos_b3;
  array[n_combos_b3] int need_12_t_gp_combos_b3;
  array[n_combos_b3] int enable_ms_03_combos_b3;
  array[n_combos_b3] int enable_ms_32_combos_b3;

  // ── B4: compute_ms_obs_visit_covar_size ─────────────────────────────────────
  int<lower=0> n_combos_b4;
  int<lower=0> n_visits_b4;
  array[n_combos_b4] int enable_ms_visit_gated_01_combos_b4;
  array[n_combos_b4] int enable_ms_visit_gated_latent_01_combos_b4;
}
generated quantities {
  // ── B1 outputs ─────────────────────────────────────────────────────────────
  array[n_combos_b1] int out_b1_is_valid;
  array[n_combos_b1] int out_b1_need_s;
  array[n_combos_b1] int out_b1_need_t;
  array[n_combos_b1] int out_b1_t_int;

  for (ci in 1:n_combos_b1) {
    int is_valid; int need_s; int need_t; int t_int;
    (is_valid, need_s, need_t, t_int) = compute_ms_time_scale_flags(
      enable_ms_12_combos[ci],
      ms_time_scale_12_combos[ci],
      0  // strict=0: don't fatal-error on bad input
    );
    out_b1_is_valid[ci] = is_valid;
    out_b1_need_s[ci]   = need_s;
    out_b1_need_t[ci]   = need_t;
    out_b1_t_int[ci]    = t_int;
  }

  // ── B2 outputs ─────────────────────────────────────────────────────────────
  array[n_combos_b2] int out_b2_is_valid;
  array[n_combos_b2] int out_b2_any_re;
  array[n_combos_b2, n_levels_b2] int out_b2_is_gp;
  array[n_combos_b2] int out_b2_n_gp;
  array[n_combos_b2, n_levels_b2 + 1] int out_b2_gp_pos;
  array[n_combos_b2, n_levels_b2 + 1] int out_b2_enabled_pos;

  for (ci in 1:n_combos_b2) {
    int is_valid; int any_re;
    array[n_levels_b2] int is_gp;
    int n_gp;
    array[n_levels_b2 + 1] int gp_pos;
    array[n_levels_b2 + 1] int enabled_pos;

    array[n_levels_b2] int flags_ci;
    for (lv in 1:n_levels_b2) flags_ci[lv] = enable_ms_level_baseline_combos[ci, lv];

    (is_valid, any_re, is_gp, n_gp, gp_pos, enabled_pos) =
      compute_ms_level_baseline_flags(
        n_levels_b2,
        n_forecast_groups_b2,
        flags_ci,
        0  // strict=0
      );

    out_b2_is_valid[ci] = is_valid;
    out_b2_any_re[ci]   = any_re;
    out_b2_n_gp[ci]     = n_gp;
    for (lv in 1:n_levels_b2) {
      out_b2_is_gp[ci, lv] = is_gp[lv];
    }
    for (k in 1:(n_levels_b2 + 1)) {
      out_b2_gp_pos[ci, k]     = gp_pos[k];
      out_b2_enabled_pos[ci, k] = enabled_pos[k];
    }
  }

  // ── B3 outputs ─────────────────────────────────────────────────────────────
  array[n_combos_b3] int out_b3_is_valid;
  array[n_combos_b3] int out_b3_n_en_01;
  array[n_combos_b3] int out_b3_n_gp_01;
  array[n_combos_b3] int out_b3_n_en_02;
  array[n_combos_b3] int out_b3_n_gp_02;
  array[n_combos_b3] int out_b3_n_en_12s;
  array[n_combos_b3] int out_b3_n_gp_12s;
  array[n_combos_b3] int out_b3_n_en_12t;
  array[n_combos_b3] int out_b3_n_gp_12t;
  array[n_combos_b3] int out_b3_n_en_03;
  array[n_combos_b3] int out_b3_n_gp_03;
  array[n_combos_b3] int out_b3_n_en_32;
  array[n_combos_b3] int out_b3_n_gp_32;

  for (ci in 1:n_combos_b3) {
    int is_valid;
    int n_en_01; int n_gp_01;
    int n_en_02; int n_gp_02;
    int n_en_12s; int n_gp_12s;
    int n_en_12t; int n_gp_12t;
    int n_en_03; int n_gp_03;
    int n_en_32; int n_gp_32;

    (is_valid,
     n_en_01, n_gp_01,
     n_en_02, n_gp_02,
     n_en_12s, n_gp_12s,
     n_en_12t, n_gp_12t,
     n_en_03, n_gp_03,
     n_en_32, n_gp_32) = compute_ms_transition_group_counts(
       enable_ms_01_combos_b3[ci],
       enable_ms_02_combos_b3[ci],
       need_12_s_gp_combos_b3[ci],
       need_12_t_gp_combos_b3[ci],
       enable_ms_03_combos_b3[ci],
       enable_ms_32_combos_b3[ci],
       n_enabled_b3,
       n_gp_b3,
       0  // strict=0
     );

    out_b3_is_valid[ci] = is_valid;
    out_b3_n_en_01[ci]  = n_en_01;
    out_b3_n_gp_01[ci]  = n_gp_01;
    out_b3_n_en_02[ci]  = n_en_02;
    out_b3_n_gp_02[ci]  = n_gp_02;
    out_b3_n_en_12s[ci] = n_en_12s;
    out_b3_n_gp_12s[ci] = n_gp_12s;
    out_b3_n_en_12t[ci] = n_en_12t;
    out_b3_n_gp_12t[ci] = n_gp_12t;
    out_b3_n_en_03[ci]  = n_en_03;
    out_b3_n_gp_03[ci]  = n_gp_03;
    out_b3_n_en_32[ci]  = n_en_32;
    out_b3_n_gp_32[ci]  = n_gp_32;
  }

  // ── B4 outputs ─────────────────────────────────────────────────────────────
  array[n_combos_b4] int out_b4_size;

  for (ci in 1:n_combos_b4) {
    out_b4_size[ci] = compute_ms_obs_visit_covar_size(
      enable_ms_visit_gated_01_combos_b4[ci],
      enable_ms_visit_gated_latent_01_combos_b4[ci],
      n_visits_b4
    );
  }
}
