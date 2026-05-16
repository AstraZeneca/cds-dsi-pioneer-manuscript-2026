// tr/_sd_subhierarchy_transformed_parameters.stan
// Build per-group SD vector tr_sd_intercept_pergroup[lv] for each level.
// When sub-hierarchy inactive at lv: constant fill with tr_sd_level_intercept[lv] → bit-exact.
// When active: sum log-deviations additively, then exp once.

array[n_levels] vector[max(n_forecast_groups_per_level)] tr_sd_intercept_pergroup;
for (lv in 1:n_levels) {
  int ngroups_lv = n_forecast_groups_per_level[lv];
  if (!has_sd_subhierarchy_tr_intercept[lv]) {
    // Inactive: constant fill with existing scalar → bit-exact with pre-feature code.
    tr_sd_intercept_pergroup[lv][1:ngroups_lv] =
      rep_vector(tr_sd_level_intercept[lv], ngroups_lv);
  } else {
    // Active: additive log-scale assembly.
    vector[ngroups_lv] log_sd_g = rep_vector(
      tr_log_sd_level_intercept_pop[pop_idx_tr_intercept[lv]],
      ngroups_lv
    );
    for (sub_lv in 1:(lv - 1)) {
      int m = enable_sd_level_intercept_mode_tr[lv, sub_lv];
      if (m == LEVEL_MODE_NONE) continue;

      int n_sub_groups = n_forecast_groups_per_level[sub_lv];
      vector[n_sub_groups] log_sd_delta;

      if (m == LEVEL_MODE_RE_CP) {
        int c_lo = cp_pos_tr_log_sd_intercept[lv, sub_lv];
        int c_hi = cp_pos_tr_log_sd_intercept[lv, sub_lv + 1] - 1;
        log_sd_delta = tr_cp_log_sd_level_intercept[c_lo:c_hi];
      } else {
        // FE or RE
        int r_lo = raw_pos_tr_log_sd_intercept[lv, sub_lv];
        int r_hi = raw_pos_tr_log_sd_intercept[lv, sub_lv + 1] - 1;
        real hyperscale = (m == LEVEL_MODE_FE)
          ? tr_fe_sd_hyperscale_level_intercept[lv, sub_lv]
          : tr_sd_hyperscale_level_intercept_raw[hyperscale_idx_tr_intercept[lv, sub_lv]];
        log_sd_delta = hyperscale * tr_raw_log_sd_level_intercept[r_lo:r_hi];
      }

      // Scatter log_sd_delta from sub-level-ℓ groups to location-level-L groups via subgroup_idx.
      for (g_lv in 1:ngroups_lv) {
        int g_sub = subgroup_idx_tr_intercept[lv, sub_lv, g_lv];
        log_sd_g[g_lv] += log_sd_delta[g_sub];
      }
    }
    tr_sd_intercept_pergroup[lv][1:ngroups_lv] = exp(log_sd_g);
  }
}
