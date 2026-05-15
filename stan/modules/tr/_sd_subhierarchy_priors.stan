// tr/_sd_subhierarchy_priors.stan
// Priors for SD sub-hierarchy parameters (issue #110).
// All statements are vacuously true when the mode matrix is all-NONE (size-0 vectors).

// Population log-SD per location level
for (lv in 1:n_levels) {
  if (has_sd_subhierarchy_tr_intercept[lv]) {
    int idx = pop_idx_tr_intercept[lv];
    tr_log_sd_level_intercept_pop[idx]
      ~ normal(tr_log_sd_level_intercept_pop_mean[lv], tr_log_sd_level_intercept_pop_sd[lv]);
  }
}

// Free hyperscales (half-normal via <lower=0>)
for (lv in 1:n_levels) {
  for (sub_lv in 1:(lv - 1)) {
    int m = enable_sd_level_intercept_mode_tr[lv, sub_lv];
    if (m == LEVEL_MODE_RE || m == LEVEL_MODE_RE_CP) {
      int idx = hyperscale_idx_tr_intercept[lv, sub_lv];
      tr_sd_hyperscale_level_intercept_raw[idx]
        ~ normal(0, tr_sd_hyperscale_level_intercept_sd[lv, sub_lv]);
    }
  }
}

// NCP raw deviations (FE and RE share raw bucket)
if (n_raw_groups_tr_log_sd_intercept > 0) {
  tr_raw_log_sd_level_intercept ~ std_normal();
}

// CP deviations — sampled at hyperscale directly
for (lv in 1:n_levels) {
  if (!has_sd_subhierarchy_tr_intercept[lv]) continue;
  for (sub_lv in 1:(lv - 1)) {
    if (enable_sd_level_intercept_mode_tr[lv, sub_lv] == LEVEL_MODE_RE_CP) {
      int c_lo = cp_pos_tr_log_sd_intercept[lv, sub_lv];
      int c_hi = cp_pos_tr_log_sd_intercept[lv, sub_lv + 1] - 1;
      real hs = tr_sd_hyperscale_level_intercept_raw[hyperscale_idx_tr_intercept[lv, sub_lv]];
      tr_cp_log_sd_level_intercept[c_lo:c_hi] ~ normal(0, hs);
    }
  }
}
