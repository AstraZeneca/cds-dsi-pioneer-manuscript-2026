// tr/priors.stan — active priors for total rate module
// Multi-level hierarchy: loop over all levels for unified prior structure
// Uses enabled position arrays for efficient indexing into compacted parameter arrays

// Population intercept
tr_loc_pop ~ normal(tr_loc_pop_mean, tr_loc_pop_sd);

// Population covariate effects
if (enable_pop_cov_tr) {
  tr_coef_qr_pop ~ normal(tr_coef_qr_pop_mean, tr_coef_qr_pop_sd);
}

// ===== UNIFIED LOOP OVER ALL LEVELS =====
{
  int sd_idx = 0;
  for (lv in 1:n_levels) {
    // SD prior only for RE levels (mode=2); FE levels use fixed hyperparameter
    if (enable_level_intercept_tr[lv] == LEVEL_MODE_RE) {
      sd_idx += 1;
      tr_sd_level_intercept_raw[sd_idx] ~ normal(0, tr_sd_level_intercept_sd[lv]);
    }
    if (n_covar > 0) {
      tr_sd_level_slope[lv] ~ normal(0, tr_sd_level_slope_sd[lv]);
    }

    // Intercept raw effects - only apply prior to enabled levels (FE and RE both)
    // (parameter array is sized by enabled groups only)
    if (enable_level_intercept_tr[lv]) {
      int lv_start = enabled_level_pos_tr_intercept[lv];
      int lv_end = enabled_level_pos_tr_intercept[lv + 1] - 1;
      tr_raw_level_intercept[lv_start:lv_end] ~ std_normal();
    }

    // Slope raw effects - only apply prior to enabled levels
    if (enable_level_cov_tr[lv] && n_covar > 0) {
      int lv_start = enabled_level_pos_tr_slope[lv];
      int lv_end = enabled_level_pos_tr_slope[lv + 1] - 1;
      to_vector(tr_raw_level_slope[lv_start:lv_end, :]) ~ std_normal();
    }
  }
}

// Patient-level process noise priors - only when feature is enabled
if (enable_patient_process_noise_tr) {
	to_vector(tr_raw_patient_process_noise) ~ std_normal();
	tr_log_sd_pop_process_noise[1] ~ normal(tr_log_sd_pop_process_noise_mean, tr_log_sd_pop_process_noise_sd);
	tr_logit_phi_pop_process_noise[1] ~ normal(tr_logit_phi_pop_process_noise_mean, tr_logit_phi_pop_process_noise_sd);
	tr_sd_patient_log_sd_process_noise[1] ~ normal(0, tr_log_sd_patient_process_noise_sd);
	tr_sd_patient_phi_process_noise[1] ~ normal(0, tr_phi_patient_process_noise_sd);
	tr_raw_patient_log_sd_process_noise ~ std_normal();
	tr_raw_patient_phi_process_noise ~ std_normal();
}

// Population-level time-varying process noise priors
if (enable_pop_process_noise_tr) {
	tr_raw_pop_process_noise ~ std_normal();
	tr_log_sd_pop_process_noise_pop[1] ~ normal(tr_log_sd_pop_process_noise_pop_mean, tr_log_sd_pop_process_noise_pop_sd);
	tr_logit_phi_pop_process_noise_pop[1] ~ normal(tr_logit_phi_pop_process_noise_pop_mean, tr_logit_phi_pop_process_noise_pop_sd);
}
