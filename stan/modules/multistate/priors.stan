// ============================================================================
// Multistate Hazard Model Priors
// ============================================================================

// ============================================================================
// 0→1 TRANSITION PRIORS
// ============================================================================

if (enable_ms_01) {
  // Population-level baseline hazard GP
  log_lambda_gp_01_pop_alpha[1] ~ inv_gamma(
    log_lambda_gp_01_pop_alpha_alpha, log_lambda_gp_01_pop_alpha_beta
  );
  log_lambda_gp_01_pop_rho[1] ~ inv_gamma(
    log_lambda_gp_01_pop_rho_alpha, log_lambda_gp_01_pop_rho_beta
  );
  log_lambda_gp_01_pop_intercept[1] ~ normal(
    log_lambda_gp_01_pop_intercept_mean, log_lambda_gp_01_pop_intercept_sd
  );
  to_vector(log_lambda_gp_01_pop_eta) ~ std_normal();

  // Level-level baseline hazard GP
  // Always set priors for all levels (even disabled ones) to avoid improper posteriors
  for (lv in 1:n_levels) {
    log_lambda_gp_01_level_alpha[lv] ~ inv_gamma(
      log_lambda_gp_01_level_alpha_alpha[lv], log_lambda_gp_01_level_alpha_beta[lv]
    );
    log_lambda_gp_01_level_rho[lv] ~ inv_gamma(
      log_lambda_gp_01_level_rho_alpha[lv], log_lambda_gp_01_level_rho_beta[lv]
    );
    // SD prior for all levels when the SD array is allocated. Without this,
    // non-RE slots of log_lambda_gp_*_level_intercept_sd are free parameters
    // with no prior and drift to numerical infinity, distorting mass-matrix
    // adaptation for the rest of the model.
    if (any_re_level) {
      log_lambda_gp_01_level_intercept_sd[lv] ~ normal(
        0, log_lambda_gp_01_level_intercept_sd_sd[lv]
      );
    }
  }

  // Group-level intercepts: per-level dispatch on mode
  for (lv in 1:n_levels) {
    int mode = enable_ms_level_baseline_hazard[lv];
    if (!mode) continue;

    // RAW path: FE, RE, RE_GP sample from std_normal / student_t(nu, 0, 1)
    if (mode == LEVEL_MODE_FE || mode == LEVEL_MODE_RE || mode == LEVEL_MODE_RE_GP) {
      int r_lo = raw_level_pos_ms_baseline_shared[lv];
      int r_hi = raw_level_pos_ms_baseline_shared[lv + 1] - 1;
      if (r_hi >= r_lo) {
        if (enable_student_t_hierarchy)
          raw_log_lambda_gp_01_level_intercept[r_lo:r_hi]
            ~ student_t(ms_nu_baseline_level[lv], 0, 1);
        else
          raw_log_lambda_gp_01_level_intercept[r_lo:r_hi] ~ std_normal();
      }
    }

    // CP path: RE_CP samples at natural scale
    if (mode == LEVEL_MODE_RE_CP) {
      int c_lo = cp_level_pos_ms_baseline_shared[lv];
      int c_hi = cp_level_pos_ms_baseline_shared[lv + 1] - 1;
      if (c_hi >= c_lo) {
        if (enable_student_t_hierarchy)
          cp_log_lambda_gp_01_level_intercept[c_lo:c_hi]
            ~ student_t(ms_nu_baseline_level[lv], 0, log_lambda_gp_01_level_intercept_sd[lv]);
        else
          cp_log_lambda_gp_01_level_intercept[c_lo:c_hi]
            ~ normal(0, log_lambda_gp_01_level_intercept_sd[lv]);
      }
    }
  }
  // GP eta (GP mode only)
  if (n_gp_groups_ms_baseline_01 > 0) {
    to_vector(log_lambda_gp_01_level_eta) ~ std_normal();
  }

  // Time-varying covariate coefficients
  // 0->1 time-varying covariate prior (size 1 in visit-gated, n_time_varying_covar otherwise)
  if (enable_ms_pop_time_varying_cov && (enable_ms_visit_gated_01 || n_time_varying_covar > 0)) {
    time_varying_coef_01 ~ normal(time_varying_coef_01_mean, time_varying_coef_01_sd);
  }

  // Time-invariant covariate coefficients
  if (enable_ms_pop_time_invariant_cov && n_time_invariant_covar > 0) {
    time_invariant_coef_qr_01 ~ normal(time_invariant_coef_01_mean, time_invariant_coef_01_sd);
  }

  // Multi-level random slopes
  for (lv in 1:n_levels) {
    if (n_time_invariant_covar > 0) {
      sd_level_slope_01[lv] ~ normal(0, sd_level_slope_01_sd[lv]');
    }
    if (enable_ms_level_cov[lv] && n_time_invariant_covar > 0) {
      int mode = enable_ms_level_baseline_hazard[lv];
      // RAW path: FE, RE, RE_GP sample from std_normal / student_t(nu, 0, 1)
      if (mode == LEVEL_MODE_FE || mode == LEVEL_MODE_RE || mode == LEVEL_MODE_RE_GP) {
        int r_lo = raw_level_pos_ms_slope_shared[lv];
        int r_hi = raw_level_pos_ms_slope_shared[lv + 1] - 1;
        if (r_hi >= r_lo) {
          if (enable_student_t_hierarchy)
            to_vector(raw_level_slope_01[r_lo:r_hi, :]) ~ student_t(ms_nu_slope_level[lv], 0, 1);
          else
            to_vector(raw_level_slope_01[r_lo:r_hi, :]) ~ std_normal();
        }
      }
      // CP path: RE_CP samples at natural scale
      if (mode == LEVEL_MODE_RE_CP) {
        int c_lo = cp_level_pos_ms_slope_shared[lv];
        int c_hi = cp_level_pos_ms_slope_shared[lv + 1] - 1;
        if (c_hi >= c_lo) {
          for (k in 1:n_time_invariant_covar) {
            if (enable_student_t_hierarchy)
              cp_level_slope_01[c_lo:c_hi, k]
                ~ student_t(ms_nu_slope_level[lv], 0, sd_level_slope_01[lv, k]);
            else
              cp_level_slope_01[c_lo:c_hi, k]
                ~ normal(0, sd_level_slope_01[lv, k]);
          }
        }
      }
    }
  }
}

// ============================================================================
// 0→2 TRANSITION PRIORS
// ============================================================================

if (enable_ms_02) {
  // Population-level baseline hazard GP
  if (!share_dead_gp_shape) {
    log_lambda_gp_02_pop_alpha[1] ~ inv_gamma(
      log_lambda_gp_02_pop_alpha_alpha, log_lambda_gp_02_pop_alpha_beta
    );
    log_lambda_gp_02_pop_rho[1] ~ inv_gamma(
      log_lambda_gp_02_pop_rho_alpha, log_lambda_gp_02_pop_rho_beta
    );
    to_vector(log_lambda_gp_02_pop_eta) ~ std_normal();
  }
  log_lambda_gp_02_pop_intercept[1] ~ normal(
    log_lambda_gp_02_pop_intercept_mean, log_lambda_gp_02_pop_intercept_sd
  );

  // Level-level baseline hazard GP
  for (lv in 1:n_levels) {
    log_lambda_gp_02_level_alpha[lv] ~ inv_gamma(
      log_lambda_gp_02_level_alpha_alpha[lv], log_lambda_gp_02_level_alpha_beta[lv]
    );
    log_lambda_gp_02_level_rho[lv] ~ inv_gamma(
      log_lambda_gp_02_level_rho_alpha[lv], log_lambda_gp_02_level_rho_beta[lv]
    );
    if (any_re_level) {
      log_lambda_gp_02_level_intercept_sd[lv] ~ normal(
        0, log_lambda_gp_02_level_intercept_sd_sd[lv]
      );
    }
  }
  for (lv in 1:n_levels) {
    int mode = enable_ms_level_baseline_hazard[lv];
    if (!mode) continue;
    if (mode == LEVEL_MODE_FE || mode == LEVEL_MODE_RE || mode == LEVEL_MODE_RE_GP) {
      int r_lo = raw_level_pos_ms_baseline_shared[lv];
      int r_hi = raw_level_pos_ms_baseline_shared[lv + 1] - 1;
      if (r_hi >= r_lo) {
        if (enable_student_t_hierarchy)
          raw_log_lambda_gp_02_level_intercept[r_lo:r_hi]
            ~ student_t(ms_nu_baseline_level[lv], 0, 1);
        else
          raw_log_lambda_gp_02_level_intercept[r_lo:r_hi] ~ std_normal();
      }
    }
    if (mode == LEVEL_MODE_RE_CP) {
      int c_lo = cp_level_pos_ms_baseline_shared[lv];
      int c_hi = cp_level_pos_ms_baseline_shared[lv + 1] - 1;
      if (c_hi >= c_lo) {
        if (enable_student_t_hierarchy)
          cp_log_lambda_gp_02_level_intercept[c_lo:c_hi]
            ~ student_t(ms_nu_baseline_level[lv], 0, log_lambda_gp_02_level_intercept_sd[lv]);
        else
          cp_log_lambda_gp_02_level_intercept[c_lo:c_hi]
            ~ normal(0, log_lambda_gp_02_level_intercept_sd[lv]);
      }
    }
  }
  if (n_gp_groups_ms_baseline_02 > 0) {
    to_vector(log_lambda_gp_02_level_eta) ~ std_normal();
  }

  // Time-varying covariate coefficients
  // 0->2 time-varying covariate prior
  if (enable_ms_pop_time_varying_cov && enable_ms_02_time_varying_cov && n_time_varying_covar > 0) {
    time_varying_coef_02 ~ normal(time_varying_coef_02_mean, time_varying_coef_02_sd);
  }

  // Time-invariant covariate coefficients
  if (enable_ms_pop_time_invariant_cov && n_time_invariant_covar > 0) {
    time_invariant_coef_qr_02 ~ normal(time_invariant_coef_02_mean, time_invariant_coef_02_sd);
  }

  // Multi-level random slopes
  for (lv in 1:n_levels) {
    if (n_time_invariant_covar > 0) {
      sd_level_slope_02[lv] ~ normal(0, sd_level_slope_02_sd[lv]');
    }
    if (enable_ms_level_cov[lv] && n_time_invariant_covar > 0) {
      int mode = enable_ms_level_baseline_hazard[lv];
      if (mode == LEVEL_MODE_FE || mode == LEVEL_MODE_RE || mode == LEVEL_MODE_RE_GP) {
        int r_lo = raw_level_pos_ms_slope_shared[lv];
        int r_hi = raw_level_pos_ms_slope_shared[lv + 1] - 1;
        if (r_hi >= r_lo) {
          if (enable_student_t_hierarchy)
            to_vector(raw_level_slope_02[r_lo:r_hi, :]) ~ student_t(ms_nu_slope_level[lv], 0, 1);
          else
            to_vector(raw_level_slope_02[r_lo:r_hi, :]) ~ std_normal();
        }
      }
      if (mode == LEVEL_MODE_RE_CP) {
        int c_lo = cp_level_pos_ms_slope_shared[lv];
        int c_hi = cp_level_pos_ms_slope_shared[lv + 1] - 1;
        if (c_hi >= c_lo) {
          for (k in 1:n_time_invariant_covar) {
            if (enable_student_t_hierarchy)
              cp_level_slope_02[c_lo:c_hi, k]
                ~ student_t(ms_nu_slope_level[lv], 0, sd_level_slope_02[lv, k]);
            else
              cp_level_slope_02[c_lo:c_hi, k]
                ~ normal(0, sd_level_slope_02[lv, k]);
          }
        }
      }
    }
  }
}

// ============================================================================
// 1→2 TRANSITION PRIORS (Sojourn Time GP)
// ============================================================================

if (need_12_s_gp) {
  log_lambda_gp_12_s_pop_alpha[1] ~ inv_gamma(
    log_lambda_gp_12_s_pop_alpha_alpha, log_lambda_gp_12_s_pop_alpha_beta
  );
  log_lambda_gp_12_s_pop_rho[1] ~ inv_gamma(
    log_lambda_gp_12_s_pop_rho_alpha, log_lambda_gp_12_s_pop_rho_beta
  );
  log_lambda_gp_12_s_pop_intercept[1] ~ normal(
    log_lambda_gp_12_s_pop_intercept_mean, log_lambda_gp_12_s_pop_intercept_sd
  );
  to_vector(log_lambda_gp_12_s_pop_eta) ~ std_normal();

  for (lv in 1:n_levels) {
    log_lambda_gp_12_s_level_alpha[lv] ~ inv_gamma(
      log_lambda_gp_12_s_level_alpha_alpha[lv], log_lambda_gp_12_s_level_alpha_beta[lv]
    );
    log_lambda_gp_12_s_level_rho[lv] ~ inv_gamma(
      log_lambda_gp_12_s_level_rho_alpha[lv], log_lambda_gp_12_s_level_rho_beta[lv]
    );
    if (any_re_level) {
      log_lambda_gp_12_s_level_intercept_sd[lv] ~ normal(
        0, log_lambda_gp_12_s_level_intercept_sd_sd[lv]
      );
    }
  }
  for (lv in 1:n_levels) {
    int mode = enable_ms_level_baseline_hazard[lv];
    if (!mode) continue;
    if (mode == LEVEL_MODE_FE || mode == LEVEL_MODE_RE || mode == LEVEL_MODE_RE_GP) {
      int r_lo = raw_level_pos_ms_baseline_shared[lv];
      int r_hi = raw_level_pos_ms_baseline_shared[lv + 1] - 1;
      if (r_hi >= r_lo) {
        if (enable_student_t_hierarchy)
          raw_log_lambda_gp_12_s_level_intercept[r_lo:r_hi]
            ~ student_t(ms_nu_baseline_level[lv], 0, 1);
        else
          raw_log_lambda_gp_12_s_level_intercept[r_lo:r_hi] ~ std_normal();
      }
    }
    if (mode == LEVEL_MODE_RE_CP) {
      int c_lo = cp_level_pos_ms_baseline_shared[lv];
      int c_hi = cp_level_pos_ms_baseline_shared[lv + 1] - 1;
      if (c_hi >= c_lo) {
        if (enable_student_t_hierarchy)
          cp_log_lambda_gp_12_s_level_intercept[c_lo:c_hi]
            ~ student_t(ms_nu_baseline_level[lv], 0, log_lambda_gp_12_s_level_intercept_sd[lv]);
        else
          cp_log_lambda_gp_12_s_level_intercept[c_lo:c_hi]
            ~ normal(0, log_lambda_gp_12_s_level_intercept_sd[lv]);
      }
    }
  }
  if (n_gp_groups_ms_baseline_12_s > 0) {
    to_vector(log_lambda_gp_12_s_level_eta) ~ std_normal();
  }
}

// ============================================================================
// 1→2 TRANSITION PRIORS (Clock-forward Time GP)
// ============================================================================

if (need_12_t_gp) {
  if (!share_dead_gp_shape) {
    log_lambda_gp_12_t_pop_alpha[1] ~ inv_gamma(
      log_lambda_gp_12_t_pop_alpha_alpha, log_lambda_gp_12_t_pop_alpha_beta
    );
    log_lambda_gp_12_t_pop_rho[1] ~ inv_gamma(
      log_lambda_gp_12_t_pop_rho_alpha, log_lambda_gp_12_t_pop_rho_beta
    );
    to_vector(log_lambda_gp_12_t_pop_eta) ~ std_normal();
  }
  log_lambda_gp_12_t_pop_intercept[1] ~ normal(
    log_lambda_gp_12_t_pop_intercept_mean, log_lambda_gp_12_t_pop_intercept_sd
  );

  for (lv in 1:n_levels) {
    log_lambda_gp_12_t_level_alpha[lv] ~ inv_gamma(
      log_lambda_gp_12_t_level_alpha_alpha[lv], log_lambda_gp_12_t_level_alpha_beta[lv]
    );
    log_lambda_gp_12_t_level_rho[lv] ~ inv_gamma(
      log_lambda_gp_12_t_level_rho_alpha[lv], log_lambda_gp_12_t_level_rho_beta[lv]
    );
    if (any_re_level) {
      log_lambda_gp_12_t_level_intercept_sd[lv] ~ normal(
        0, log_lambda_gp_12_t_level_intercept_sd_sd[lv]
      );
    }
  }
  for (lv in 1:n_levels) {
    int mode = enable_ms_level_baseline_hazard[lv];
    if (!mode) continue;
    if (mode == LEVEL_MODE_FE || mode == LEVEL_MODE_RE || mode == LEVEL_MODE_RE_GP) {
      int r_lo = raw_level_pos_ms_baseline_shared[lv];
      int r_hi = raw_level_pos_ms_baseline_shared[lv + 1] - 1;
      if (r_hi >= r_lo) {
        if (enable_student_t_hierarchy)
          raw_log_lambda_gp_12_t_level_intercept[r_lo:r_hi]
            ~ student_t(ms_nu_baseline_level[lv], 0, 1);
        else
          raw_log_lambda_gp_12_t_level_intercept[r_lo:r_hi] ~ std_normal();
      }
    }
    if (mode == LEVEL_MODE_RE_CP) {
      int c_lo = cp_level_pos_ms_baseline_shared[lv];
      int c_hi = cp_level_pos_ms_baseline_shared[lv + 1] - 1;
      if (c_hi >= c_lo) {
        if (enable_student_t_hierarchy)
          cp_log_lambda_gp_12_t_level_intercept[c_lo:c_hi]
            ~ student_t(ms_nu_baseline_level[lv], 0, log_lambda_gp_12_t_level_intercept_sd[lv]);
        else
          cp_log_lambda_gp_12_t_level_intercept[c_lo:c_hi]
            ~ normal(0, log_lambda_gp_12_t_level_intercept_sd[lv]);
      }
    }
  }
  if (n_gp_groups_ms_baseline_12_t > 0) {
    to_vector(log_lambda_gp_12_t_level_eta) ~ std_normal();
  }
}

// ============================================================================
// SHARED DEAD GP SHAPE PRIORS (when share_dead_gp_shape=1)
// ============================================================================
if (share_dead_gp_shape) {
  log_lambda_gp_dead_pop_alpha[1] ~ inv_gamma(
    log_lambda_gp_dead_pop_alpha_alpha, log_lambda_gp_dead_pop_alpha_beta
  );
  log_lambda_gp_dead_pop_rho[1] ~ inv_gamma(
    log_lambda_gp_dead_pop_rho_alpha, log_lambda_gp_dead_pop_rho_beta
  );
  to_vector(log_lambda_gp_dead_pop_eta) ~ std_normal();
}

// ============================================================================
// 0→3 TRANSITION PRIORS (Dropout, GP baseline hazard with N-level hierarchy)
// ============================================================================
if (enable_ms_03) {
  log_lambda_gp_03_pop_alpha[1] ~ inv_gamma(
    log_lambda_gp_03_pop_alpha_alpha, log_lambda_gp_03_pop_alpha_beta
  );
  log_lambda_gp_03_pop_rho[1] ~ inv_gamma(
    log_lambda_gp_03_pop_rho_alpha, log_lambda_gp_03_pop_rho_beta
  );
  log_lambda_gp_03_pop_intercept[1] ~ normal(
    log_lambda_gp_03_pop_intercept_mean, log_lambda_gp_03_pop_intercept_sd
  );
  to_vector(log_lambda_gp_03_pop_eta) ~ std_normal();

  for (lv in 1:n_levels) {
    log_lambda_gp_03_level_alpha[lv] ~ inv_gamma(
      log_lambda_gp_03_level_alpha_alpha[lv], log_lambda_gp_03_level_alpha_beta[lv]
    );
    log_lambda_gp_03_level_rho[lv] ~ inv_gamma(
      log_lambda_gp_03_level_rho_alpha[lv], log_lambda_gp_03_level_rho_beta[lv]
    );
    if (any_re_level) {
      log_lambda_gp_03_level_intercept_sd[lv] ~ normal(
        0, log_lambda_gp_03_level_intercept_sd_sd[lv]
      );
    }
  }

  for (lv in 1:n_levels) {
    int mode = enable_ms_level_baseline_hazard[lv];
    if (!mode) continue;
    if (mode == LEVEL_MODE_FE || mode == LEVEL_MODE_RE || mode == LEVEL_MODE_RE_GP) {
      int r_lo = raw_level_pos_ms_baseline_shared[lv];
      int r_hi = raw_level_pos_ms_baseline_shared[lv + 1] - 1;
      if (r_hi >= r_lo) {
        if (enable_student_t_hierarchy)
          raw_log_lambda_gp_03_level_intercept[r_lo:r_hi]
            ~ student_t(ms_nu_baseline_level[lv], 0, 1);
        else
          raw_log_lambda_gp_03_level_intercept[r_lo:r_hi] ~ std_normal();
      }
    }
    if (mode == LEVEL_MODE_RE_CP) {
      int c_lo = cp_level_pos_ms_baseline_shared[lv];
      int c_hi = cp_level_pos_ms_baseline_shared[lv + 1] - 1;
      if (c_hi >= c_lo) {
        if (enable_student_t_hierarchy)
          cp_log_lambda_gp_03_level_intercept[c_lo:c_hi]
            ~ student_t(ms_nu_baseline_level[lv], 0, log_lambda_gp_03_level_intercept_sd[lv]);
        else
          cp_log_lambda_gp_03_level_intercept[c_lo:c_hi]
            ~ normal(0, log_lambda_gp_03_level_intercept_sd[lv]);
      }
    }
  }
  if (n_gp_groups_ms_baseline_03 > 0) {
    to_vector(log_lambda_gp_03_level_eta) ~ std_normal();
  }
}

// ============================================================================
// 3→2 TRANSITION PRIORS (Off-trial death, sojourn time GP baseline hazard)
// ============================================================================
if (enable_ms_32) {
  // Population-level baseline hazard GP
  log_lambda_gp_32_s_pop_alpha[1] ~ inv_gamma(
    log_lambda_gp_32_s_pop_alpha_alpha, log_lambda_gp_32_s_pop_alpha_beta
  );
  log_lambda_gp_32_s_pop_rho[1] ~ inv_gamma(
    log_lambda_gp_32_s_pop_rho_alpha, log_lambda_gp_32_s_pop_rho_beta
  );
  log_lambda_gp_32_s_pop_intercept[1] ~ normal(
    log_lambda_gp_32_s_pop_intercept_mean, log_lambda_gp_32_s_pop_intercept_sd
  );
  to_vector(log_lambda_gp_32_s_pop_eta) ~ std_normal();

  // Level-level baseline hazard GP
  for (lv in 1:n_levels) {
    log_lambda_gp_32_s_level_alpha[lv] ~ inv_gamma(
      log_lambda_gp_32_s_level_alpha_alpha[lv], log_lambda_gp_32_s_level_alpha_beta[lv]
    );
    log_lambda_gp_32_s_level_rho[lv] ~ inv_gamma(
      log_lambda_gp_32_s_level_rho_alpha[lv], log_lambda_gp_32_s_level_rho_beta[lv]
    );
    if (any_re_level) {
      log_lambda_gp_32_s_level_intercept_sd[lv] ~ normal(
        0, log_lambda_gp_32_s_level_intercept_sd_sd[lv]
      );
    }
  }
  for (lv in 1:n_levels) {
    int mode = enable_ms_level_baseline_hazard[lv];
    if (!mode) continue;
    if (mode == LEVEL_MODE_FE || mode == LEVEL_MODE_RE || mode == LEVEL_MODE_RE_GP) {
      int r_lo = raw_level_pos_ms_baseline_shared[lv];
      int r_hi = raw_level_pos_ms_baseline_shared[lv + 1] - 1;
      if (r_hi >= r_lo) {
        if (enable_student_t_hierarchy)
          raw_log_lambda_gp_32_s_level_intercept[r_lo:r_hi]
            ~ student_t(ms_nu_baseline_level[lv], 0, 1);
        else
          raw_log_lambda_gp_32_s_level_intercept[r_lo:r_hi] ~ std_normal();
      }
    }
    if (mode == LEVEL_MODE_RE_CP) {
      int c_lo = cp_level_pos_ms_baseline_shared[lv];
      int c_hi = cp_level_pos_ms_baseline_shared[lv + 1] - 1;
      if (c_hi >= c_lo) {
        if (enable_student_t_hierarchy)
          cp_log_lambda_gp_32_s_level_intercept[c_lo:c_hi]
            ~ student_t(ms_nu_baseline_level[lv], 0, log_lambda_gp_32_s_level_intercept_sd[lv]);
        else
          cp_log_lambda_gp_32_s_level_intercept[c_lo:c_hi]
            ~ normal(0, log_lambda_gp_32_s_level_intercept_sd[lv]);
      }
    }
  }
  if (n_gp_groups_ms_baseline_32 > 0) {
    to_vector(log_lambda_gp_32_s_level_eta) ~ std_normal();
  }
}

// Student-t nu priors for multistate (only when enabled)
if (enable_student_t_hierarchy) {
  for (lv in 1:n_levels) {
    ms_nu_baseline_level[lv] ~ gamma(ms_nu_baseline_level_prior_alpha[lv],
                                      ms_nu_baseline_level_prior_beta[lv]);
    ms_nu_slope_level[lv] ~ gamma(ms_nu_slope_level_prior_alpha[lv],
                                   ms_nu_slope_level_prior_beta[lv]);
  }
}

// PSA-at-entry covariate priors
if (enable_ms_12 && enable_ms_12_entry_covar) {
  coef_log_entry_covar_12[1] ~ normal(coef_log_entry_covar_12_mean, coef_log_entry_covar_12_sd);
}
if (enable_ms_32 && enable_ms_32_entry_covar) {
  coef_log_entry_covar_32[1] ~ normal(coef_log_entry_covar_32_mean, coef_log_entry_covar_32_sd);
}

// 1→2 covariate priors
if (enable_ms_12) {
  if (enable_ms_pop_time_varying_cov && n_time_varying_covar > 0) {
    time_varying_coef_12 ~ normal(time_varying_coef_12_mean, time_varying_coef_12_sd);
  }

  if (enable_ms_pop_time_invariant_cov && n_time_invariant_covar > 0) {
    time_invariant_coef_qr_12 ~ normal(time_invariant_coef_12_mean, time_invariant_coef_12_sd);
  }

  for (lv in 1:n_levels) {
    if (n_time_invariant_covar > 0) {
      sd_level_slope_12[lv] ~ normal(0, sd_level_slope_12_sd[lv]');
    }
    if (enable_ms_level_cov[lv] && n_time_invariant_covar > 0) {
      int mode = enable_ms_level_baseline_hazard[lv];
      if (mode == LEVEL_MODE_FE || mode == LEVEL_MODE_RE || mode == LEVEL_MODE_RE_GP) {
        int r_lo = raw_level_pos_ms_slope_shared[lv];
        int r_hi = raw_level_pos_ms_slope_shared[lv + 1] - 1;
        if (r_hi >= r_lo) {
          if (enable_student_t_hierarchy)
            to_vector(raw_level_slope_12[r_lo:r_hi, :]) ~ student_t(ms_nu_slope_level[lv], 0, 1);
          else
            to_vector(raw_level_slope_12[r_lo:r_hi, :]) ~ std_normal();
        }
      }
      if (mode == LEVEL_MODE_RE_CP) {
        int c_lo = cp_level_pos_ms_slope_shared[lv];
        int c_hi = cp_level_pos_ms_slope_shared[lv + 1] - 1;
        if (c_hi >= c_lo) {
          for (k in 1:n_time_invariant_covar) {
            if (enable_student_t_hierarchy)
              cp_level_slope_12[c_lo:c_hi, k]
                ~ student_t(ms_nu_slope_level[lv], 0, sd_level_slope_12[lv, k]);
            else
              cp_level_slope_12[c_lo:c_hi, k]
                ~ normal(0, sd_level_slope_12[lv, k]);
          }
        }
      }
    }
  }
}
