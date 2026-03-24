// ============================================================================
// Multistate Hazard Model Transformed Parameters
// ============================================================================

// ============================================================================
// 0→1 TRANSITION: Log Conditional Survival
// ============================================================================

// Population-level baseline hazard
row_vector[enable_ms_01 ? max_all_t : 0] log_pop_lambda_01;

// Level-level baseline hazard (additive to population)
matrix[n_enabled_groups_ms_baseline_01, enable_ms_01 ? max_all_t : 0] log_level_lambda_01_residual;
vector[n_enabled_groups_ms_baseline_01] log_lambda_gp_01_level_intercept;

// Patient-level log conditional survival probability
matrix[enable_ms_01 ? n_hmc_patients : 0, enable_ms_01 ? max_all_t : 0] log_cond_surv_01;

if (enable_ms_01) {
  // Compute population GP on coarse grid then expand to weekly
  log_pop_lambda_01 = calc_gp_pred(
    ms_gp_cal_t,
    log_lambda_gp_01_pop_intercept[1],
    log_lambda_gp_01_pop_alpha[1],
    log_lambda_gp_01_pop_rho[1],
    delta,
    log_lambda_gp_01_pop_eta
  )[knot_of_cal];

  // Initialize patient hazards with population baseline
  log_cond_surv_01 = rep_matrix(log_pop_lambda_01, n_hmc_patients);

  // Add level-level residuals (intercept-only or full GP)
  for (lv in 1:n_levels) {
    if (enable_ms_level_baseline_hazard[lv]) {
      int lv_start, lv_end;
      (lv_start, lv_end) = get_pos(enabled_level_pos_ms_baseline, lv);

      // Scale intercepts (shared by both intercept-only and GP modes)
      log_lambda_gp_01_level_intercept[lv_start:lv_end] =
        raw_log_lambda_gp_01_level_intercept[lv_start:lv_end] *
        log_lambda_gp_01_level_intercept_sd[lv];

      if (enable_ms_level_baseline_hazard[lv] == 2) {
        // GP mode: full time-varying residual via calc_gp_pred
        int gp_start, gp_end;
        (gp_start, gp_end) = get_pos(gp_level_pos_ms_baseline, lv);
        for (g in lv_start:lv_end) {
          int g_gp = gp_start + (g - lv_start);
          log_level_lambda_01_residual[g] = calc_gp_pred(
            ms_gp_cal_t,
            log_lambda_gp_01_level_intercept[g],
            log_lambda_gp_01_level_alpha[lv],
            log_lambda_gp_01_level_rho[lv],
            delta,
            log_lambda_gp_01_level_eta[g_gp]
          )[knot_of_cal];
        }
      } else {
        // Intercept-only mode: constant shift across all time points
        for (g in lv_start:lv_end) {
          log_level_lambda_01_residual[g] = rep_row_vector(
            log_lambda_gp_01_level_intercept[g], max_all_t
          );
        }
      }

      // Add level residuals to patient hazards
      for (j in 1:n_hmc_patients) {
        int p = hmc_patient_idx[j];
        log_cond_surv_01[j] += log_level_lambda_01_residual[patient_ms_baseline_flat_idx[p, lv]];
      }
    }
  }

  // -------------------------------------------------------------------------
  // Time-varying covariate effects (generic - no feature knowledge)
  // -------------------------------------------------------------------------
  // 0->1 time-varying covariate: add to log_cond_surv_01 ONLY in continuous mode.
  // In visit-gated mode, the covariate is applied per-visit in multistate_lpmf.
  if (enable_ms_pop_time_varying_cov && n_time_varying_covar > 0 && !enable_ms_visit_gated_01) {
    for (k in 1:n_time_varying_covar) {
      log_cond_surv_01 += time_varying_coef_01[k] * ms_time_varying_covar_01[k];
    }
  }

  // -------------------------------------------------------------------------
  // Time-invariant covariate effects (QR space with multi-level random slopes)
  // -------------------------------------------------------------------------
  if (enable_ms_pop_time_invariant_cov && n_time_invariant_covar > 0) {
    // Population-level covariate effects (QR space)
    vector[n_hmc_patients] linpred_pop_01 = Q_covar_design_matrix[hmc_patient_idx, :] * time_invariant_coef_qr_01;

    // Multi-level random slopes (if enabled)
    if (n_enabled_groups_ms_slope > 0) {
      // Step 1: Scale all raw slope effects at once
      matrix[n_enabled_groups_ms_slope, n_time_invariant_covar] ms_scaled_level_slope_01;
      for (lv in 1:n_levels) {
        if (enable_ms_level_cov[lv]) {
          int lv_start, lv_end;
          (lv_start, lv_end) = get_pos(enabled_level_pos_ms_slope, lv);
          ms_scaled_level_slope_01[lv_start:lv_end, :] =
            raw_level_slope_01[lv_start:lv_end, :] .*
            rep_matrix(sd_level_slope_01[lv]', lv_end - lv_start + 1);
        }
      }

      // Step 2: Gather and compute dot products using pre-computed flat indices
      for (lv in 1:n_levels) {
        if (enable_ms_level_cov[lv]) {
          linpred_pop_01 += rows_dot_product(
            Q_covar_design_matrix[hmc_patient_idx, :],
            ms_scaled_level_slope_01[patient_ms_slope_flat_idx[hmc_patient_idx, lv], :]
          );
        }
      }
    }

    // Add time-invariant effects (broadcast to all times)
    log_cond_surv_01 += rep_matrix(linpred_pop_01, max_all_t);
  }

  // Transform log-hazard to log conditional survival probability
  // log P(survive interval t) = -exp(log_hazard[t]) = -hazard[t]
  // Clamp log-hazard to [-20, 10] before exp():
  //   upper 10: exp(10)≈22000/week already means instantaneous death; tighter than
  //             the old 35 cap (exp(35)≈1.5e15/week) which created lp≈-1e15 with
  //             zero gradient (fmin is piecewise-constant) trapping HMC. Consistent
  //             with the laplace.stanfunctions cap of 5–10 for non-target patients.
  //   lower -20: prevents exp() underflow to 0 → log_cond_surv = 0 →
  //              log1m_exp(0) = -Inf at event times.
  log_cond_surv_01 = -exp(fmax(log_cond_surv_01, -20.0));
}

// ============================================================================
// 0→2 TRANSITION: Log Conditional Survival
// ============================================================================

row_vector[enable_ms_02 ? max_all_t : 0] log_pop_lambda_02;
matrix[n_enabled_groups_ms_baseline_02, enable_ms_02 ? max_all_t : 0] log_level_lambda_02_residual;
vector[n_enabled_groups_ms_baseline_02] log_lambda_gp_02_level_intercept;
matrix[enable_ms_02 ? n_hmc_patients : 0, enable_ms_02 ? max_all_t : 0] log_cond_surv_02;

if (enable_ms_02) {
  // Compute population GP on coarse grid then expand to weekly
  log_pop_lambda_02 = calc_gp_pred(
    ms_gp_cal_t,
    log_lambda_gp_02_pop_intercept[1],
    log_lambda_gp_02_pop_alpha[1],
    log_lambda_gp_02_pop_rho[1],
    delta,
    log_lambda_gp_02_pop_eta
  )[knot_of_cal];

  // Initialize with population baseline
  log_cond_surv_02 = rep_matrix(log_pop_lambda_02, n_hmc_patients);

  // Add level-level residuals (intercept-only or full GP, same pattern as 0→1)
  for (lv in 1:n_levels) {
    if (enable_ms_level_baseline_hazard[lv]) {
      int lv_start, lv_end;
      (lv_start, lv_end) = get_pos(enabled_level_pos_ms_baseline, lv);

      log_lambda_gp_02_level_intercept[lv_start:lv_end] =
        raw_log_lambda_gp_02_level_intercept[lv_start:lv_end] *
        log_lambda_gp_02_level_intercept_sd[lv];

      if (enable_ms_level_baseline_hazard[lv] == 2) {
        int gp_start, gp_end;
        (gp_start, gp_end) = get_pos(gp_level_pos_ms_baseline, lv);
        for (g in lv_start:lv_end) {
          int g_gp = gp_start + (g - lv_start);
          log_level_lambda_02_residual[g] = calc_gp_pred(
            ms_gp_cal_t,
            log_lambda_gp_02_level_intercept[g],
            log_lambda_gp_02_level_alpha[lv],
            log_lambda_gp_02_level_rho[lv],
            delta,
            log_lambda_gp_02_level_eta[g_gp]
          )[knot_of_cal];
        }
      } else {
        for (g in lv_start:lv_end) {
          log_level_lambda_02_residual[g] = rep_row_vector(
            log_lambda_gp_02_level_intercept[g], max_all_t
          );
        }
      }

      for (j in 1:n_hmc_patients) {
        int p = hmc_patient_idx[j];
        log_cond_surv_02[j] += log_level_lambda_02_residual[patient_ms_baseline_flat_idx[p, lv]];
      }
    }
  }

  // -------------------------------------------------------------------------
  // Time-varying covariate effects (same pattern as 0→1)
  // -------------------------------------------------------------------------
  // 0->2 time-varying covariate: uses modeled PSA (ms_time_varying_covar_01)
  // Controlled by enable_ms_02_time_varying_cov, independent of 0->1 mode.
  if (enable_ms_pop_time_varying_cov && enable_ms_02_time_varying_cov && n_time_varying_covar > 0) {
    for (k in 1:n_time_varying_covar) {
      log_cond_surv_02 += time_varying_coef_02[k] * ms_time_varying_covar_01[k];
    }
  }

  // -------------------------------------------------------------------------
  // Time-invariant covariate effects
  // -------------------------------------------------------------------------
  if (enable_ms_pop_time_invariant_cov && n_time_invariant_covar > 0) {
    vector[n_hmc_patients] linpred_pop_02 = Q_covar_design_matrix[hmc_patient_idx, :] * time_invariant_coef_qr_02;

    if (n_enabled_groups_ms_slope > 0) {
      matrix[n_enabled_groups_ms_slope, n_time_invariant_covar] ms_scaled_level_slope_02;
      for (lv in 1:n_levels) {
        if (enable_ms_level_cov[lv]) {
          int lv_start, lv_end;
          (lv_start, lv_end) = get_pos(enabled_level_pos_ms_slope, lv);
          ms_scaled_level_slope_02[lv_start:lv_end, :] =
            raw_level_slope_02[lv_start:lv_end, :] .*
            rep_matrix(sd_level_slope_02[lv]', lv_end - lv_start + 1);
        }
      }

      for (lv in 1:n_levels) {
        if (enable_ms_level_cov[lv]) {
          linpred_pop_02 += rows_dot_product(
            Q_covar_design_matrix[hmc_patient_idx, :],
            ms_scaled_level_slope_02[patient_ms_slope_flat_idx[hmc_patient_idx, lv], :]
          );
        }
      }
    }

    log_cond_surv_02 += rep_matrix(linpred_pop_02, max_all_t);
  }

  // Transform to log conditional survival (clamped, see 0→1 comment)
  log_cond_surv_02 = -exp(fmax(log_cond_surv_02, -20.0));
}

// ============================================================================
// 1→2 TRANSITION: Log Conditional Survival (on sojourn and/or clock-forward time)
// ============================================================================

// Sojourn time GP (semi-Markov or extended)
row_vector[need_12_s_gp ? ms_max_sojourn_t : 0] log_pop_lambda_12_s;
matrix[n_enabled_groups_ms_baseline_12_s, need_12_s_gp ? ms_max_sojourn_t : 0] log_level_lambda_12_s_residual;
vector[n_enabled_groups_ms_baseline_12_s] log_lambda_gp_12_s_level_intercept;
matrix[need_12_s_gp ? n_hmc_patients : 0, need_12_s_gp ? ms_max_sojourn_t : 0] log_cond_surv_12_s;

// Clock-forward time GP (Markov or extended)
row_vector[need_12_t_gp ? max_all_t : 0] log_pop_lambda_12_t;
matrix[n_enabled_groups_ms_baseline_12_t, need_12_t_gp ? max_all_t : 0] log_level_lambda_12_t_residual;
vector[n_enabled_groups_ms_baseline_12_t] log_lambda_gp_12_t_level_intercept;
matrix[need_12_t_gp ? n_hmc_patients : 0, need_12_t_gp ? max_all_t : 0] log_cond_surv_12_t;

if (need_12_s_gp) {
  // Sojourn time GP on coarse grid then expand to weekly
  log_pop_lambda_12_s = calc_gp_pred(
    ms_gp_sojourn_t,
    log_lambda_gp_12_s_pop_intercept[1],
    log_lambda_gp_12_s_pop_alpha[1],
    log_lambda_gp_12_s_pop_rho[1],
    delta,
    log_lambda_gp_12_s_pop_eta
  )[knot_of_sojourn];

  log_cond_surv_12_s = rep_matrix(log_pop_lambda_12_s, n_hmc_patients);

  // Add level residuals (intercept-only or full GP)
  for (lv in 1:n_levels) {
    if (enable_ms_level_baseline_hazard[lv]) {
      int lv_start, lv_end;
      (lv_start, lv_end) = get_pos(enabled_level_pos_ms_baseline, lv);

      log_lambda_gp_12_s_level_intercept[lv_start:lv_end] =
        raw_log_lambda_gp_12_s_level_intercept[lv_start:lv_end] *
        log_lambda_gp_12_s_level_intercept_sd[lv];

      if (enable_ms_level_baseline_hazard[lv] == 2) {
        int gp_start, gp_end;
        (gp_start, gp_end) = get_pos(gp_level_pos_ms_baseline, lv);
        for (g in lv_start:lv_end) {
          int g_gp = gp_start + (g - lv_start);
          log_level_lambda_12_s_residual[g] = calc_gp_pred(
            ms_gp_sojourn_t,
            log_lambda_gp_12_s_level_intercept[g],
            log_lambda_gp_12_s_level_alpha[lv],
            log_lambda_gp_12_s_level_rho[lv],
            delta,
            log_lambda_gp_12_s_level_eta[g_gp]
          )[knot_of_sojourn];
        }
      } else {
        for (g in lv_start:lv_end) {
          log_level_lambda_12_s_residual[g] = rep_row_vector(
            log_lambda_gp_12_s_level_intercept[g], ms_max_sojourn_t  // output is weekly-size
          );
        }
      }

      for (j in 1:n_hmc_patients) {
        int p = hmc_patient_idx[j];
        log_cond_surv_12_s[j] += log_level_lambda_12_s_residual[patient_ms_baseline_flat_idx[p, lv]];
      }
    }
  }

  // -------------------------------------------------------------------------
  // Time-invariant covariate effects for 1→2 sojourn
  // -------------------------------------------------------------------------
  // Note: Time-varying covariates for 1→2 would need different indexing (sojourn time)
  // For now, only time-invariant covariates are supported for 1→2
  if (enable_ms_pop_time_invariant_cov && n_time_invariant_covar > 0) {
    vector[n_hmc_patients] linpred_pop_12 = Q_covar_design_matrix[hmc_patient_idx, :] * time_invariant_coef_qr_12;

    if (n_enabled_groups_ms_slope > 0) {
      matrix[n_enabled_groups_ms_slope, n_time_invariant_covar] ms_scaled_level_slope_12;
      for (lv in 1:n_levels) {
        if (enable_ms_level_cov[lv]) {
          int lv_start, lv_end;
          (lv_start, lv_end) = get_pos(enabled_level_pos_ms_slope, lv);
          ms_scaled_level_slope_12[lv_start:lv_end, :] =
            raw_level_slope_12[lv_start:lv_end, :] .*
            rep_matrix(sd_level_slope_12[lv]', lv_end - lv_start + 1);
        }
      }

      for (lv in 1:n_levels) {
        if (enable_ms_level_cov[lv]) {
          linpred_pop_12 += rows_dot_product(
            Q_covar_design_matrix[hmc_patient_idx, :],
            ms_scaled_level_slope_12[patient_ms_slope_flat_idx[hmc_patient_idx, lv], :]
          );
        }
      }
    }

    log_cond_surv_12_s += rep_matrix(linpred_pop_12, ms_max_sojourn_t);
  }

  log_cond_surv_12_s = -exp(fmax(log_cond_surv_12_s, -20.0));
}

if (need_12_t_gp) {
  // Clock-forward time GP on coarse grid then expand to weekly
  // In extended mode, intercept is zero (sojourn GP carries it) to avoid non-identifiability
  log_pop_lambda_12_t = calc_gp_pred(
    ms_gp_cal_t,
    ms_12_t_has_intercept ? log_lambda_gp_12_t_pop_intercept[1] : 0.0,
    log_lambda_gp_12_t_pop_alpha[1],
    log_lambda_gp_12_t_pop_rho[1],
    delta,
    log_lambda_gp_12_t_pop_eta
  )[knot_of_cal];

  log_cond_surv_12_t = rep_matrix(log_pop_lambda_12_t, n_hmc_patients);

  // Add level residuals (intercept-only or full GP)
  for (lv in 1:n_levels) {
    if (enable_ms_level_baseline_hazard[lv]) {
      int lv_start, lv_end;
      (lv_start, lv_end) = get_pos(enabled_level_pos_ms_baseline, lv);

      log_lambda_gp_12_t_level_intercept[lv_start:lv_end] =
        raw_log_lambda_gp_12_t_level_intercept[lv_start:lv_end] *
        log_lambda_gp_12_t_level_intercept_sd[lv];

      if (enable_ms_level_baseline_hazard[lv] == 2) {
        int gp_start, gp_end;
        (gp_start, gp_end) = get_pos(gp_level_pos_ms_baseline, lv);
        for (g in lv_start:lv_end) {
          int g_gp = gp_start + (g - lv_start);
          // In extended mode, level intercept is zero (sojourn GP carries it)
          log_level_lambda_12_t_residual[g] = calc_gp_pred(
            ms_gp_cal_t,
            ms_12_t_has_intercept ? log_lambda_gp_12_t_level_intercept[g] : 0.0,
            log_lambda_gp_12_t_level_alpha[lv],
            log_lambda_gp_12_t_level_rho[lv],
            delta,
            log_lambda_gp_12_t_level_eta[g_gp]
          )[knot_of_cal];
        }
      } else {
        for (g in lv_start:lv_end) {
          log_level_lambda_12_t_residual[g] = rep_row_vector(
            ms_12_t_has_intercept ? log_lambda_gp_12_t_level_intercept[g] : 0.0, max_all_t
          );
        }
      }

      for (j in 1:n_hmc_patients) {
        int p = hmc_patient_idx[j];
        log_cond_surv_12_t[j] += log_level_lambda_12_t_residual[patient_ms_baseline_flat_idx[p, lv]];
      }
    }
  }

  // -------------------------------------------------------------------------
  // Time-invariant covariate effects for 1→2 clock-forward
  // -------------------------------------------------------------------------
  if (enable_ms_pop_time_invariant_cov && n_time_invariant_covar > 0) {
    vector[n_hmc_patients] linpred_pop_12 = Q_covar_design_matrix[hmc_patient_idx, :] * time_invariant_coef_qr_12;

    if (n_enabled_groups_ms_slope > 0) {
      matrix[n_enabled_groups_ms_slope, n_time_invariant_covar] ms_scaled_level_slope_12;
      for (lv in 1:n_levels) {
        if (enable_ms_level_cov[lv]) {
          int lv_start, lv_end;
          (lv_start, lv_end) = get_pos(enabled_level_pos_ms_slope, lv);
          ms_scaled_level_slope_12[lv_start:lv_end, :] =
            raw_level_slope_12[lv_start:lv_end, :] .*
            rep_matrix(sd_level_slope_12[lv]', lv_end - lv_start + 1);
        }
      }

      for (lv in 1:n_levels) {
        if (enable_ms_level_cov[lv]) {
          linpred_pop_12 += rows_dot_product(
            Q_covar_design_matrix[hmc_patient_idx, :],
            ms_scaled_level_slope_12[patient_ms_slope_flat_idx[hmc_patient_idx, lv], :]
          );
        }
      }
    }

    log_cond_surv_12_t += rep_matrix(linpred_pop_12, max_all_t);
  }

  log_cond_surv_12_t = -exp(fmax(log_cond_surv_12_t, -20.0));
}

// ============================================================================
// 0→3 TRANSITION: GP baseline hazard (clock-forward time) with N-level hierarchy
// ============================================================================
row_vector[enable_ms_03 ? max_all_t : 0] log_pop_lambda_03;
matrix[n_enabled_groups_ms_baseline_03, enable_ms_03 ? max_all_t : 0] log_level_lambda_03_residual;
vector[n_enabled_groups_ms_baseline_03] log_lambda_gp_03_level_intercept;
matrix[enable_ms_03 ? n_hmc_patients : 0, enable_ms_03 ? max_all_t : 0] log_cond_surv_03;

if (enable_ms_03) {
  // Compute population GP on coarse grid then expand to weekly
  log_pop_lambda_03 = calc_gp_pred(
    ms_gp_cal_t,
    log_lambda_gp_03_pop_intercept[1],
    log_lambda_gp_03_pop_alpha[1],
    log_lambda_gp_03_pop_rho[1],
    delta,
    log_lambda_gp_03_pop_eta
  )[knot_of_cal];

  // Initialize patient hazards with population baseline
  log_cond_surv_03 = rep_matrix(log_pop_lambda_03, n_hmc_patients);

  // Add level-level residuals (intercept-only or full GP)
  for (lv in 1:n_levels) {
    if (enable_ms_level_baseline_hazard[lv]) {
      int lv_start, lv_end;
      (lv_start, lv_end) = get_pos(enabled_level_pos_ms_baseline, lv);

      log_lambda_gp_03_level_intercept[lv_start:lv_end] =
        raw_log_lambda_gp_03_level_intercept[lv_start:lv_end] *
        log_lambda_gp_03_level_intercept_sd[lv];

      if (enable_ms_level_baseline_hazard[lv] == 2) {
        int gp_start, gp_end;
        (gp_start, gp_end) = get_pos(gp_level_pos_ms_baseline, lv);
        for (g in lv_start:lv_end) {
          int g_gp = gp_start + (g - lv_start);
          log_level_lambda_03_residual[g] = calc_gp_pred(
            ms_gp_cal_t,
            log_lambda_gp_03_level_intercept[g],
            log_lambda_gp_03_level_alpha[lv],
            log_lambda_gp_03_level_rho[lv],
            delta,
            log_lambda_gp_03_level_eta[g_gp]
          )[knot_of_cal];
        }
      } else {
        for (g in lv_start:lv_end) {
          log_level_lambda_03_residual[g] = rep_row_vector(
            log_lambda_gp_03_level_intercept[g], max_all_t
          );
        }
      }

      for (j in 1:n_hmc_patients) {
        int p = hmc_patient_idx[j];
        log_cond_surv_03[j] += log_level_lambda_03_residual[patient_ms_baseline_flat_idx[p, lv]];
      }
    }
  }

  // Transform log-hazard to log conditional survival probability (clamped)
  log_cond_surv_03 = -exp(fmax(log_cond_surv_03, -20.0));
}

// ============================================================================
// 3→2 TRANSITION: Sojourn time GP baseline hazard (semi-Markov)
// ============================================================================

row_vector[enable_ms_32 ? ms_max_sojourn_t_32 : 0] log_pop_lambda_32;
matrix[n_enabled_groups_ms_baseline_32, enable_ms_32 ? ms_max_sojourn_t_32 : 0] log_level_lambda_32_residual;
vector[n_enabled_groups_ms_baseline_32] log_lambda_gp_32_s_level_intercept;
matrix[enable_ms_32 ? n_hmc_patients : 0, enable_ms_32 ? ms_max_sojourn_t_32 : 0] log_cond_surv_32;

if (enable_ms_32) {
  // Compute population GP on coarse grid then expand to weekly
  log_pop_lambda_32 = calc_gp_pred(
    ms_gp_sojourn_32_t,
    log_lambda_gp_32_s_pop_intercept[1],
    log_lambda_gp_32_s_pop_alpha[1],
    log_lambda_gp_32_s_pop_rho[1],
    delta,
    log_lambda_gp_32_s_pop_eta
  )[knot_of_sojourn_32];

  // Initialize with population baseline
  log_cond_surv_32 = rep_matrix(log_pop_lambda_32, n_hmc_patients);

  // Add level-level residuals (intercept-only or full GP)
  for (lv in 1:n_levels) {
    if (enable_ms_level_baseline_hazard[lv]) {
      int lv_start, lv_end;
      (lv_start, lv_end) = get_pos(enabled_level_pos_ms_baseline, lv);

      log_lambda_gp_32_s_level_intercept[lv_start:lv_end] =
        raw_log_lambda_gp_32_s_level_intercept[lv_start:lv_end] *
        log_lambda_gp_32_s_level_intercept_sd[lv];

      if (enable_ms_level_baseline_hazard[lv] == 2) {
        int gp_start, gp_end;
        (gp_start, gp_end) = get_pos(gp_level_pos_ms_baseline, lv);
        for (g in lv_start:lv_end) {
          int g_gp = gp_start + (g - lv_start);
          log_level_lambda_32_residual[g] = calc_gp_pred(
            ms_gp_sojourn_32_t,
            log_lambda_gp_32_s_level_intercept[g],
            log_lambda_gp_32_s_level_alpha[lv],
            log_lambda_gp_32_s_level_rho[lv],
            delta,
            log_lambda_gp_32_s_level_eta[g_gp]
          )[knot_of_sojourn_32];
        }
      } else {
        for (g in lv_start:lv_end) {
          log_level_lambda_32_residual[g] = rep_row_vector(
            log_lambda_gp_32_s_level_intercept[g], ms_max_sojourn_t_32  // output is weekly-size
          );
        }
      }

      for (j in 1:n_hmc_patients) {
        int p = hmc_patient_idx[j];
        log_cond_surv_32[j] += log_level_lambda_32_residual[patient_ms_baseline_flat_idx[p, lv]];
      }
    }
  }

  // No covariates for 3→2 in v1

  // Transform to log conditional survival
  log_cond_surv_32 = -exp(fmax(log_cond_surv_32, -20.0));
}

