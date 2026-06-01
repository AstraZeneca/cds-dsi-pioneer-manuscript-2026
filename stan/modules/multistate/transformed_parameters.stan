// ============================================================================
// Multistate Hazard Model Transformed Parameters
// ============================================================================
// When ms_needs_inline_burden = TRUE, burden models compute time-varying covariates
// inline in _ms_burden_inline_tv_covar.stan instead of via the grid-based
// ms_time_varying_covar_01 matrix. Guards on !ms_needs_inline_burden below prevent
// double-counting.

// ============================================================================
// CORRELATED INTERCEPT BLOCKS (Phase 2 — cross-transition frailty)
// ============================================================================
// Assemble the MVN-correlated level intercepts once, up front. Row m of block b
// holds the correlated intercept for that member's transition, one column per
// GROUP at the block's level. Each member-row is scattered into its transition's
// own log_lambda_gp_<slot>_level_intercept slice further below (replacing the
// scalar sigma*raw path for correlated members only).
//   u_b = diag_pre_multiply(sigma_members, L_b) * z_b
// sigma_members[m] reuses the existing per-transition level-lv intercept SD
// (no new scale hyperparameters). When n_ms_corr_blocks == 0 this block is empty.
array[n_ms_corr_blocks] matrix[ms_corr_dim, ms_corr_n_groups] ms_corr_u;
for (b in 1:n_ms_corr_blocks) {
  int lv = ms_corr_block_level[b];
  vector[ms_corr_dim] sigma_members;
  for (m in 1:ms_corr_dim) {
    int k = ms_corr_block_member_slot[b, m];
    // Gather the member transition's level-lv intercept SD. Dispatch on slot;
    // the SD arrays are the same ones the scalar NCP path multiplies by.
    real s;
    if (k == MS_SLOT_01)        s = log_lambda_gp_01_level_intercept_sd[lv];
    else if (k == MS_SLOT_02)   s = log_lambda_gp_02_level_intercept_sd[lv];
    else if (k == MS_SLOT_03)   s = log_lambda_gp_03_level_intercept_sd[lv];
    else if (k == MS_SLOT_12_S) s = log_lambda_gp_12_s_level_intercept_sd[lv];
    else if (k == MS_SLOT_12_T) s = log_lambda_gp_12_t_level_intercept_sd[lv];
    else                        s = log_lambda_gp_32_s_level_intercept_sd[lv];
    sigma_members[m] = s;
  }
  ms_corr_u[b] = diag_pre_multiply(sigma_members, L_ms_intercept_corr[b])
                 * z_ms_intercept[b];
}

// ============================================================================
// 0→1 TRANSITION: Log Conditional Survival
// ============================================================================

// Population-level baseline hazard
row_vector[enable_ms_01 ? max_all_t : 0] log_pop_lambda_01;

// Level-level baseline hazard (additive to population)
matrix[n_enabled_groups_ms_baseline_01, enable_ms_01 ? max_all_t : 0] log_level_lambda_01_residual;
vector[n_enabled_groups_ms_baseline_01] log_lambda_gp_01_level_intercept;

// Patient-level log conditional survival probability
matrix[enable_ms_01 ? n_forecast_patients : 0, enable_ms_01 ? max_all_t : 0] log_cond_surv_01;

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
  log_cond_surv_01 = rep_matrix(log_pop_lambda_01, n_forecast_patients);

  // Add level-level residuals (intercept-only or full GP)
  for (lv in 1:n_levels) {
    if (ms_legacy_mode[MS_SLOT_01, lv]) {
      int lv_start, lv_end;
      (lv_start, lv_end) = get_pos(enabled_level_pos_ms_baseline_slot[MS_SLOT_01], lv);

      // Scale intercepts: route by mode
      int mode_01 = ms_legacy_mode[MS_SLOT_01, lv];
      if (mode_01 == LEVEL_MODE_RE_CP) {
        // CP path: centered parameters are already at natural scale
        int c_lo = cp_level_pos_ms_baseline_slot[MS_SLOT_01][lv];
        int c_hi = cp_level_pos_ms_baseline_slot[MS_SLOT_01][lv + 1] - 1;
        log_lambda_gp_01_level_intercept[lv_start:lv_end] =
          cp_log_lambda_gp_01_level_intercept[c_lo:c_hi];
      } else {
        int r_lo = raw_level_pos_ms_baseline_slot[MS_SLOT_01][lv];
        int r_hi = raw_level_pos_ms_baseline_slot[MS_SLOT_01][lv + 1] - 1;
        if (mode_01 == LEVEL_MODE_FE) {
          // Fixed effects: no pooling
          log_lambda_gp_01_level_intercept[lv_start:lv_end] =
            raw_log_lambda_gp_01_level_intercept[r_lo:r_hi] *
            fe_log_lambda_gp_01_level_intercept_sd[lv];
        } else if (ms_corr_member_block[MS_SLOT_01, lv] > 0) {
          // Correlated RE member: take this transition's MVN row (already scaled
          // by sigma in ms_corr_u). One column per group at level lv.
          int b = ms_corr_member_block[MS_SLOT_01, lv];
          int row = ms_corr_member_row[MS_SLOT_01, lv];
          log_lambda_gp_01_level_intercept[lv_start:lv_end] =
            ms_corr_u[b][row, ]';
        } else {
          // Random effects (2 or 3): hierarchical pooling
          log_lambda_gp_01_level_intercept[lv_start:lv_end] =
            raw_log_lambda_gp_01_level_intercept[r_lo:r_hi] *
            log_lambda_gp_01_level_intercept_sd[lv];
        }
      }

      if (ms_legacy_mode[MS_SLOT_01, lv] == 3) {
        // GP mode: full time-varying residual via calc_gp_pred
        int gp_start, gp_end;
        (gp_start, gp_end) = get_pos(gp_level_pos_ms_baseline_slot[MS_SLOT_01], lv);
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
      for (j in 1:n_forecast_patients) {
        int p = forecast_patient_idx[j];
        log_cond_surv_01[j] += log_level_lambda_01_residual[patient_ms_baseline_flat_idx_slot[MS_SLOT_01, p, lv]];
      }
    }
  }

  // -------------------------------------------------------------------------
  // Time-varying covariate effects
  // -------------------------------------------------------------------------
  // Continuous mode: dense matrix addition across all weeks.
  if (enable_ms_pop_time_varying_cov && n_time_varying_covar > 0 && !enable_ms_visit_gated_01
      && size(ms_time_varying_covar_01) > 0) {
    for (k in 1:n_time_varying_covar) {
      log_cond_surv_01 += time_varying_coef_01[k] * ms_time_varying_covar_01[k];
    }
  }
  // Visit-gated mode: sparse update at observed visit weeks only (before -exp,
  // same log-hazard-level addition as continuous mode — no special handling needed).
  // Two sub-paths:
  //   Latent: read all n_time_varying_covar components of the modeled trajectory
  //     (ms_time_varying_covar_01) at the visit week.
  //   Observed: read the single biomarker stored in ms_obs_visit_covar_flat at
  //     the visit position v.
  // When ms_needs_inline_burden && enable_ms_visit_gated_latent_01, the latent
  // path is handled by the model's inline include — skip here. Observed path
  // still runs.
  if (enable_ms_pop_time_varying_cov && enable_ms_visit_gated_01 &&
      !(ms_needs_inline_burden && enable_ms_visit_gated_latent_01) &&
      (!enable_ms_visit_gated_latent_01 || size(ms_time_varying_covar_01) > 0)) {
    for (j in 1:n_forecast_patients) {
      int p = forecast_patient_idx[j];
      int v_start; int v_end;
      (v_start, v_end) = get_pos(patient_visit_pos, p);
      for (v in v_start:v_end) {
        int wk = t_patient_visits[v];
        if (wk >= 1 && wk <= max_all_t) {
          if (enable_ms_visit_gated_latent_01) {
            // Latent: dot-product over all modeled time-varying covariates
            for (k in 1:n_time_varying_covar) {
              log_cond_surv_01[j, wk] += time_varying_coef_01[k] * ms_time_varying_covar_01[k][j, wk];
            }
          } else {
            // Observed: single biomarker per visit
            log_cond_surv_01[j, wk] += time_varying_coef_01[1] * ms_obs_visit_covar_flat[v];
          }
        }
      }
    }
  }

  // -------------------------------------------------------------------------
  // Time-invariant covariate effects (QR space with multi-level random slopes)
  // -------------------------------------------------------------------------
  if (enable_ms_pop_time_invariant_cov && n_time_invariant_covar > 0) {
    // Population-level covariate effects (QR space)
    vector[n_forecast_patients] linpred_pop_01 = Q_covar_design_matrix[forecast_patient_idx, :] * time_invariant_coef_qr_01;

    // Multi-level random slopes (if enabled)
    if (n_enabled_groups_ms_slope > 0) {
      // Step 1: Scale all slope effects at once, dispatching by mode
      matrix[n_enabled_groups_ms_slope, n_time_invariant_covar] ms_scaled_level_slope_01;
      for (lv in 1:n_levels) {
        if (enable_ms_level_cov[lv]) {
          int lv_start, lv_end;
          (lv_start, lv_end) = get_pos(enabled_level_pos_ms_slope, lv);
          int mode = ms_legacy_mode[MS_SLOT_01, lv];
          if (mode == LEVEL_MODE_RE_CP) {
            int c_lo = cp_level_pos_ms_slope_shared[lv];
            int c_hi = cp_level_pos_ms_slope_shared[lv + 1] - 1;
            ms_scaled_level_slope_01[lv_start:lv_end, :] = cp_level_slope_01[c_lo:c_hi, :];
          } else {
            int r_lo = raw_level_pos_ms_slope_shared[lv];
            int r_hi = raw_level_pos_ms_slope_shared[lv + 1] - 1;
            ms_scaled_level_slope_01[lv_start:lv_end, :] =
              raw_level_slope_01[r_lo:r_hi, :] .*
              rep_matrix(sd_level_slope_01[lv]', lv_end - lv_start + 1);
          }
        }
      }

      // Step 2: Gather and compute dot products using pre-computed flat indices
      for (lv in 1:n_levels) {
        if (enable_ms_level_cov[lv]) {
          linpred_pop_01 += rows_dot_product(
            Q_covar_design_matrix[forecast_patient_idx, :],
            ms_scaled_level_slope_01[patient_ms_slope_flat_idx[forecast_patient_idx, lv], :]
          );
        }
      }
    }

    // Add time-invariant effects (broadcast to all times)
    log_cond_surv_01 += rep_matrix(linpred_pop_01, max_all_t);
  }

  // NOTE: -exp() transform for log_cond_surv_01 is deferred to
  // modules/multistate/cond_surv_transform.stan to allow inline burden
  // TV-covariate insertion (e.g. PSA inline path) before the transform.
}

// ============================================================================
// 0→2 TRANSITION: Log Conditional Survival
// ============================================================================

row_vector[enable_ms_02 ? max_all_t : 0] log_pop_lambda_02;
matrix[n_enabled_groups_ms_baseline_02, enable_ms_02 ? max_all_t : 0] log_level_lambda_02_residual;
vector[n_enabled_groups_ms_baseline_02] log_lambda_gp_02_level_intercept;
matrix[enable_ms_02 ? n_forecast_patients : 0, enable_ms_02 ? max_all_t : 0] log_cond_surv_02;

if (enable_ms_02) {
  // Compute population GP on coarse grid then expand to weekly
  if (share_dead_gp_shape) {
    log_pop_lambda_02 = calc_gp_pred(
      ms_gp_cal_t,
      log_lambda_gp_02_pop_intercept[1],
      log_lambda_gp_dead_pop_alpha[1],
      log_lambda_gp_dead_pop_rho[1],
      delta,
      log_lambda_gp_dead_pop_eta
    )[knot_of_cal];
  } else {
    log_pop_lambda_02 = calc_gp_pred(
      ms_gp_cal_t,
      log_lambda_gp_02_pop_intercept[1],
      log_lambda_gp_02_pop_alpha[1],
      log_lambda_gp_02_pop_rho[1],
      delta,
      log_lambda_gp_02_pop_eta
    )[knot_of_cal];
  }

  // Initialize with population baseline
  log_cond_surv_02 = rep_matrix(log_pop_lambda_02, n_forecast_patients);

  // Add level-level residuals (intercept-only or full GP, same pattern as 0→1)
  for (lv in 1:n_levels) {
    if (ms_legacy_mode[MS_SLOT_02, lv]) {
      int lv_start, lv_end;
      (lv_start, lv_end) = get_pos(enabled_level_pos_ms_baseline_slot[MS_SLOT_02], lv);

      // Scale intercepts: route by mode
      int mode_02 = ms_legacy_mode[MS_SLOT_02, lv];
      if (mode_02 == LEVEL_MODE_RE_CP) {
        int c_lo = cp_level_pos_ms_baseline_slot[MS_SLOT_02][lv];
        int c_hi = cp_level_pos_ms_baseline_slot[MS_SLOT_02][lv + 1] - 1;
        log_lambda_gp_02_level_intercept[lv_start:lv_end] =
          cp_log_lambda_gp_02_level_intercept[c_lo:c_hi];
      } else {
        int r_lo = raw_level_pos_ms_baseline_slot[MS_SLOT_02][lv];
        int r_hi = raw_level_pos_ms_baseline_slot[MS_SLOT_02][lv + 1] - 1;
        if (mode_02 == LEVEL_MODE_FE) {
          log_lambda_gp_02_level_intercept[lv_start:lv_end] =
            raw_log_lambda_gp_02_level_intercept[r_lo:r_hi] *
            fe_log_lambda_gp_02_level_intercept_sd[lv];
        } else if (ms_corr_member_block[MS_SLOT_02, lv] > 0) {
          int b = ms_corr_member_block[MS_SLOT_02, lv];
          int row = ms_corr_member_row[MS_SLOT_02, lv];
          log_lambda_gp_02_level_intercept[lv_start:lv_end] =
            ms_corr_u[b][row, ]';
        } else {
          log_lambda_gp_02_level_intercept[lv_start:lv_end] =
            raw_log_lambda_gp_02_level_intercept[r_lo:r_hi] *
            log_lambda_gp_02_level_intercept_sd[lv];
        }
      }

      if (ms_legacy_mode[MS_SLOT_02, lv] == 3) {
        int gp_start, gp_end;
        (gp_start, gp_end) = get_pos(gp_level_pos_ms_baseline_slot[MS_SLOT_02], lv);
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

      for (j in 1:n_forecast_patients) {
        int p = forecast_patient_idx[j];
        log_cond_surv_02[j] += log_level_lambda_02_residual[patient_ms_baseline_flat_idx_slot[MS_SLOT_02, p, lv]];
      }
    }
  }

  // -------------------------------------------------------------------------
  // Time-varying covariate effects (same pattern as 0→1)
  // -------------------------------------------------------------------------
  // 0->2 time-varying covariate: uses modeled burden (ms_time_varying_covar_01)
  // Controlled by enable_ms_02_time_varying_cov, independent of 0->1 mode.
  // When ms_needs_inline_burden, the inline burden path handles this instead.
  if (!ms_needs_inline_burden && enable_ms_pop_time_varying_cov && enable_ms_02_time_varying_cov
      && n_time_varying_covar > 0 && size(ms_time_varying_covar_01) > 0) {
    for (k in 1:n_time_varying_covar) {
      log_cond_surv_02 += time_varying_coef_02[k] * ms_time_varying_covar_01[k];
    }
  }

  // -------------------------------------------------------------------------
  // Time-invariant covariate effects
  // -------------------------------------------------------------------------
  if (enable_ms_pop_time_invariant_cov && n_time_invariant_covar > 0) {
    vector[n_forecast_patients] linpred_pop_02 = Q_covar_design_matrix[forecast_patient_idx, :] * time_invariant_coef_qr_02;

    if (n_enabled_groups_ms_slope > 0) {
      matrix[n_enabled_groups_ms_slope, n_time_invariant_covar] ms_scaled_level_slope_02;
      for (lv in 1:n_levels) {
        if (enable_ms_level_cov[lv]) {
          int lv_start, lv_end;
          (lv_start, lv_end) = get_pos(enabled_level_pos_ms_slope, lv);
          int mode = ms_legacy_mode[MS_SLOT_02, lv];
          if (mode == LEVEL_MODE_RE_CP) {
            int c_lo = cp_level_pos_ms_slope_shared[lv];
            int c_hi = cp_level_pos_ms_slope_shared[lv + 1] - 1;
            ms_scaled_level_slope_02[lv_start:lv_end, :] = cp_level_slope_02[c_lo:c_hi, :];
          } else {
            int r_lo = raw_level_pos_ms_slope_shared[lv];
            int r_hi = raw_level_pos_ms_slope_shared[lv + 1] - 1;
            ms_scaled_level_slope_02[lv_start:lv_end, :] =
              raw_level_slope_02[r_lo:r_hi, :] .*
              rep_matrix(sd_level_slope_02[lv]', lv_end - lv_start + 1);
          }
        }
      }

      for (lv in 1:n_levels) {
        if (enable_ms_level_cov[lv]) {
          linpred_pop_02 += rows_dot_product(
            Q_covar_design_matrix[forecast_patient_idx, :],
            ms_scaled_level_slope_02[patient_ms_slope_flat_idx[forecast_patient_idx, lv], :]
          );
        }
      }
    }

    log_cond_surv_02 += rep_matrix(linpred_pop_02, max_all_t);
  }

  // NOTE: -exp() transform for log_cond_surv_02 is deferred to
  // modules/multistate/cond_surv_transform.stan (same reason as 0→1).
}

// ============================================================================
// 1→2 TRANSITION: Log Conditional Survival (on sojourn and/or clock-forward time)
// ============================================================================

// Sojourn time GP (semi-Markov or extended)
row_vector[need_12_s_gp ? ms_max_sojourn_t : 0] log_pop_lambda_12_s;
matrix[n_enabled_groups_ms_baseline_12_s, need_12_s_gp ? ms_max_sojourn_t : 0] log_level_lambda_12_s_residual;
vector[n_enabled_groups_ms_baseline_12_s] log_lambda_gp_12_s_level_intercept;
matrix[need_12_s_gp ? n_forecast_patients : 0, need_12_s_gp ? ms_max_sojourn_t : 0] log_cond_surv_12_s;

// Clock-forward time GP (Markov or extended)
row_vector[need_12_t_gp ? max_all_t : 0] log_pop_lambda_12_t;
matrix[n_enabled_groups_ms_baseline_12_t, need_12_t_gp ? max_all_t : 0] log_level_lambda_12_t_residual;
vector[n_enabled_groups_ms_baseline_12_t] log_lambda_gp_12_t_level_intercept;
matrix[need_12_t_gp ? n_forecast_patients : 0, need_12_t_gp ? max_all_t : 0] log_cond_surv_12_t;

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

  log_cond_surv_12_s = rep_matrix(log_pop_lambda_12_s, n_forecast_patients);

  // Add level residuals (intercept-only or full GP)
  for (lv in 1:n_levels) {
    if (ms_legacy_mode[MS_SLOT_12_S, lv]) {
      int lv_start, lv_end;
      (lv_start, lv_end) = get_pos(enabled_level_pos_ms_baseline_slot[MS_SLOT_12_S], lv);

      // Scale intercepts: route by mode
      int mode_12_s = ms_legacy_mode[MS_SLOT_12_S, lv];
      if (mode_12_s == LEVEL_MODE_RE_CP) {
        int c_lo = cp_level_pos_ms_baseline_slot[MS_SLOT_12_S][lv];
        int c_hi = cp_level_pos_ms_baseline_slot[MS_SLOT_12_S][lv + 1] - 1;
        log_lambda_gp_12_s_level_intercept[lv_start:lv_end] =
          cp_log_lambda_gp_12_s_level_intercept[c_lo:c_hi];
      } else {
        int r_lo = raw_level_pos_ms_baseline_slot[MS_SLOT_12_S][lv];
        int r_hi = raw_level_pos_ms_baseline_slot[MS_SLOT_12_S][lv + 1] - 1;
        if (mode_12_s == LEVEL_MODE_FE) {
          log_lambda_gp_12_s_level_intercept[lv_start:lv_end] =
            raw_log_lambda_gp_12_s_level_intercept[r_lo:r_hi] *
            fe_log_lambda_gp_12_s_level_intercept_sd[lv];
        } else if (ms_corr_member_block[MS_SLOT_12_S, lv] > 0) {
          int b = ms_corr_member_block[MS_SLOT_12_S, lv];
          int row = ms_corr_member_row[MS_SLOT_12_S, lv];
          log_lambda_gp_12_s_level_intercept[lv_start:lv_end] =
            ms_corr_u[b][row, ]';
        } else {
          log_lambda_gp_12_s_level_intercept[lv_start:lv_end] =
            raw_log_lambda_gp_12_s_level_intercept[r_lo:r_hi] *
            log_lambda_gp_12_s_level_intercept_sd[lv];
        }
      }

      if (ms_legacy_mode[MS_SLOT_12_S, lv] == 3) {
        int gp_start, gp_end;
        (gp_start, gp_end) = get_pos(gp_level_pos_ms_baseline_slot[MS_SLOT_12_S], lv);
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

      for (j in 1:n_forecast_patients) {
        int p = forecast_patient_idx[j];
        log_cond_surv_12_s[j] += log_level_lambda_12_s_residual[patient_ms_baseline_flat_idx_slot[MS_SLOT_12_S, p, lv]];
      }
    }
  }

  // -------------------------------------------------------------------------
  // Time-invariant covariate effects for 1→2 sojourn
  // -------------------------------------------------------------------------
  // Note: Time-varying covariates for 1→2 would need different indexing (sojourn time)
  // For now, only time-invariant covariates are supported for 1→2
  if (enable_ms_pop_time_invariant_cov && n_time_invariant_covar > 0) {
    vector[n_forecast_patients] linpred_pop_12 = Q_covar_design_matrix[forecast_patient_idx, :] * time_invariant_coef_qr_12;

    if (n_enabled_groups_ms_slope > 0) {
      matrix[n_enabled_groups_ms_slope, n_time_invariant_covar] ms_scaled_level_slope_12;
      for (lv in 1:n_levels) {
        if (enable_ms_level_cov[lv]) {
          int lv_start, lv_end;
          (lv_start, lv_end) = get_pos(enabled_level_pos_ms_slope, lv);
          int mode = ms_legacy_mode[MS_SLOT_12_S, lv];
          if (mode == LEVEL_MODE_RE_CP) {
            int c_lo = cp_level_pos_ms_slope_shared[lv];
            int c_hi = cp_level_pos_ms_slope_shared[lv + 1] - 1;
            ms_scaled_level_slope_12[lv_start:lv_end, :] = cp_level_slope_12[c_lo:c_hi, :];
          } else {
            int r_lo = raw_level_pos_ms_slope_shared[lv];
            int r_hi = raw_level_pos_ms_slope_shared[lv + 1] - 1;
            ms_scaled_level_slope_12[lv_start:lv_end, :] =
              raw_level_slope_12[r_lo:r_hi, :] .*
              rep_matrix(sd_level_slope_12[lv]', lv_end - lv_start + 1);
          }
        }
      }

      for (lv in 1:n_levels) {
        if (enable_ms_level_cov[lv]) {
          linpred_pop_12 += rows_dot_product(
            Q_covar_design_matrix[forecast_patient_idx, :],
            ms_scaled_level_slope_12[patient_ms_slope_flat_idx[forecast_patient_idx, lv], :]
          );
        }
      }
    }

    log_cond_surv_12_s += rep_matrix(linpred_pop_12, ms_max_sojourn_t);
  }

  // Burden-at-entry covariate: shift sojourn hazard per patient based on the
  // (PSA or other) burden value at the moment of progression (entry into state 1).
  if (enable_ms_12_entry_covar) {
    log_cond_surv_12_s += rep_matrix(
      coef_log_entry_covar_12[1] * to_vector(entry_covar_12), ms_max_sojourn_t
    );
  }

  log_cond_surv_12_s = -exp(log_cond_surv_12_s);
}

if (need_12_t_gp) {
  // Clock-forward time GP on coarse grid then expand to weekly
  // In extended mode, intercept is zero (sojourn GP carries it) to avoid non-identifiability
  if (share_dead_gp_shape) {
    log_pop_lambda_12_t = calc_gp_pred(
      ms_gp_cal_t,
      ms_12_t_has_intercept ? log_lambda_gp_12_t_pop_intercept[1] : 0.0,
      log_lambda_gp_dead_pop_alpha[1],
      log_lambda_gp_dead_pop_rho[1],
      delta,
      log_lambda_gp_dead_pop_eta
    )[knot_of_cal];
  } else {
    log_pop_lambda_12_t = calc_gp_pred(
      ms_gp_cal_t,
      ms_12_t_has_intercept ? log_lambda_gp_12_t_pop_intercept[1] : 0.0,
      log_lambda_gp_12_t_pop_alpha[1],
      log_lambda_gp_12_t_pop_rho[1],
      delta,
      log_lambda_gp_12_t_pop_eta
    )[knot_of_cal];
  }

  log_cond_surv_12_t = rep_matrix(log_pop_lambda_12_t, n_forecast_patients);

  // Add level residuals (intercept-only or full GP)
  for (lv in 1:n_levels) {
    if (ms_legacy_mode[MS_SLOT_12_T, lv]) {
      int lv_start, lv_end;
      (lv_start, lv_end) = get_pos(enabled_level_pos_ms_baseline_slot[MS_SLOT_12_T], lv);

      // Scale intercepts: route by mode
      int mode_12_t = ms_legacy_mode[MS_SLOT_12_T, lv];
      if (mode_12_t == LEVEL_MODE_RE_CP) {
        int c_lo = cp_level_pos_ms_baseline_slot[MS_SLOT_12_T][lv];
        int c_hi = cp_level_pos_ms_baseline_slot[MS_SLOT_12_T][lv + 1] - 1;
        log_lambda_gp_12_t_level_intercept[lv_start:lv_end] =
          cp_log_lambda_gp_12_t_level_intercept[c_lo:c_hi];
      } else {
        int r_lo = raw_level_pos_ms_baseline_slot[MS_SLOT_12_T][lv];
        int r_hi = raw_level_pos_ms_baseline_slot[MS_SLOT_12_T][lv + 1] - 1;
        if (mode_12_t == LEVEL_MODE_FE) {
          log_lambda_gp_12_t_level_intercept[lv_start:lv_end] =
            raw_log_lambda_gp_12_t_level_intercept[r_lo:r_hi] *
            fe_log_lambda_gp_12_t_level_intercept_sd[lv];
        } else if (ms_corr_member_block[MS_SLOT_12_T, lv] > 0) {
          int b = ms_corr_member_block[MS_SLOT_12_T, lv];
          int row = ms_corr_member_row[MS_SLOT_12_T, lv];
          log_lambda_gp_12_t_level_intercept[lv_start:lv_end] =
            ms_corr_u[b][row, ]';
        } else {
          log_lambda_gp_12_t_level_intercept[lv_start:lv_end] =
            raw_log_lambda_gp_12_t_level_intercept[r_lo:r_hi] *
            log_lambda_gp_12_t_level_intercept_sd[lv];
        }
      }

      if (ms_legacy_mode[MS_SLOT_12_T, lv] == 3) {
        int gp_start, gp_end;
        (gp_start, gp_end) = get_pos(gp_level_pos_ms_baseline_slot[MS_SLOT_12_T], lv);
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

      for (j in 1:n_forecast_patients) {
        int p = forecast_patient_idx[j];
        log_cond_surv_12_t[j] += log_level_lambda_12_t_residual[patient_ms_baseline_flat_idx_slot[MS_SLOT_12_T, p, lv]];
      }
    }
  }

  // -------------------------------------------------------------------------
  // Time-invariant covariate effects for 1→2 clock-forward
  // -------------------------------------------------------------------------
  if (enable_ms_pop_time_invariant_cov && n_time_invariant_covar > 0) {
    vector[n_forecast_patients] linpred_pop_12 = Q_covar_design_matrix[forecast_patient_idx, :] * time_invariant_coef_qr_12;

    if (n_enabled_groups_ms_slope > 0) {
      matrix[n_enabled_groups_ms_slope, n_time_invariant_covar] ms_scaled_level_slope_12;
      for (lv in 1:n_levels) {
        if (enable_ms_level_cov[lv]) {
          int lv_start, lv_end;
          (lv_start, lv_end) = get_pos(enabled_level_pos_ms_slope, lv);
          int mode = ms_legacy_mode[MS_SLOT_12_T, lv];
          if (mode == LEVEL_MODE_RE_CP) {
            int c_lo = cp_level_pos_ms_slope_shared[lv];
            int c_hi = cp_level_pos_ms_slope_shared[lv + 1] - 1;
            ms_scaled_level_slope_12[lv_start:lv_end, :] = cp_level_slope_12[c_lo:c_hi, :];
          } else {
            int r_lo = raw_level_pos_ms_slope_shared[lv];
            int r_hi = raw_level_pos_ms_slope_shared[lv + 1] - 1;
            ms_scaled_level_slope_12[lv_start:lv_end, :] =
              raw_level_slope_12[r_lo:r_hi, :] .*
              rep_matrix(sd_level_slope_12[lv]', lv_end - lv_start + 1);
          }
        }
      }

      for (lv in 1:n_levels) {
        if (enable_ms_level_cov[lv]) {
          linpred_pop_12 += rows_dot_product(
            Q_covar_design_matrix[forecast_patient_idx, :],
            ms_scaled_level_slope_12[patient_ms_slope_flat_idx[forecast_patient_idx, lv], :]
          );
        }
      }
    }

    log_cond_surv_12_t += rep_matrix(linpred_pop_12, max_all_t);
  }

  log_cond_surv_12_t = -exp(log_cond_surv_12_t);
}

// ============================================================================
// 0→3 TRANSITION: GP baseline hazard (clock-forward time) with N-level hierarchy
// ============================================================================
row_vector[enable_ms_03 ? max_all_t : 0] log_pop_lambda_03;
matrix[n_enabled_groups_ms_baseline_03, enable_ms_03 ? max_all_t : 0] log_level_lambda_03_residual;
vector[n_enabled_groups_ms_baseline_03] log_lambda_gp_03_level_intercept;
matrix[enable_ms_03 ? n_forecast_patients : 0, enable_ms_03 ? max_all_t : 0] log_cond_surv_03;

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
  log_cond_surv_03 = rep_matrix(log_pop_lambda_03, n_forecast_patients);

  // Add level-level residuals (intercept-only or full GP)
  for (lv in 1:n_levels) {
    if (ms_legacy_mode[MS_SLOT_03, lv]) {
      int lv_start, lv_end;
      (lv_start, lv_end) = get_pos(enabled_level_pos_ms_baseline_slot[MS_SLOT_03], lv);

      // Scale intercepts: route by mode
      int mode_03 = ms_legacy_mode[MS_SLOT_03, lv];
      if (mode_03 == LEVEL_MODE_RE_CP) {
        int c_lo = cp_level_pos_ms_baseline_slot[MS_SLOT_03][lv];
        int c_hi = cp_level_pos_ms_baseline_slot[MS_SLOT_03][lv + 1] - 1;
        log_lambda_gp_03_level_intercept[lv_start:lv_end] =
          cp_log_lambda_gp_03_level_intercept[c_lo:c_hi];
      } else {
        int r_lo = raw_level_pos_ms_baseline_slot[MS_SLOT_03][lv];
        int r_hi = raw_level_pos_ms_baseline_slot[MS_SLOT_03][lv + 1] - 1;
        if (mode_03 == LEVEL_MODE_FE) {
          log_lambda_gp_03_level_intercept[lv_start:lv_end] =
            raw_log_lambda_gp_03_level_intercept[r_lo:r_hi] *
            fe_log_lambda_gp_03_level_intercept_sd[lv];
        } else if (ms_corr_member_block[MS_SLOT_03, lv] > 0) {
          // Correlated RE member: take this transition's MVN row.
          int b = ms_corr_member_block[MS_SLOT_03, lv];
          int row = ms_corr_member_row[MS_SLOT_03, lv];
          log_lambda_gp_03_level_intercept[lv_start:lv_end] =
            ms_corr_u[b][row, ]';
        } else {
          log_lambda_gp_03_level_intercept[lv_start:lv_end] =
            raw_log_lambda_gp_03_level_intercept[r_lo:r_hi] *
            log_lambda_gp_03_level_intercept_sd[lv];
        }
      }

      if (ms_legacy_mode[MS_SLOT_03, lv] == 3) {
        int gp_start, gp_end;
        (gp_start, gp_end) = get_pos(gp_level_pos_ms_baseline_slot[MS_SLOT_03], lv);
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

      for (j in 1:n_forecast_patients) {
        int p = forecast_patient_idx[j];
        log_cond_surv_03[j] += log_level_lambda_03_residual[patient_ms_baseline_flat_idx_slot[MS_SLOT_03, p, lv]];
      }
    }
  }

  // -------------------------------------------------------------------------
  // Time-varying covariate effects for 0->3 (modeled burden trajectory)
  // -------------------------------------------------------------------------
  // Uses ms_time_varying_covar_01 (modeled biomarker matrix; PSA for pioneer,
  // SLD-derived for tumor models), independent of 0->1 mode.
  if (!ms_needs_inline_burden && enable_ms_pop_time_varying_cov && enable_ms_03_time_varying_cov
      && n_time_varying_covar > 0 && size(ms_time_varying_covar_01) > 0) {
    for (k in 1:n_time_varying_covar) {
      log_cond_surv_03 += time_varying_coef_03[k] * ms_time_varying_covar_01[k];
    }
  }

  // -------------------------------------------------------------------------
  // Time-invariant covariate effects (QR space with multi-level random slopes)
  // -------------------------------------------------------------------------
  if (enable_ms_pop_time_invariant_cov && enable_ms_03_time_invariant_cov && n_time_invariant_covar > 0) {
    vector[n_forecast_patients] linpred_pop_03 = Q_covar_design_matrix[forecast_patient_idx, :] * time_invariant_coef_qr_03;

    if (n_enabled_groups_ms_slope > 0) {
      matrix[n_enabled_groups_ms_slope, n_time_invariant_covar] ms_scaled_level_slope_03;
      for (lv in 1:n_levels) {
        if (enable_ms_level_cov[lv]) {
          int lv_start, lv_end;
          (lv_start, lv_end) = get_pos(enabled_level_pos_ms_slope, lv);
          int mode = ms_legacy_mode[MS_SLOT_03, lv];
          if (mode == LEVEL_MODE_RE_CP) {
            int c_lo = cp_level_pos_ms_slope_shared[lv];
            int c_hi = cp_level_pos_ms_slope_shared[lv + 1] - 1;
            ms_scaled_level_slope_03[lv_start:lv_end, :] = cp_level_slope_03[c_lo:c_hi, :];
          } else {
            int r_lo = raw_level_pos_ms_slope_shared[lv];
            int r_hi = raw_level_pos_ms_slope_shared[lv + 1] - 1;
            ms_scaled_level_slope_03[lv_start:lv_end, :] =
              raw_level_slope_03[r_lo:r_hi, :] .*
              rep_matrix(sd_level_slope_03[lv]', lv_end - lv_start + 1);
          }
        }
      }

      for (lv in 1:n_levels) {
        if (enable_ms_level_cov[lv]) {
          linpred_pop_03 += rows_dot_product(
            Q_covar_design_matrix[forecast_patient_idx, :],
            ms_scaled_level_slope_03[patient_ms_slope_flat_idx[forecast_patient_idx, lv], :]
          );
        }
      }
    }

    log_cond_surv_03 += rep_matrix(linpred_pop_03, max_all_t);
  }
  // log_cond_surv_03 is transformed by cond_surv_transform.stan (after any
  // inline burden TV-covariate contributions are added).
}

// ============================================================================
// 3→2 TRANSITION: Sojourn time GP baseline hazard (semi-Markov)
// ============================================================================

row_vector[enable_ms_32 ? ms_max_sojourn_t_32 : 0] log_pop_lambda_32;
matrix[n_enabled_groups_ms_baseline_32, enable_ms_32 ? ms_max_sojourn_t_32 : 0] log_level_lambda_32_residual;
vector[n_enabled_groups_ms_baseline_32] log_lambda_gp_32_s_level_intercept;
matrix[enable_ms_32 ? n_forecast_patients : 0, enable_ms_32 ? ms_max_sojourn_t_32 : 0] log_cond_surv_32;

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
  log_cond_surv_32 = rep_matrix(log_pop_lambda_32, n_forecast_patients);

  // Add level-level residuals (intercept-only or full GP)
  for (lv in 1:n_levels) {
    if (ms_legacy_mode[MS_SLOT_32, lv]) {
      int lv_start, lv_end;
      (lv_start, lv_end) = get_pos(enabled_level_pos_ms_baseline_slot[MS_SLOT_32], lv);

      // Scale intercepts: route by mode
      int mode_32 = ms_legacy_mode[MS_SLOT_32, lv];
      if (mode_32 == LEVEL_MODE_RE_CP) {
        int c_lo = cp_level_pos_ms_baseline_slot[MS_SLOT_32][lv];
        int c_hi = cp_level_pos_ms_baseline_slot[MS_SLOT_32][lv + 1] - 1;
        log_lambda_gp_32_s_level_intercept[lv_start:lv_end] =
          cp_log_lambda_gp_32_s_level_intercept[c_lo:c_hi];
      } else {
        int r_lo = raw_level_pos_ms_baseline_slot[MS_SLOT_32][lv];
        int r_hi = raw_level_pos_ms_baseline_slot[MS_SLOT_32][lv + 1] - 1;
        if (mode_32 == LEVEL_MODE_FE) {
          log_lambda_gp_32_s_level_intercept[lv_start:lv_end] =
            raw_log_lambda_gp_32_s_level_intercept[r_lo:r_hi] *
            fe_log_lambda_gp_32_s_level_intercept_sd[lv];
        } else if (ms_corr_member_block[MS_SLOT_32, lv] > 0) {
          int b = ms_corr_member_block[MS_SLOT_32, lv];
          int row = ms_corr_member_row[MS_SLOT_32, lv];
          log_lambda_gp_32_s_level_intercept[lv_start:lv_end] =
            ms_corr_u[b][row, ]';
        } else {
          log_lambda_gp_32_s_level_intercept[lv_start:lv_end] =
            raw_log_lambda_gp_32_s_level_intercept[r_lo:r_hi] *
            log_lambda_gp_32_s_level_intercept_sd[lv];
        }
      }

      if (ms_legacy_mode[MS_SLOT_32, lv] == 3) {
        int gp_start, gp_end;
        (gp_start, gp_end) = get_pos(gp_level_pos_ms_baseline_slot[MS_SLOT_32], lv);
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

      for (j in 1:n_forecast_patients) {
        int p = forecast_patient_idx[j];
        log_cond_surv_32[j] += log_level_lambda_32_residual[patient_ms_baseline_flat_idx_slot[MS_SLOT_32, p, lv]];
      }
    }
  }

  // -------------------------------------------------------------------------
  // Time-invariant covariate effects (QR space with multi-level random slopes)
  // -------------------------------------------------------------------------
  // Time-varying covariates not implemented for 3->2: would require sojourn-clock
  // indexing of the modeled burden matrix.
  if (enable_ms_pop_time_invariant_cov && enable_ms_32_time_invariant_cov && n_time_invariant_covar > 0) {
    vector[n_forecast_patients] linpred_pop_32 = Q_covar_design_matrix[forecast_patient_idx, :] * time_invariant_coef_qr_32;

    if (n_enabled_groups_ms_slope > 0) {
      matrix[n_enabled_groups_ms_slope, n_time_invariant_covar] ms_scaled_level_slope_32;
      for (lv in 1:n_levels) {
        if (enable_ms_level_cov[lv]) {
          int lv_start, lv_end;
          (lv_start, lv_end) = get_pos(enabled_level_pos_ms_slope, lv);
          int mode = ms_legacy_mode[MS_SLOT_32, lv];
          if (mode == LEVEL_MODE_RE_CP) {
            int c_lo = cp_level_pos_ms_slope_shared[lv];
            int c_hi = cp_level_pos_ms_slope_shared[lv + 1] - 1;
            ms_scaled_level_slope_32[lv_start:lv_end, :] = cp_level_slope_32[c_lo:c_hi, :];
          } else {
            int r_lo = raw_level_pos_ms_slope_shared[lv];
            int r_hi = raw_level_pos_ms_slope_shared[lv + 1] - 1;
            ms_scaled_level_slope_32[lv_start:lv_end, :] =
              raw_level_slope_32[r_lo:r_hi, :] .*
              rep_matrix(sd_level_slope_32[lv]', lv_end - lv_start + 1);
          }
        }
      }

      for (lv in 1:n_levels) {
        if (enable_ms_level_cov[lv]) {
          linpred_pop_32 += rows_dot_product(
            Q_covar_design_matrix[forecast_patient_idx, :],
            ms_scaled_level_slope_32[patient_ms_slope_flat_idx[forecast_patient_idx, lv], :]
          );
        }
      }
    }

    log_cond_surv_32 += rep_matrix(linpred_pop_32, ms_max_sojourn_t_32);
  }

  // Burden-at-entry covariate: shift sojourn hazard per patient based on the
  // (PSA or other) burden value at the moment of dropout (entry into state 3).
  if (enable_ms_32_entry_covar) {
    log_cond_surv_32 += rep_matrix(
      coef_log_entry_covar_32[1] * to_vector(entry_covar_32), ms_max_sojourn_t_32
    );
  }

  // Transform to log conditional survival
  log_cond_surv_32 = -exp(log_cond_surv_32);
}

