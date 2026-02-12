// ============================================================================
// Multistate Hazard Model Transformed Parameters
// ============================================================================

// ============================================================================
// 0→1 TRANSITION: Log Conditional Survival
// ============================================================================

// Population-level baseline hazard
row_vector[enable_ms_01 ? max_all_t : 0] ms_log_pop_lambda_01;

// Level-level baseline hazard (additive to population)
matrix[n_enabled_groups_ms_baseline_01, enable_ms_01 ? max_all_t : 0] ms_log_level_lambda_01_residual;
vector[n_enabled_groups_ms_baseline_01] ms_log_lambda_gp_01_level_intercept;

// Patient-level log conditional survival probability
matrix[enable_ms_01 ? n_patients : 0, enable_ms_01 ? max_all_t : 0] ms_log_cond_surv_01;

if (enable_ms_01) {
  // Compute population GP
  ms_log_pop_lambda_01 = calc_gp_pred(
    all_tumor_measure_t,
    ms_log_lambda_gp_01_pop_intercept[1],
    ms_log_lambda_gp_01_pop_alpha[1],
    ms_log_lambda_gp_01_pop_rho[1],
    delta,
    ms_log_lambda_gp_01_pop_eta
  );

  // Initialize patient hazards with population baseline
  ms_log_cond_surv_01 = rep_matrix(ms_log_pop_lambda_01, n_patients);

  // Add level-level GP residuals
  for (lv in 1:n_levels) {
    if (enable_ms_level_baseline_hazard[lv]) {
      int lv_start, lv_end;
      (lv_start, lv_end) = get_pos(enabled_level_pos_ms_baseline, lv);

      // Scale intercepts
      ms_log_lambda_gp_01_level_intercept[lv_start:lv_end] =
        ms_raw_log_lambda_gp_01_level_intercept[lv_start:lv_end] *
        ms_log_lambda_gp_01_level_intercept_sd[lv];

      // Compute level GP residuals
      for (g in lv_start:lv_end) {
        ms_log_level_lambda_01_residual[g] = calc_gp_pred(
          all_tumor_measure_t,
          ms_log_lambda_gp_01_level_intercept[g],
          ms_log_lambda_gp_01_level_alpha[lv],
          ms_log_lambda_gp_01_level_rho[lv],
          delta,
          ms_log_lambda_gp_01_level_eta[g]
        );
      }

      // Add level residuals to patient hazards
      for (i in 1:n_patients) {
        ms_log_cond_surv_01[i] += ms_log_level_lambda_01_residual[patient_ms_baseline_flat_idx[i, lv]];
      }
    }
  }

  // TODO: Add time-varying and time-invariant covariate effects
  // For Phase 1, we just use baseline hazard

  // Transform log-hazard to log conditional survival probability
  // log P(survive interval t) = -exp(log_hazard[t]) = -hazard[t]
  ms_log_cond_surv_01 = -exp(ms_log_cond_surv_01);
}

// ============================================================================
// 0→2 TRANSITION: Log Conditional Survival
// ============================================================================

row_vector[enable_ms_02 ? max_all_t : 0] ms_log_pop_lambda_02;
matrix[n_enabled_groups_ms_baseline_02, enable_ms_02 ? max_all_t : 0] ms_log_level_lambda_02_residual;
vector[n_enabled_groups_ms_baseline_02] ms_log_lambda_gp_02_level_intercept;
matrix[enable_ms_02 ? n_patients : 0, enable_ms_02 ? max_all_t : 0] ms_log_cond_surv_02;

if (enable_ms_02) {
  // Compute population GP
  ms_log_pop_lambda_02 = calc_gp_pred(
    all_tumor_measure_t,
    ms_log_lambda_gp_02_pop_intercept[1],
    ms_log_lambda_gp_02_pop_alpha[1],
    ms_log_lambda_gp_02_pop_rho[1],
    delta,
    ms_log_lambda_gp_02_pop_eta
  );

  // Initialize with population baseline
  ms_log_cond_surv_02 = rep_matrix(ms_log_pop_lambda_02, n_patients);

  // Add level-level GP residuals (same pattern as 0→1)
  for (lv in 1:n_levels) {
    if (enable_ms_level_baseline_hazard[lv]) {
      int lv_start, lv_end;
      (lv_start, lv_end) = get_pos(enabled_level_pos_ms_baseline, lv);

      ms_log_lambda_gp_02_level_intercept[lv_start:lv_end] =
        ms_raw_log_lambda_gp_02_level_intercept[lv_start:lv_end] *
        ms_log_lambda_gp_02_level_intercept_sd[lv];

      for (g in lv_start:lv_end) {
        ms_log_level_lambda_02_residual[g] = calc_gp_pred(
          all_tumor_measure_t,
          ms_log_lambda_gp_02_level_intercept[g],
          ms_log_lambda_gp_02_level_alpha[lv],
          ms_log_lambda_gp_02_level_rho[lv],
          delta,
          ms_log_lambda_gp_02_level_eta[g]
        );
      }

      for (i in 1:n_patients) {
        ms_log_cond_surv_02[i] += ms_log_level_lambda_02_residual[patient_ms_baseline_flat_idx[i, lv]];
      }
    }
  }

  // Transform to log conditional survival
  ms_log_cond_surv_02 = -exp(ms_log_cond_surv_02);
}

// ============================================================================
// 1→2 TRANSITION: Log Conditional Survival (on sojourn and/or clock-forward time)
// ============================================================================

// Sojourn time GP (semi-Markov or extended)
row_vector[need_12_s_gp ? ms_max_sojourn_t : 0] ms_log_pop_lambda_12_s;
matrix[n_enabled_groups_ms_baseline_12_s, need_12_s_gp ? ms_max_sojourn_t : 0] ms_log_level_lambda_12_s_residual;
vector[n_enabled_groups_ms_baseline_12_s] ms_log_lambda_gp_12_s_level_intercept;
matrix[need_12_s_gp ? n_patients : 0, need_12_s_gp ? ms_max_sojourn_t : 0] ms_log_cond_surv_12_s;

// Clock-forward time GP (Markov or extended)
row_vector[need_12_t_gp ? max_all_t : 0] ms_log_pop_lambda_12_t;
matrix[n_enabled_groups_ms_baseline_12_t, need_12_t_gp ? max_all_t : 0] ms_log_level_lambda_12_t_residual;
vector[n_enabled_groups_ms_baseline_12_t] ms_log_lambda_gp_12_t_level_intercept;
matrix[need_12_t_gp ? n_patients : 0, need_12_t_gp ? max_all_t : 0] ms_log_cond_surv_12_t;

if (need_12_s_gp) {
  // Sojourn time GP
  // Create sojourn time grid
  array[ms_max_sojourn_t] int sojourn_time_grid;
  for (s in 1:ms_max_sojourn_t) sojourn_time_grid[s] = s;

  ms_log_pop_lambda_12_s = calc_gp_pred(
    sojourn_time_grid,
    ms_log_lambda_gp_12_s_pop_intercept[1],
    ms_log_lambda_gp_12_s_pop_alpha[1],
    ms_log_lambda_gp_12_s_pop_rho[1],
    delta,
    ms_log_lambda_gp_12_s_pop_eta
  );

  ms_log_cond_surv_12_s = rep_matrix(ms_log_pop_lambda_12_s, n_patients);

  // Add level residuals (same pattern)
  for (lv in 1:n_levels) {
    if (enable_ms_level_baseline_hazard[lv]) {
      int lv_start, lv_end;
      (lv_start, lv_end) = get_pos(enabled_level_pos_ms_baseline, lv);

      ms_log_lambda_gp_12_s_level_intercept[lv_start:lv_end] =
        ms_raw_log_lambda_gp_12_s_level_intercept[lv_start:lv_end] *
        ms_log_lambda_gp_12_s_level_intercept_sd[lv];

      for (g in lv_start:lv_end) {
        ms_log_level_lambda_12_s_residual[g] = calc_gp_pred(
          sojourn_time_grid,
          ms_log_lambda_gp_12_s_level_intercept[g],
          ms_log_lambda_gp_12_s_level_alpha[lv],
          ms_log_lambda_gp_12_s_level_rho[lv],
          delta,
          ms_log_lambda_gp_12_s_level_eta[g]
        );
      }

      for (i in 1:n_patients) {
        ms_log_cond_surv_12_s[i] += ms_log_level_lambda_12_s_residual[patient_ms_baseline_flat_idx[i, lv]];
      }
    }
  }

  ms_log_cond_surv_12_s = -exp(ms_log_cond_surv_12_s);
}

if (need_12_t_gp) {
  // Clock-forward time GP
  ms_log_pop_lambda_12_t = calc_gp_pred(
    all_tumor_measure_t,
    ms_log_lambda_gp_12_t_pop_intercept[1],
    ms_log_lambda_gp_12_t_pop_alpha[1],
    ms_log_lambda_gp_12_t_pop_rho[1],
    delta,
    ms_log_lambda_gp_12_t_pop_eta
  );

  ms_log_cond_surv_12_t = rep_matrix(ms_log_pop_lambda_12_t, n_patients);

  // Add level residuals
  for (lv in 1:n_levels) {
    if (enable_ms_level_baseline_hazard[lv]) {
      int lv_start, lv_end;
      (lv_start, lv_end) = get_pos(enabled_level_pos_ms_baseline, lv);

      ms_log_lambda_gp_12_t_level_intercept[lv_start:lv_end] =
        ms_raw_log_lambda_gp_12_t_level_intercept[lv_start:lv_end] *
        ms_log_lambda_gp_12_t_level_intercept_sd[lv];

      for (g in lv_start:lv_end) {
        ms_log_level_lambda_12_t_residual[g] = calc_gp_pred(
          all_tumor_measure_t,
          ms_log_lambda_gp_12_t_level_intercept[g],
          ms_log_lambda_gp_12_t_level_alpha[lv],
          ms_log_lambda_gp_12_t_level_rho[lv],
          delta,
          ms_log_lambda_gp_12_t_level_eta[g]
        );
      }

      for (i in 1:n_patients) {
        ms_log_cond_surv_12_t[i] += ms_log_level_lambda_12_t_residual[patient_ms_baseline_flat_idx[i, lv]];
      }
    }
  }

  ms_log_cond_surv_12_t = -exp(ms_log_cond_surv_12_t);
}
