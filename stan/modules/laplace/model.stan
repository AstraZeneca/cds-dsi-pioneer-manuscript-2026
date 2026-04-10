// laplace/model.stan — Hand-coded Laplace marginal for non-target patients
// Included in the model block after existing priors and likelihood.
// Pre-computes per-patient offsets and base hazards, then calls per-patient
// Laplace function that uses Newton solver with analytical Jacobians.

if (enable_laplace_nontarget && n_background_patients > 0) {
  profile("laplace nontarget") {
    // =====================================================================
    // (A) Pre-compute per-patient PSA offsets: pop + group-level intercepts
    // =====================================================================
    vector[n_background_patients] lap_tr_offset = rep_vector(tr_loc_pop, n_background_patients);
    vector[n_background_patients] lap_frac_offset = rep_vector(frac_logit_loc_pop, n_background_patients);
    vector[n_background_patients] lap_init_offset = rep_vector(init_logit_loc_pop, n_background_patients);

    // Add group-level intercepts (trial, arm, etc.) — skip patient level (n_levels)
    for (lv in 1:(n_levels - 1)) {
      if (enable_level_intercept_tr[lv]) {
        for (i in 1:n_background_patients)
          lap_tr_offset[i] += tr_scaled_level_intercept[laplace_tr_intercept_flat_idx[i, lv]];
      }
      if (enable_level_intercept_frac[lv]) {
        for (i in 1:n_background_patients)
          lap_frac_offset[i] += frac_scaled_level_intercept[laplace_frac_intercept_flat_idx[i, lv]];
      }
      if (enable_level_intercept_init[lv]) {
        for (i in 1:n_background_patients)
          lap_init_offset[i] += init_scaled_level_intercept[laplace_init_intercept_flat_idx[i, lv]];
      }
    }

    // Add population covariate effects (QR-transformed)
    if (n_covar > 0) {
      if (enable_pop_cov_tr)
        lap_tr_offset += laplace_Q_covar * tr_coef_qr_pop;
      if (enable_pop_cov_frac)
        lap_frac_offset += laplace_Q_covar * frac_coef_qr_pop;
      if (enable_pop_cov_init)
        lap_init_offset += laplace_Q_covar * init_coef_qr_pop;
    }

    // Add group-level covariate slopes (skip patient level)
    if (n_covar > 0) {
      for (lv in 1:(n_levels - 1)) {
        if (enable_level_cov_tr[lv]) {
          lap_tr_offset += rows_dot_product(
            laplace_Q_covar,
            tr_scaled_level_slope[laplace_tr_slope_flat_idx[, lv], :]
          );
        }
        if (enable_level_cov_frac[lv]) {
          lap_frac_offset += rows_dot_product(
            laplace_Q_covar,
            frac_scaled_level_slope[laplace_frac_slope_flat_idx[, lv], :]
          );
        }
        if (enable_level_cov_init[lv]) {
          lap_init_offset += rows_dot_product(
            laplace_Q_covar,
            init_scaled_level_slope[laplace_init_slope_flat_idx[, lv], :]
          );
        }
      }
    }

    // =====================================================================
    // (B) Pre-compute MS base log-hazard for non-target patients (0→1, 0→2)
    // =====================================================================
    // Base = pop GP + level residuals + time-invariant covariates
    // Time-varying covariates are z-dependent → computed inside Newton solver

    matrix[n_background_patients, enable_ms_01 ? max_all_t : 0] lap_base_log_hazard_01;
    if (enable_ms_01) {
      lap_base_log_hazard_01 = rep_matrix(log_pop_lambda_01, n_background_patients);

      // Level residuals
      for (lv in 1:n_levels) {
        if (enable_ms_level_baseline_hazard[lv]) {
          for (i in 1:n_background_patients)
            lap_base_log_hazard_01[i] +=
              log_level_lambda_01_residual[laplace_ms_baseline_flat_idx[i, lv]];
        }
      }

      // Time-invariant covariates (QR space + multi-level random slopes)
      if (enable_ms_pop_time_invariant_cov && n_time_invariant_covar > 0) {
        vector[n_background_patients] linpred_01 = laplace_Q_covar * time_invariant_coef_qr_01;

        if (n_enabled_groups_ms_slope > 0) {
          matrix[n_enabled_groups_ms_slope, n_time_invariant_covar] lap_ms_scaled_slope_01;
          for (lv in 1:n_levels) {
            if (enable_ms_level_cov[lv]) {
              int lv_start, lv_end;
              (lv_start, lv_end) = get_pos(enabled_level_pos_ms_slope, lv);
              lap_ms_scaled_slope_01[lv_start:lv_end, :] =
                raw_level_slope_01[lv_start:lv_end, :] .*
                rep_matrix(sd_level_slope_01[lv]', lv_end - lv_start + 1);
            }
          }
          for (lv in 1:n_levels) {
            if (enable_ms_level_cov[lv]) {
              linpred_01 += rows_dot_product(
                laplace_Q_covar,
                lap_ms_scaled_slope_01[laplace_ms_slope_flat_idx[, lv], :]
              );
            }
          }
        }

        lap_base_log_hazard_01 += rep_matrix(linpred_01, max_all_t);
      }

    }

    // Mirror for 0→2
    matrix[n_background_patients, enable_ms_02 ? max_all_t : 0] lap_base_log_hazard_02;
    if (enable_ms_02) {
      lap_base_log_hazard_02 = rep_matrix(log_pop_lambda_02, n_background_patients);

      for (lv in 1:n_levels) {
        if (enable_ms_level_baseline_hazard[lv]) {
          for (i in 1:n_background_patients)
            lap_base_log_hazard_02[i] +=
              log_level_lambda_02_residual[laplace_ms_baseline_flat_idx[i, lv]];
        }
      }

      if (enable_ms_pop_time_invariant_cov && n_time_invariant_covar > 0) {
        vector[n_background_patients] linpred_02 = laplace_Q_covar * time_invariant_coef_qr_02;

        if (n_enabled_groups_ms_slope > 0) {
          matrix[n_enabled_groups_ms_slope, n_time_invariant_covar] lap_ms_scaled_slope_02;
          for (lv in 1:n_levels) {
            if (enable_ms_level_cov[lv]) {
              int lv_start, lv_end;
              (lv_start, lv_end) = get_pos(enabled_level_pos_ms_slope, lv);
              lap_ms_scaled_slope_02[lv_start:lv_end, :] =
                raw_level_slope_02[lv_start:lv_end, :] .*
                rep_matrix(sd_level_slope_02[lv]', lv_end - lv_start + 1);
            }
          }
          for (lv in 1:n_levels) {
            if (enable_ms_level_cov[lv]) {
              linpred_02 += rows_dot_product(
                laplace_Q_covar,
                lap_ms_scaled_slope_02[laplace_ms_slope_flat_idx[, lv], :]
              );
            }
          }
        }

        lap_base_log_hazard_02 += rep_matrix(linpred_02, max_all_t);
      }

    }

    // =====================================================================
    // (C) z-independent MS loglik: 1→2 transition for non-target patients
    // =====================================================================
    // 1→2 hazard doesn't depend on PSA z_i → pre-compute and add as constant
    real lap_zi_independent_ll = 0;
    if (enable_ms_12) {
      // Build log_cond_surv_12_s for non-target patients
      // Uses the SAME GP arrays as target patients (already computed)
      int lap_sojourn_cols = need_12_s_gp ? ms_max_sojourn_t : 0;
      matrix[need_12_s_gp ? n_background_patients : 0, lap_sojourn_cols] lap_log_cond_surv_12_s;

      if (need_12_s_gp) {
        lap_log_cond_surv_12_s = rep_matrix(log_pop_lambda_12_s, n_background_patients);

        for (lv in 1:n_levels) {
          if (enable_ms_level_baseline_hazard[lv]) {
            for (i in 1:n_background_patients)
              lap_log_cond_surv_12_s[i] +=
                log_level_lambda_12_s_residual[laplace_ms_baseline_flat_idx[i, lv]];
          }
        }

        // Time-invariant covariates for 1→2
        if (enable_ms_pop_time_invariant_cov && n_time_invariant_covar > 0) {
          vector[n_background_patients] linpred_12 = laplace_Q_covar * time_invariant_coef_qr_12;

          if (n_enabled_groups_ms_slope > 0) {
            matrix[n_enabled_groups_ms_slope, n_time_invariant_covar] lap_ms_scaled_slope_12;
            for (lv in 1:n_levels) {
              if (enable_ms_level_cov[lv]) {
                int lv_start, lv_end;
                (lv_start, lv_end) = get_pos(enabled_level_pos_ms_slope, lv);
                lap_ms_scaled_slope_12[lv_start:lv_end, :] =
                  raw_level_slope_12[lv_start:lv_end, :] .*
                  rep_matrix(sd_level_slope_12[lv]', lv_end - lv_start + 1);
              }
            }
            for (lv in 1:n_levels) {
              if (enable_ms_level_cov[lv]) {
                linpred_12 += rows_dot_product(
                  laplace_Q_covar,
                  lap_ms_scaled_slope_12[laplace_ms_slope_flat_idx[, lv], :]
                );
              }
            }
          }

          lap_log_cond_surv_12_s += rep_matrix(linpred_12, ms_max_sojourn_t);
        }

        // Floor at -20: prevents log1m_exp(0) = -Inf (gradient = -Inf → NaN)
        lap_log_cond_surv_12_s = -exp(fmax(lap_log_cond_surv_12_s, -20.0));
      }

      lap_zi_independent_ll = laplace_ms_zi_independent_ll(
        n_background_patients,
        ms_final_state[background_patient_idx],
        ms_time_01[background_patient_idx], ms_time_12[background_patient_idx],
        lap_log_cond_surv_12_s,
        enable_ms_12, ms_time_scale_12
      );
    }

    target += lap_zi_independent_ll;

    // =====================================================================
    // (D) Main Laplace loop — parallelised via reduce_sum (TBB threads)
    // =====================================================================
    // grainsize=1: each patient is an independent Newton solve; TBB auto-balances.
    // threads_per_chain set at runtime (see pioneer_targets.R n_shards).
    target += reduce_sum(
      laplace_partial_sum,
      laplace_patient_range,
      1,  // grainsize
      normalized_psa,
      patient_visit_pos, n_patient_screening_visits, t_patient_visit_idx,
      background_patient_idx,
      lap_tr_offset, lap_frac_offset, lap_init_offset,
      tr_sd_level_intercept[n_levels],
      frac_sd_level_intercept[n_levels],
      init_sd_level_intercept[n_levels],
      measure_sd_psa,
      lap_base_log_hazard_01, lap_base_log_hazard_02,
      enable_ms_01 && enable_ms_pop_time_varying_cov
        ? time_varying_coef_01 : rep_vector(0, 0),
      enable_ms_02 && enable_ms_pop_time_varying_cov
        ? time_varying_coef_02 : rep_vector(0, 0),
      enable_ms_pop_time_varying_cov ? n_time_varying_covar : 0,
      log_baseline_psa,
      median_log_psa_obs, iqr_log_psa_obs,
      ms_final_state, ms_time_01, ms_time_02,
      ms_censored_01, ms_prog_deterministic,
      enable_ms_01, enable_ms_02,
      laplace_newton_tol, laplace_newton_max_iter
    );
  }
}
