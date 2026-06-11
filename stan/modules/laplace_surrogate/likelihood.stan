// --- Laplace surrogate: backgrounded-patient marginalized SLD ---------------
// Backgrounded patients inform population params only; their per-patient latents
// are integrated out via the validated quadratic surrogate. Forecast patients
// (and the current default behavior) are untouched.
if (enable_background_surrogate == 1 && n_background_patients > 0) {
  // Patient-level total marginal SD = sqrt(sum over enabled levels of sd^2).
  // (frac/init always; tr guaranteed simple here by the transformed_data guard.)
  real tr_patient_sd = 0;
  real frac_patient_sd = 0;
  real init_patient_sd = 0;
  for (lv in 1:n_levels) {
    if (enable_level_intercept_tr[lv]   != 0) tr_patient_sd   += square(tr_sd_level_intercept[lv]);
    if (enable_level_intercept_frac[lv] != 0) frac_patient_sd += square(frac_sd_level_intercept[lv]);
    if (enable_level_intercept_init[lv] != 0) init_patient_sd += square(init_sd_level_intercept[lv]);
  }
  tr_patient_sd   = sqrt(tr_patient_sd);
  frac_patient_sd = sqrt(frac_patient_sd);
  init_patient_sd = sqrt(init_patient_sd);

  // Compact background-only views (surrogate_bg_obs / surrogate_bg_time /
  // surrogate_bg_pos / surrogate_bg_log_lod) are assembled in the module's
  // transformed_data.stan so they satisfy surrogate_ll's data-only qualifiers.

  vector[3] bg_beta_pop;
  matrix[2, 2] bg_Sigma;
  (bg_beta_pop, bg_Sigma) = surrogate_bridge(
    tr_loc_pop, frac_logit_loc_pop, init_logit_loc_pop,
    tr_patient_sd, frac_patient_sd, init_patient_sd,
    surrogate_Vinv, surrogate_anchor_times,
    surrogate_gh_x, surrogate_gh_w, surrogate_jitter);

  if (surrogate_n_frailty_slots > 0) {
    // --- per-bg static baselines (theta-independent), row-major flattened ---
    vector[n_background_patients * surrogate_n_wk] base01_static_flat;
    vector[n_background_patients * surrogate_n_wk] base03_static_flat;
    vector[n_background_patients] log_baseline_burden_bg;
    // TI-cov linear predictors in QR space on the BACKGROUND rows (mirror
    // transformed_parameters.stan:216,787). Publication: pop-level term only
    // (enable_ms_level_cov both FALSE), so no level-slope path. GUARD with the
    // SAME conditions the forecast path uses — time_invariant_coef_qr_01/03 are
    // vector[0] when their flags are off (parameters.stan:47,154), so an
    // unconditional multiply would be a (n_bg × n_covar)·vector[0] mismatch.
    vector[n_background_patients] linpred_bg_01 = zeros_vector(n_background_patients);
    vector[n_background_patients] linpred_bg_03 = zeros_vector(n_background_patients);
    if (enable_ms_pop_time_invariant_cov && n_time_invariant_covar > 0)
      linpred_bg_01 = Q_covar_design_matrix[background_patient_idx, ] * time_invariant_coef_qr_01;
    if (enable_ms_pop_time_invariant_cov && enable_ms_03_time_invariant_cov && n_time_invariant_covar > 0)
      linpred_bg_03 = Q_covar_design_matrix[background_patient_idx, ] * time_invariant_coef_qr_03;
    for (j in 1:n_background_patients) {
      int p = background_patient_idx[j];
      log_baseline_burden_bg[j] = log_baseline_burden[p];
      row_vector[surrogate_n_wk] b01 = log_pop_lambda_01;
      row_vector[surrogate_n_wk] b03 = log_pop_lambda_03;
      for (lv in 1:n_levels) {
        if (lv < n_levels) {  // trial-level baseline residual (patient level is frailty, in theta)
          if (enable_ms_01)
            b01 += log_level_lambda_01_residual[patient_ms_baseline_flat_idx_slot[MS_SLOT_01, p, lv]];
          if (enable_ms_03)
            b03 += log_level_lambda_03_residual[patient_ms_baseline_flat_idx_slot[MS_SLOT_03, p, lv]];
        }
      }
      b01 += linpred_bg_01[j];
      b03 += linpred_bg_03[j];
      for (w in 1:surrogate_n_wk) {
        base01_static_flat[(j - 1) * surrogate_n_wk + w] = b01[w];
        base03_static_flat[(j - 1) * surrogate_n_wk + w] = b03[w];
      }
    }

    // --- Sigma_u: diag(sigma) * (L L') * diag(sigma) (matches ms_corr_u, tp:37) ---
    vector[surrogate_n_frailty_slots] sigma_u;
    sigma_u[1] = (surrogate_frailty_slot[1] == MS_SLOT_01)
               ? log_lambda_gp_01_level_intercept_sd[surrogate_patient_lv]
               : log_lambda_gp_03_level_intercept_sd[surrogate_patient_lv];
    if (surrogate_n_frailty_slots >= 2)
      sigma_u[2] = log_lambda_gp_03_level_intercept_sd[surrogate_patient_lv];
    matrix[surrogate_n_frailty_slots, surrogate_n_frailty_slots] Lu =
      diag_pre_multiply(sigma_u, L_ms_intercept_corr[surrogate_frailty_block]);
    matrix[surrogate_n_frailty_slots, surrogate_n_frailty_slots] Sigma_u =
      Lu * Lu' + diag_matrix(rep_vector(surrogate_jitter, surrogate_n_frailty_slots));

    target += laplace_marginal_tol(
      surrogate_ll,
      (bg_beta_pop, measure_sd_sld, surrogate_bg_log_lod, n_background_patients,
       surrogate_bg_obs, surrogate_bg_pos, surrogate_bg_time,
       surrogate_d, surrogate_n_frailty_slots, surrogate_n_wk,
       base01_static_flat, base03_static_flat, log_baseline_burden_bg,
       time_varying_coef_01, time_varying_coef_03,
       median_log_burden_obs, iqr_log_burden_obs, median_velocity_obs, iqr_velocity_obs,
       surrogate_bg_event_wk_01, surrogate_bg_censored_01,
       surrogate_bg_event_wk_03, surrogate_bg_censored_03,
       surrogate_bg_baseline_week, surrogate_bg_visit_wk_01),
      surrogate_hessian_block_size,
      surrogate_K_fn,
      (bg_Sigma, Sigma_u, n_background_patients, surrogate_d),
      (surrogate_theta_0, surrogate_tolerance, surrogate_max_num_steps,
       surrogate_solver, surrogate_max_steps_line_search, surrogate_allow_fallback)
    );
  } else {
    target += laplace_marginal_tol(
      surrogate_ll,
      (bg_beta_pop, measure_sd_sld, surrogate_bg_log_lod, n_background_patients,
       surrogate_bg_obs, surrogate_bg_pos, surrogate_bg_time),
      surrogate_hessian_block_size,
      surrogate_K_fn,
      (bg_Sigma, n_background_patients),
      (surrogate_theta_0, surrogate_tolerance, surrogate_max_num_steps,
       surrogate_solver, surrogate_max_steps_line_search, surrogate_allow_fallback)
    );
  }
}
