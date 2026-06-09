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

  // Compact background-only views; per-patient LOD offset (burden normalized to
  // each patient's baseline => LOD shifts by log_baseline_sld[p]).
  array[n_background_patients + 1] int bg_pos;
  bg_pos[1] = 1;
  for (j in 1:n_background_patients) {
    int p = background_patient_idx[j];
    int vs, ve;
    (vs, ve) = get_pos(patient_visit_pos, p);
    bg_pos[j + 1] = bg_pos[j] + (ve - vs + 1);
  }
  int n_bg_visits = bg_pos[n_background_patients + 1] - 1;
  vector[n_bg_visits] bg_obs;
  array[n_bg_visits] int bg_time;
  // Per-patient LOD offset replicated per visit (surrogate_ll uses one global
  // log_lod, so we pass each patient's already-offset value via a vector form).
  vector[n_bg_visits] bg_log_lod;
  {
    int w = 1;
    for (j in 1:n_background_patients) {
      int p = background_patient_idx[j];
      int vs, ve;
      (vs, ve) = get_pos(patient_visit_pos, p);
      for (v in vs:ve) {
        bg_obs[w]     = normalized_sld[v];
        bg_time[w]    = t_patient_visit_idx[v];
        bg_log_lod[w] = log_lod - log_baseline_sld[p];
        w += 1;
      }
    }
  }

  vector[3] bg_beta_pop;
  matrix[2, 2] bg_Sigma;
  (bg_beta_pop, bg_Sigma) = surrogate_bridge(
    tr_loc_pop, frac_logit_loc_pop, init_logit_loc_pop,
    tr_patient_sd, frac_patient_sd, init_patient_sd,
    surrogate_Vinv, surrogate_anchor_times,
    surrogate_gh_x, surrogate_gh_w, surrogate_jitter);

  target += laplace_marginal_tol(
    surrogate_ll,
    (bg_beta_pop, measure_sd_sld, bg_log_lod, n_background_patients,
     bg_obs, bg_pos, bg_time),
    surrogate_hessian_block_size,
    surrogate_K_fn,
    (bg_Sigma, n_background_patients),
    (surrogate_theta_0, surrogate_tolerance, surrogate_max_num_steps,
     surrogate_solver, surrogate_max_steps_line_search, surrogate_allow_fallback)
  );
}
