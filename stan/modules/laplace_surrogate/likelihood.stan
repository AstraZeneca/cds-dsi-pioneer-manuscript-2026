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
