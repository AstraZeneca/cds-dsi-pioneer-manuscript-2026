// init/generated_quantities.stan
// Back-transform QR coefficients to original covariate space

vector[enable_pop_cov_init && n_covar > 0 ? n_covar : 0] init_coef_pop;
if (enable_pop_cov_init && n_covar > 0) {
  init_coef_pop = mdivide_left_tri_low(R_covar_design_matrix', init_coef_qr_pop);
}

// Population-level baseline composition (diagnostics).
// pi_decrease = inv_logit(loc_pop); rest = 1 - pi_decrease.
real init_pi_decrease_pop = inv_logit(init_logit_loc_pop);
real init_pi_static_pop = enable_static_init
  ? (1 - init_pi_decrease_pop) * inv_logit(init_logit_static_loc_pop[1])
  : 0.0;
real init_pi_growth_pop = (1 - init_pi_decrease_pop) - init_pi_static_pop;
