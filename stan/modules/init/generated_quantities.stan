// init/generated_quantities.stan
// Back-transform QR coefficients to original covariate space

vector[enable_pop_cov_init && n_covar > 0 ? n_covar : 0] init_coef_pop;
if (enable_pop_cov_init && n_covar > 0) {
  init_coef_pop = mdivide_left_tri_low(R_covar_design_matrix', init_coef_qr_pop);
}
