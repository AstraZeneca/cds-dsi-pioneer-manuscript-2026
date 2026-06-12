// gr_decay/priors.stan — flag-gated weakly-informative prior on log(kappa).
// Entire block omitted when off (zero sampling cost).

if (enable_gr_decay) {
  gr_decay_log_loc_pop[1] ~ normal(gr_decay_log_loc_pop_mean, gr_decay_log_loc_pop_sd);
  if (enable_pop_cov_gr_decay) {
    gr_decay_coef_qr_pop ~ normal(gr_decay_coef_qr_pop_mean, gr_decay_coef_qr_pop_sd);
  }
}
