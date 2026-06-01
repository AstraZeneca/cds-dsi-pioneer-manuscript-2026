// modules/multistate/cond_surv_transform.stan
// Apply the -exp() transform to all from-state-0 log conditional survival matrices.
//
// Separated from transformed_parameters.stan so that model-specific includes
// (e.g., inline burden TV covariates) can add to log_cond_surv before the
// transform is applied. Include this AFTER _ms_burden_inline_tv_covar.stan.

if (enable_ms_01) {
  log_cond_surv_01 = -exp(log_cond_surv_01);
}

if (enable_ms_02) {
  log_cond_surv_02 = -exp(log_cond_surv_02);
}

if (enable_ms_03) {
  log_cond_surv_03 = -exp(log_cond_surv_03);
}
