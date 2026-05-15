// modules/multistate/cond_surv_transform.stan
// Apply the -exp() transform to log_cond_surv_01 and log_cond_surv_02.
//
// Separated from transformed_parameters.stan so that model-specific includes
// (e.g., inline PSA TV covariates) can add to log_cond_surv before the
// transform is applied. Include this AFTER any such model-specific code.

if (enable_ms_01) {
  log_cond_surv_01 = -exp(log_cond_surv_01);
}

if (enable_ms_02) {
  log_cond_surv_02 = -exp(log_cond_surv_02);
}
