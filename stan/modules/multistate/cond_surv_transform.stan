// modules/multistate/cond_surv_transform.stan
// Apply the -exp() transform to log_cond_surv_01 and log_cond_surv_02.
//
// Separated from transformed_parameters.stan so that model-specific includes
// (e.g., inline PSA TV covariates) can add to log_cond_surv before the
// transform is applied. Include this AFTER any such model-specific code.
//
// Clamping: log-hazard ∈ [-20, 10] before exp().
//   upper 10: exp(10)≈22000/week already means instantaneous death
//   lower -20: prevents exp() underflow to 0 → log_cond_surv = 0 → log1m_exp(0) = -Inf

if (enable_ms_01) {
  log_cond_surv_01 = -exp(fmax(log_cond_surv_01, -20.0));
}

if (enable_ms_02) {
  log_cond_surv_02 = -exp(fmax(log_cond_surv_02, -20.0));
}
