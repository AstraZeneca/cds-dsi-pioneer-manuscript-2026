// shared/helpers.stan
// Inline utility functions & common transforms (to be populated as refactor progresses)

/**
 * Stable logit-derived logs
 */
real logit_to_log_prob(real x_logit) {
  return -log1p_exp(-x_logit);
}
real logit_to_log_comp(real x_logit) {
  return -log1p_exp(x_logit);
}
