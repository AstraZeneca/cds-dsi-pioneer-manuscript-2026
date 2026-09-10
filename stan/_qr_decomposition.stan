// ============================================================================
// QR Decomposition of Covariate Design Matrix
// Shared across all models that use time-invariant covariates.
// Requires: n_patients, n_covar, covar_design_matrix in scope.
// ============================================================================

matrix[n_patients, n_covar] Q_covar_design_matrix =
  rep_matrix(0, n_patients, n_covar);
matrix[n_covar, n_covar] R_covar_design_matrix =
  rep_matrix(0, n_covar, n_covar);

if (n_covar > 0) {
  Q_covar_design_matrix = qr_thin_Q(covar_design_matrix) * sqrt(n_patients - 1);
  R_covar_design_matrix = qr_thin_R(covar_design_matrix) / sqrt(n_patients - 1);
}
