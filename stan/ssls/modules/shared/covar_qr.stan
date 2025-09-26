// shared/covar_qr.stan
// Placeholder for shared covariate design matrices and QR components.
/*
int<lower=0> K_tr;     // number of total rate covariates
int<lower=0> K_frac;   // number of fraction covariates
int<lower=0> K_init;   // number of initial proportion covariates

matrix[n_train_patients, K_tr]   Q_tr;
matrix[n_train_patients, K_frac] Q_frac;
matrix[n_train_patients, K_init] Q_init;
// Upper-triangular R matrices (if distinct per module); else share one R if covariate sets identical.
matrix[K_tr, K_tr]   R_tr;
matrix[K_frac, K_frac] R_frac;
matrix[K_init, K_init] R_init;
*/
