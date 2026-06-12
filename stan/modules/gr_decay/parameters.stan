// gr_decay/parameters.stan
// Lean pop-level parameters for log(kappa). Master-flag gated: sized 0 when off
// (an off model samples/stores nothing), mirroring propensity/parameters.stan.

// Population intercept on log scale => kappa = exp(.) > 0. Array size 1 when on, 0 when off
// (mirrors init_logit_static_loc_pop).
array[enable_gr_decay ? 1 : 0] real gr_decay_log_loc_pop;

// Population covariate coefficients (QR space). Inner gate: 0 unless both the master
// flag and the covariate flag are on.
vector[(enable_gr_decay && enable_pop_cov_gr_decay) ? n_covar : 0] gr_decay_coef_qr_pop;
