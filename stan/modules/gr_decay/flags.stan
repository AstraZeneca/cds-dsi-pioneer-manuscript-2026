// gr_decay/flags.stan
// Gompertz growth-rate decay module — gating flags (optional module, propensity-style).
// enable_gr_decay is the MASTER gate on all gr_decay parameters/transforms/priors.

int<lower=0,upper=1> enable_gr_decay;          // master: warp growth-time when 1
int<lower=0,upper=1> enable_pop_cov_gr_decay;  // pop-level baseline-covariate slopes on log(kappa)
