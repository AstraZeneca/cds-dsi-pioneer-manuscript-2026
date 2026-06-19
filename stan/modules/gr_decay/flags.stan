// gr_decay/flags.stan
// Gompertz growth-rate decay module — gating flags (optional module, propensity-style).
// enable_gr_decay is the MASTER gate on all gr_decay parameters/transforms/priors.

int<lower=0,upper=1> enable_gr_decay;          // master: warp growth-time when 1
int<lower=0,upper=1> enable_pop_cov_gr_decay;  // pop-level baseline-covariate slopes on log(kappa)

// Per-level intercept modes on log(kappa): index 1..n_levels.
// 0=NONE, 1=FE, 2=RE, 3=RE_GP (rejected here), 4=RE_CP. Mirrors frac.
array[n_levels] int<lower=0,upper=4> enable_level_intercept_gr_decay;
// Per-level covariate-slope flags on log(kappa). Default all 0 (off).
array[n_levels] int<lower=0,upper=1> enable_level_cov_gr_decay;
