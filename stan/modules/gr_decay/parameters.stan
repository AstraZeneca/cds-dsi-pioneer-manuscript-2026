// gr_decay/parameters.stan
// Lean pop-level parameters for log(kappa). Master-flag gated: sized 0 when off
// (an off model samples/stores nothing), mirroring propensity/parameters.stan.

// Population intercept on log scale => kappa = exp(.) > 0. Array size 1 when on, 0 when off
// (mirrors init_logit_static_loc_pop).
array[enable_gr_decay ? 1 : 0] real gr_decay_log_loc_pop;

// Population covariate coefficients (QR space). Inner gate: 0 unless both the master
// flag and the covariate flag are on.
vector[(enable_gr_decay && enable_pop_cov_gr_decay) ? n_covar : 0] gr_decay_coef_qr_pop;

// ===== UNIFIED LEVEL STRUCTURE (mirrors frac) =====
// All sized by enabled-group counts from transformed_data; collapse to 0 when
// every level mode is NONE (pop-only path) or when enable_gr_decay is off
// (the R side passes all-NONE flag arrays when off — see Task 9 wiring — so the
// counts are 0 and these params vanish, preserving "off = zero params").

// Intercept SD free params — one per RE/RE_CP level.
array[n_re_levels_gr_decay_intercept] real<lower=0> gr_decay_sd_level_intercept_raw;

// Raw NCP and CP intercept draws — split by mode.
vector[n_raw_groups_gr_decay_intercept] gr_decay_raw_level_intercept;
vector[n_cp_groups_gr_decay_intercept]  gr_decay_cp_level_intercept;

// Slope SD hyperparameters — one vector per level.
array[n_levels] vector<lower=0>[n_covar] gr_decay_sd_level_slope;

// Raw NCP and CP slope draws — split by mode.
matrix[n_raw_groups_gr_decay_slope, n_covar] gr_decay_raw_level_slope;
matrix[n_cp_groups_gr_decay_slope,  n_covar] gr_decay_cp_level_slope;
