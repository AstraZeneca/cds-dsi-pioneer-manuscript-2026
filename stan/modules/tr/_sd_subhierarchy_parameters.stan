// tr/_sd_subhierarchy_parameters.stan
// Parameters for the SD sub-hierarchy (issue #110).
// All sizes are zero when the mode matrix is all-NONE → no new params allocated → bit-exact.

// Population log-SD per location level (one entry per level with an active sub-hierarchy)
vector[n_subhier_active_tr_intercept] tr_log_sd_level_intercept_pop;

// Free hyperscales for RE / RE_CP sub-levels
vector<lower=0>[n_sd_hyperscales_tr_intercept] tr_sd_hyperscale_level_intercept_raw;

// Flat raw bucket for FE/RE sub-levels (indexed by raw_pos_tr_log_sd_intercept)
vector[n_raw_groups_tr_log_sd_intercept] tr_raw_log_sd_level_intercept;

// Flat cp bucket for RE_CP sub-levels (indexed by cp_pos_tr_log_sd_intercept)
vector[n_cp_groups_tr_log_sd_intercept] tr_cp_log_sd_level_intercept;
