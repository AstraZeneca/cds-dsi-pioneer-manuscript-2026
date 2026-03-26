// modules/propensity/flags.stan
// Feature flags for propensity-weighted borrowing.
// propensity_split_level = 0 disables the module entirely.
int<lower=0, upper=1> enable_propensity_weighting;
int<lower=0> propensity_split_level;     // >0: which hierarchy level to split on
int<lower=0> propensity_target_group;    // group ID at that level (weight = 1.0)
int<lower=0, upper=1> fit_propensity_data;  // 0 during prior predictive runs
