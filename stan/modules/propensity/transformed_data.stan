// modules/propensity/transformed_data.stan
// Locate the contiguous range of target patients using the level-agnostic
// patient_level_groups matrix. Runs once at initialization.
int propensity_target_start = 1;
int propensity_target_end = n_patients;
int propensity_n_target = n_patients;

// Binary indicator: 1 = target patient (weight=1), 0 = non-target (weight=propensity score)
array[n_patients] int<lower=0, upper=1> propensity_is_target = rep_array(1, n_patients);

// Marginal log-odds of being a target patient: log(N_target / N_nontarget).
// Used to stabilize propensity scores into density ratios.
real propensity_log_marginal_odds = 0.0;

if (enable_propensity_weighting && propensity_split_level > 0) {
  propensity_target_start = 0;
  propensity_target_end = 0;
  int n_found = 0;
  for (i in 1:n_patients) {
    if (patient_level_groups[i, propensity_split_level] == propensity_target_group) {
      if (propensity_target_start == 0) propensity_target_start = i;
      propensity_target_end = i;
      n_found += 1;
    }
  }

  propensity_n_target = n_found;

  // Validate that target patients are contiguous (required for range assignment)
  if (propensity_target_end - propensity_target_start + 1 != n_found)
    fatal_error("Propensity target group patients are not contiguous: ",
                "found ", n_found, " patients but range [",
                propensity_target_start, ", ", propensity_target_end, "] has ",
                propensity_target_end - propensity_target_start + 1, " slots");

  propensity_is_target = zeros_int_array(n_patients);
  for (i in propensity_target_start:propensity_target_end) {
    propensity_is_target[i] = 1;
  }

  // Marginal log-odds: log(N_target / N_nontarget)
  // Stabilizes the propensity logit into a density ratio:
  //   log(P(X|trial)/P(X|RWD)) = logit(P(trial|X)) - logit(P(trial))
  propensity_log_marginal_odds = log(propensity_n_target * 1.0)
                                 - log((n_patients - propensity_n_target) * 1.0);

  print("Propensity weighting enabled:");
  print("  split_level=", propensity_split_level,
        " target_group=", propensity_target_group);
  print("  target: [", propensity_target_start, ", ", propensity_target_end,
        "] n=", propensity_n_target,
        " non-target: ", n_patients - propensity_n_target);
  print("  log_marginal_odds=", propensity_log_marginal_odds);
}
