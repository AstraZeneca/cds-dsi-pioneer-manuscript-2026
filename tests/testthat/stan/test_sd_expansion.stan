// Test harness for SD expansion logic (transformed_parameters pattern).
// Verifies that the SD array assembled from RE free params and FE hyperparams
// matches expectations for all mode combinations (none=0, fe=1, re=2).
// No functions needed — logic is self-contained.

data {
  int<lower=1> n_levels;
  // Mode per level: 0=none, 1=fe, 2=re
  array[n_levels] int<lower=0,upper=2> mode;
  // FE SD hyperparams (only used for mode=1 levels)
  array[n_levels] real<lower=0> fe_sd;
  // Number of RE levels (for sizing the free SD param array)
  int<lower=0> n_re_levels;
  // Free RE SD values (one per RE level, in level order)
  array[n_re_levels] real<lower=0> re_sd_values;
}

generated quantities {
  // Expand to full n_levels SD array — mirrors transformed_parameters logic
  array[n_levels] real<lower=0> sd_expanded;
  {
    int sd_idx = 0;
    for (lv in 1:n_levels) {
      if (mode[lv] == 1) {  // LEVEL_MODE_FE
        sd_expanded[lv] = fe_sd[lv];
      } else if (mode[lv] == 2) {  // LEVEL_MODE_RE
        sd_idx += 1;
        sd_expanded[lv] = re_sd_values[sd_idx];
      } else {
        sd_expanded[lv] = 0.0;
      }
    }
  }
}
