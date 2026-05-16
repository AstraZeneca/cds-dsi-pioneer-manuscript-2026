functions {
  #include "util.stanfunctions"
  #include "pos.stanfunctions"
  #include "hierarchy.stanfunctions"
}

data {
  int<lower=1> n_levels;
  array[n_levels] int<lower=0> n_groups_per_level;
  array[n_levels, n_levels] int<lower=0,upper=4> sd_mode;
}

generated quantities {
  int raw_total;
  array[n_levels, n_levels + 1] int raw_pos;
  int cp_total;
  array[n_levels, n_levels + 1] int cp_pos;
  (raw_total, raw_pos, cp_total, cp_pos) =
    split_sd_cp_ncp_pos(n_levels, n_groups_per_level, sd_mode);
}
