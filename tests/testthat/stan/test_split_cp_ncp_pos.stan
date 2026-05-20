functions {
  #include "util.stanfunctions"
  #include "pos.stanfunctions"
  #include "hierarchy.stanfunctions"
}

data {
  int<lower=1> n_levels;
  array[n_levels] int<lower=0,upper=4> mode;
  array[n_levels] int<lower=0> n_groups;
}

generated quantities {
  int n_raw;
  int n_cp;
  array[n_levels + 1] int raw_pos;
  array[n_levels + 1] int cp_pos;
  (n_raw, raw_pos, n_cp, cp_pos) = split_cp_ncp_pos(n_levels, n_groups, mode);
}
