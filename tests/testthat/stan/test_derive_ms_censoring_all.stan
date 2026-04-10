functions {
  #include "pfs.stanfunctions"
  #include "util.stanfunctions"
  #include "pos.stanfunctions"
  #include "multistate.stanfunctions"
}
data {
  int<lower=1> N;
  array[N] int ms_final_state;
  array[N] int ms_time_01;
  array[N] int ms_time_32;
}
generated quantities {
  array[N] int censored_02;
  array[N] int censored_12;
  array[N] int censored_32;
  (censored_02, censored_12, censored_32) = derive_ms_censoring_indicators(
    ms_final_state, ms_time_01, ms_time_32
  );
}
