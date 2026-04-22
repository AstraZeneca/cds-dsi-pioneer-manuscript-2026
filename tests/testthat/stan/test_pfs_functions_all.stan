functions {
  #include "multistate.stanfunctions"
  #include "pfs.stanfunctions"
  #include "multistate.stanfunctions"
  #include "util.stanfunctions"
  #include "pos.stanfunctions"
}

data {
  int<lower=1> N_CASES;
  int<lower=1> N;
  array[N_CASES, N] int pfs;
  array[N_CASES, N] int right_censored;
  array[N_CASES, N] int max_time;
}

generated quantities {
  array[N_CASES, N] int out_pfs;
  array[N_CASES, N] int out_right_censored;

  for (case in 1:N_CASES) {
    tuple(array[N] int, array[N] int) res =
      truncate_at_max_time(pfs[case], right_censored[case], max_time[case]);
    out_pfs[case]            = res.1;
    out_right_censored[case] = res.2;
  }
}
