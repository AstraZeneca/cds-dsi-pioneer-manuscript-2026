functions {
  #include "util.stanfunctions"
  #include "pos.stanfunctions"
  #include "hierarchy.stanfunctions"
  #include "full_model.stanfunctions"
  #include "gp.stanfunctions"
  #include "pfs.stanfunctions"
  #include "lfo.stanfunctions"
  #include "multistate.stanfunctions"
  #include "_burden.stanfunctions"
  #include "modules/state_space/sf.stanfunctions"
}

data {
  int<lower=1> n_t;
  vector[n_t] t;          // elapsed times (weeks)
  real kappa;             // decay rate
}

generated quantities {
  // Scalar warp at each time
  vector[n_t] phi;
  for (k in 1:n_t) phi[k] = growth_warp(t[k], kappa);

  // Row-vector overload must agree with scalar elementwise
  row_vector[n_t] phi_rv = growth_warp(to_row_vector(t), kappa);

  // Telescoping identity: sum of phi-differences from 0 == phi(t_n)
  // (build differences against a prepended 0)
  real phi_last = growth_warp(t[n_t], kappa);
  real telescoped = 0;
  {
    real prev = 0;  // phi(0) = 0
    for (k in 1:n_t) {
      real cur = growth_warp(t[k], kappa);
      telescoped += (cur - prev);
      prev = cur;
    }
  }
}
