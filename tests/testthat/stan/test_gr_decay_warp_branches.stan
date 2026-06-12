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
  int<lower=1> n_t;        // grid width
  real growth_rate;        // constant rate
  real kappa;
  real init_g;
}

generated quantities {
  row_vector[n_t] elapsed = linspaced_row_vector(n_t, 0, n_t - 1);

  // Direct form (B3/B4): init + rate * phi(elapsed)
  row_vector[n_t] direct;
  for (k in 1:n_t) direct[k] = init_g + growth_rate * growth_warp(elapsed[k], kappa);

  // Cumsum form (B1/B2): init + cumsum(rate * (phi(t_k)-phi(t_{k-1})))
  row_vector[n_t] cumsum_form;
  cumsum_form[1] = init_g;
  {
    row_vector[n_t] warp = growth_warp(elapsed, kappa);
    row_vector[n_t-1] warp_diff = warp[2:] - warp[:(n_t-1)];
    cumsum_form[2:] = init_g + cumulative_sum(rep_row_vector(growth_rate, n_t-1) .* warp_diff);
  }
}
