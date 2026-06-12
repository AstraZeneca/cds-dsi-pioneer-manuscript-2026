// Plateau check for the LFO forecast path: the growth arm warped by the exact
// phi-difference tv_factor must approach init_g + growth_rate/kappa as the horizon
// grows, and NEVER exceed it. Exercises sf_log_space_trajectory_ncp_decay — the
// SAME helper the LFO forecast call uses (sf-ssls-lfo.stan), with tv_factor built
// identically (exact phi-difference per step).
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
  int<lower=2> n_t;          // number of forecast time points (incl. anchor)
  array[n_t] real times;     // forecast weeks (times[1] = anchor / forecast origin)
  real init_decrease;        // anchor log-decrease state
  real init_growth;          // anchor log-growth state
  real decrease_rate;        // > 0
  real growth_rate;          // > 0
  real kappa;                // > 0 decay
}

generated quantities {
  // Build the exact phi-difference per-step factor (mirror of sf-ssls-lfo.stan:205 path).
  vector[n_t] tv_factor;
  for (t in 1:n_t) {
    if (t == 1) {
      tv_factor[t] = 1.0;
    } else {
      real e_hi = times[t] - times[1];
      real e_lo = times[t - 1] - times[1];
      real dphi = growth_warp(e_hi, kappa) - growth_warp(e_lo, kappa);
      real dt   = times[t] - times[t - 1];
      tv_factor[t] = dt > 0 ? dphi / dt : 1.0;
    }
  }

  row_vector[2] x0 = [init_decrease, init_growth];
  matrix[n_t, 2] expected_x;
  matrix[n_t, 2] traj;
  (expected_x, traj) = sf_log_space_trajectory_ncp_decay(
    x0, times, decrease_rate, growth_rate, tv_factor,
    rep_matrix(0.0, n_t - 1, 2), 0
  );

  // The growth arm (column 2). Its asymptote is init_growth + growth_rate/kappa.
  vector[n_t] growth_arm = traj[, 2];
  real plateau = init_growth + growth_rate / kappa;
  real final_growth = growth_arm[n_t];
  // Max over the whole trajectory — must never exceed the plateau (monotone bounded).
  real max_growth = max(growth_arm);
}
