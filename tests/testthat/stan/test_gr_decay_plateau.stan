// Plateau check: the growth arm warped by the exact phi-difference tv_factor must
// approach init_g + growth_rate/kappa as the horizon grows, and NEVER exceed it.
// Exercises sf_log_space_trajectory_ncp_decay (the helper the forecast uses).
//
// SCOPE LIMITATION: this harness builds tv_factor with e = times[t] - times[1], i.e.
// it anchors the warp clock at the FORECAST ORIGIN (times[1]), and the harness sets
// times[1]=0 so here cutoff == baseline == 0. It therefore does NOT distinguish a
// baseline-anchored forecast from a cutoff-anchored one, and would stay green even
// with the re-acceleration bug present. The production code anchors at baseline_week
// (NOT times[1]). Do NOT treat this test as a continuity guard: continuity at a
// cutoff != baseline needs a separate test that sets baseline_week < times[1].
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
  #include "modules/gr_decay/gr_decay.stanfunctions"
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
