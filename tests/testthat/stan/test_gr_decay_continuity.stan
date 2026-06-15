// Continuity-at-cutoff check for the warped forecast (guards BUG 2: cutoff re-anchoring).
//
// The forecast must CONTINUE the single warp clock anchored at the BASELINE week, so the
// growth arm is the exact analytic continuation of the in-sample arm:
//   growth_arm(week) = init_g + growth_rate * phi(week - baseline_week, kappa)   for ALL weeks,
// in-sample AND forecast. This test sets baseline_week < cutoff_week < horizon (all distinct)
// so it can DISTINGUISH a baseline-anchored forecast from a cutoff-anchored one — the exact
// discrimination test_gr_decay_plateau.stan cannot make (it pins times[1]=baseline=0).
//
// It builds the forecast tv_factor two ways:
//   correct_*  : warp anchored at baseline_week (the production convention, post-fix)
//   buggy_*    : warp anchored at the cutoff (forecast_time[1]) — the BUG 2 anchoring
// and exposes both so the R harness can assert (a) correct == analytic continuation, and
// (b) buggy != correct (proving the test would FAIL if the bug were reintroduced).
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
  real baseline_week;          // absolute week of the patient's baseline (last screening visit)
  real cutoff_week;            // absolute week of the last observed visit (> baseline_week)
  int<lower=2> n_forecast;     // number of forecast steps beyond the cutoff
  real horizon_week;           // absolute week of the final forecast point (> cutoff_week)
  real init_decrease;          // log-decrease state at the BASELINE
  real init_growth;            // log-growth state at the BASELINE
  real decrease_rate;          // > 0 (unwarped)
  real growth_rate;            // > 0
  real kappa;                  // > 0 decay
}

generated quantities {
  // ---- Forecast time grid (with anchor at the cutoff, mirroring production) ----
  // forecast_time[1] = cutoff_week (anchor, not advanced); then n_forecast points to horizon.
  array[n_forecast + 1] real forecast_time;
  forecast_time[1] = cutoff_week;
  {
    real step = (horizon_week - cutoff_week) / n_forecast;
    for (t in 2:(n_forecast + 1)) forecast_time[t] = cutoff_week + (t - 1) * step;
  }

  // ---- In-sample anchor state at the CUTOFF, on the baseline clock ----
  // This is what transformed_parameters produces for the last observed visit:
  //   growth arm = init_growth + growth_rate * phi(cutoff - baseline).
  real cutoff_elapsed = cutoff_week - baseline_week;
  real x0_decrease = init_decrease - decrease_rate * cutoff_elapsed;   // decrease arm is unwarped
  real x0_growth   = init_growth + growth_rate * growth_warp(cutoff_elapsed, kappa);
  row_vector[2] x0 = [x0_decrease, x0_growth];

  // ---- CORRECT factor: warp anchored at baseline_week (production, post-fix) ----
  vector[n_forecast + 1] tv_correct;
  for (t in 1:(n_forecast + 1)) {
    if (t == 1) {
      tv_correct[t] = 1.0;
    } else {
      real e_hi = forecast_time[t]     - baseline_week;
      real e_lo = forecast_time[t - 1] - baseline_week;
      real dphi = growth_warp(e_hi, kappa) - growth_warp(e_lo, kappa);
      real dt   = forecast_time[t] - forecast_time[t - 1];
      tv_correct[t] = dt > 0 ? dphi / dt : 1.0;
    }
  }

  // ---- BUGGY factor: warp anchored at the cutoff (forecast_time[1]) ----
  vector[n_forecast + 1] tv_buggy;
  for (t in 1:(n_forecast + 1)) {
    if (t == 1) {
      tv_buggy[t] = 1.0;
    } else {
      real e_hi = forecast_time[t]     - forecast_time[1];
      real e_lo = forecast_time[t - 1] - forecast_time[1];
      real dphi = growth_warp(e_hi, kappa) - growth_warp(e_lo, kappa);
      real dt   = forecast_time[t] - forecast_time[t - 1];
      tv_buggy[t] = dt > 0 ? dphi / dt : 1.0;
    }
  }

  // ---- Run both forecasts through the SAME trajectory helper the model uses ----
  matrix[n_forecast + 1, 2] exp_c, traj_c, exp_b, traj_b;
  (exp_c, traj_c) = sf_log_space_trajectory_ncp_decay(
    x0, forecast_time, decrease_rate, growth_rate, tv_correct,
    rep_matrix(0.0, n_forecast, 2), 0
  );
  (exp_b, traj_b) = sf_log_space_trajectory_ncp_decay(
    x0, forecast_time, decrease_rate, growth_rate, tv_buggy,
    rep_matrix(0.0, n_forecast, 2), 0
  );

  // ---- Analytic continuation on the one baseline clock ----
  // growth_arm(week) = init_growth + growth_rate * phi(week - baseline).
  vector[n_forecast + 1] analytic_growth;
  for (t in 1:(n_forecast + 1)) {
    analytic_growth[t] = init_growth + growth_rate * growth_warp(forecast_time[t] - baseline_week, kappa);
  }

  // ---- Exposed scalars for the R harness ----
  // Seam: first forecast point (the anchor) must equal the in-sample cutoff state exactly.
  real seam_growth_correct = traj_c[1, 2];
  real cutoff_state_growth = x0_growth;

  // Correct forecast growth arm must match the analytic continuation at every step.
  real max_abs_err_correct = max(abs(traj_c[, 2] - analytic_growth));

  // The buggy (cutoff-anchored) forecast must DIFFER from the correct one — proves the
  // test discriminates the bug. Measured at the final (longest-lever) forecast point.
  real final_growth_correct = traj_c[n_forecast + 1, 2];
  real final_growth_buggy   = traj_b[n_forecast + 1, 2];
  real final_growth_analytic = analytic_growth[n_forecast + 1];

  // The buggy variant re-accelerates => its growth arm is STRICTLY HIGHER than correct.
  real buggy_minus_correct_final = final_growth_buggy - final_growth_correct;
}
