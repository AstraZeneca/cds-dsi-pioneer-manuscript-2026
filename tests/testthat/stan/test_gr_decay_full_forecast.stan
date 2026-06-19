// Full-model forecast RNG check (guards BUG 1: full-model forecast ignored kappa).
//
// BUG 1 was that the full-model forecast (generated_quantities.stan, process-noise-OFF
// branch) called the NON-decay RNG with growth factor==1, so kappa never reached the
// forecast and the growth arm grew at the full rate forever. This test calls the SAME
// function the fix routes through — generate_patient_states_with_means_decay_rng — with a
// baseline-anchored tv_factor, and reads the DETERMINISTIC outputs: the *_mean_* channel
// (noise-free calc_log_burden_mean) and the state matrix, both independent of measure_sd.
// (measure_sd must be > 0 — the function calls normal_rng internally, and sd=0 violates
//  normal_rng's positive-scale constraint and NaN-poisons the whole gq block.)
//
// The R harness runs it at two kappa values and asserts the long-horizon forecast DIFFERS
// (smaller kappa => weaker attenuation => higher burden). A regression to the non-decay RNG
// / factor==1 makes the forecast independent of kappa => the two runs coincide => test fails.
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
  real baseline_week;          // absolute week of baseline (anchor of the warp clock)
  int<lower=2> n_forecast;     // forecast steps beyond the anchor
  real horizon_week;           // final forecast week (>> baseline so the plateau bites)
  real init_decrease;          // log-decrease state at the forecast anchor
  real init_growth;            // log-growth state at the forecast anchor
  real decrease_rate;          // > 0 (unwarped)
  real growth_rate;            // > 0
  real kappa;                  // > 0 decay
  real baseline_obs_value;     // baseline SLD (linear units) for the combine
  real<lower=0> measure_sd;    // > 0 (RNG scale); only the noise-free *_mean_* outputs are read
}

generated quantities {
  // Single observed visit at the baseline; forecast anchor is the baseline state.
  matrix[1, 2] patient_states = [[init_decrease, init_growth]];

  // Forecast time grid WITH anchor at baseline_week (so the anchor IS the baseline here).
  array[n_forecast + 1] real forecast_time;
  forecast_time[1] = baseline_week;
  {
    real step = (horizon_week - baseline_week) / n_forecast;
    for (t in 2:(n_forecast + 1)) forecast_time[t] = baseline_week + (t - 1) * step;
  }

  // Baseline-anchored Gompertz tv_factor (production convention).
  vector[n_forecast + 1] tv_factor;
  for (t in 1:(n_forecast + 1)) {
    if (t == 1) {
      tv_factor[t] = 1.0;
    } else {
      real e_hi = forecast_time[t]     - baseline_week;
      real e_lo = forecast_time[t - 1] - baseline_week;
      real dphi = growth_warp(e_hi, kappa) - growth_warp(e_lo, kappa);
      real dt   = forecast_time[t] - forecast_time[t - 1];
      tv_factor[t] = dt > 0 ? dphi / dt : 1.0;
    }
  }

  // Call the exact function the full-model fix routes through. measure_sd=0 => the
  // *_mean_* outputs are the deterministic forecast (no measurement noise).
  matrix[n_forecast, 2] fc_states;
  vector[1] rep_log_obs, rep_mean_log_obs;
  vector[n_forecast] fc_log_obs, fc_mean_log_obs;
  matrix[0, 2] obs_pn;
  // generate_patient_states_with_means_decay_rng takes log-RATES (exp'd internally) and
  // forecasts from patient_states[last]. Pass the intended rates on the log scale.
  (fc_states, rep_log_obs, rep_mean_log_obs, fc_log_obs, fc_mean_log_obs, obs_pn) =
    generate_patient_states_with_means_decay_rng(
      patient_states,
      forecast_time,
      log(decrease_rate),      // patient_log_decrease_rate
      log(growth_rate),        // patient_log_growth_rate
      baseline_obs_value,
      negative_infinity(),     // static compartment disabled
      tv_factor,
      rep_matrix(0.0, n_forecast, 2),
      measure_sd               // > 0 for normal_rng; only *_mean_* outputs are read below
    );

  // Long-horizon forecast burden (log scale), from the NOISE-FREE mean channel. Smaller
  // kappa => weaker attenuation => higher burden.
  real final_forecast_log_obs = fc_mean_log_obs[n_forecast];
  // Growth-arm state at the horizon for a direct plateau comparison.
  real final_growth_state = fc_states[n_forecast, 2];
  real growth_plateau = init_growth + growth_rate / kappa;
}
