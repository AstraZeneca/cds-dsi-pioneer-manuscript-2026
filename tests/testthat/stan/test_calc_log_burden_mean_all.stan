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
  int<lower=1> n_visits;
  matrix[n_visits, 2] patient_states;  // log-space states (decrease, growth)
  real baseline;
  real static_log_level;               // log(pi_static); -inf disables
}

generated quantities {
  // Two-arg overload (must equal current behavior)
  vector[n_visits] out_2arg = calc_log_burden_mean(patient_states, baseline);
  // Three-arg overload with the supplied static level
  vector[n_visits] out_3arg = calc_log_burden_mean(patient_states, baseline, static_log_level);
}
