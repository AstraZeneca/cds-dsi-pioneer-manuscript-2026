// Minimal harness to test the share_dead_gp_shape validation block in
// modules/multistate/transformed_data.stan.
//
// The real transformed_data.stan calls compute_ms_time_scale_flags then checks
// three fatal_error guards. This harness replicates that logic so we can test
// all three invalid configurations produce runtime errors, and that the valid
// configuration passes silently.
functions {
  #include "util.stanfunctions"
  #include "pos.stanfunctions"
  #include "pfs.stanfunctions"
  #include "multistate.stanfunctions"
}
data {
  int<lower=0, upper=1> share_dead_gp_shape;
  int<lower=0>          ms_time_scale_12;   // 0=Markov, 1=semi-Markov, 2=extended
  int<lower=0, upper=1> enable_ms_12;
  int<lower=0, upper=1> enable_ms_02;
}
transformed data {
  // Replicates the flag derivation and validation from
  // modules/multistate/transformed_data.stan lines 38-49.
  int is_valid;
  int need_12_s_gp;
  int need_12_t_gp;
  int ms_12_t_has_intercept;
  (is_valid, need_12_s_gp, need_12_t_gp, ms_12_t_has_intercept) =
    compute_ms_time_scale_flags(enable_ms_12, ms_time_scale_12, 1);

  if (share_dead_gp_shape) {
    if (ms_time_scale_12 != 0)
      fatal_error("share_dead_gp_shape requires ms_time_scale_12=0 (Markov); extended mode not supported");
    if (!need_12_t_gp)
      fatal_error("share_dead_gp_shape requires need_12_t_gp=1");
    if (!enable_ms_02)
      fatal_error("share_dead_gp_shape requires enable_ms_02=1");
  }
}
generated quantities {
  int passed = 1;
}
