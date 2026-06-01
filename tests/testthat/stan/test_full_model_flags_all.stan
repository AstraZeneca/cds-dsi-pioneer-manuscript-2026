functions {
  #include "util.stanfunctions"
  #include "pos.stanfunctions"
  #include "full_model.stanfunctions"
}
data {
  int<lower=1> n_combos;
  array[n_combos] int enable_pop_pn;
  array[n_combos] int enable_patient_pn;
  array[n_combos] int enable_states_grid;
  array[n_combos] int enable_ms_tv_cov;
  array[n_combos] int n_tv_covar;
  array[n_combos] int enable_ms_01;
  array[n_combos] int enable_visit_gated;
  array[n_combos] int enable_visit_gated_latent;
  array[n_combos] int enable_02_tv_cov;
  array[n_combos] int enable_03_tv_cov;
  array[n_combos] int has_inline_tv_covar;
}
generated quantities {
  array[n_combos] int out_any_process_noise;
  array[n_combos] int out_need_states_full_grid;
  array[n_combos] int out_ms_needs_inline_burden;

  for (c in 1:n_combos) {
    int apn;
    int nsg;
    int inl;
    (apn, nsg, inl) = compute_full_model_grid_flags(
      enable_pop_pn[c], enable_patient_pn[c],
      enable_states_grid[c],
      enable_ms_tv_cov[c], n_tv_covar[c],
      enable_ms_01[c], enable_visit_gated[c], enable_visit_gated_latent[c],
      enable_02_tv_cov[c], enable_03_tv_cov[c],
      has_inline_tv_covar[c]
    );
    out_any_process_noise[c]     = apn;
    out_need_states_full_grid[c] = nsg;
    out_ms_needs_inline_burden[c]   = inl;
  }
}
