// Phase 1b: ABLATION harness — isolate which term breaks the Laplace approx.
// ============================================================================
// Same tumor likelihood as laplace_tumor_test.stan, but two data knobs let us
// turn pieces off to see which one destroys log-concavity in z:
//
//   n_latent  (1, 2, 3): how many NCP latents per patient are integrated out.
//     1 -> z_tr only (z_frac, z_init fixed at 0 = population value)
//     2 -> z_tr, z_frac
//     3 -> z_tr, z_frac, z_init (full)
//   The below-LOD censoring (normal_lcdf) is controlled by the DATA: encode an
//   observation as <= 0 to censor it, or keep it positive to not. So a "no-LOD"
//   variant is just data with no zeros.
//
// hessian_block_size = n_latent (block-diagonal per patient).
// Compares Full HMC (mode 0) vs built-in Laplace (mode 1).
// ============================================================================

functions {
  real log_pred_burden(real dt, real init_log_dec, real init_log_gro,
                       real dec_rate, real gro_rate) {
    real state_dec = init_log_dec - dec_rate * dt;
    real state_gro = fmin(init_log_gro + gro_rate * dt, 500.0);
    return log_sum_exp(state_dec, state_gro);
  }

  real tumor_ll(vector theta,
                real tr_loc_pop, real tr_sd,
                real frac_logit_pop, real frac_sd,
                real init_logit_pop, real init_sd,
                real measure_sd, real log_lod,
                data int n_patients, data int n_latent,
                data vector normalized_obs,
                data array[] int patient_visit_pos,
                data array[] int visit_time_idx) {
    real lp = 0;
    for (i in 1:n_patients) {
      // Latents present depend on n_latent; absent ones are 0 (population).
      real z_tr   = theta[(i - 1) * n_latent + 1];
      real z_frac = n_latent >= 2 ? theta[(i - 1) * n_latent + 2] : 0.0;
      real z_init = n_latent >= 3 ? theta[(i - 1) * n_latent + 3] : 0.0;

      real tr_loc     = tr_loc_pop     + tr_sd   * z_tr;
      real frac_logit = frac_logit_pop + frac_sd * z_frac;
      real init_logit = init_logit_pop + init_sd * z_init;

      real log_dec_frac = log_inv_logit(frac_logit);
      real log_gro_frac = log1m_inv_logit(frac_logit);
      real dec_rate = exp(tr_loc + log_dec_frac);
      real gro_rate = exp(tr_loc + log_gro_frac);
      real init_log_dec = log_inv_logit(init_logit);
      real init_log_gro = log1m_inv_logit(init_logit);

      int v_start = patient_visit_pos[i];
      int v_end   = patient_visit_pos[i + 1] - 1;
      for (v in v_start:v_end) {
        real dt = visit_time_idx[v] - 1.0;
        real log_pred = log_pred_burden(dt, init_log_dec, init_log_gro,
                                        dec_rate, gro_rate);
        if (normalized_obs[v] > 0)
          lp += normal_lpdf(log(normalized_obs[v]) | log_pred, measure_sd);
        else
          lp += normal_lcdf(log_lod | log_pred, measure_sd);
      }
    }
    return lp;
  }

  matrix K_fn(int dim, int dummy) { return identity_matrix(dim); }
}

data {
  int<lower=1> n_patients;
  int<lower=1> n_total_visits;
  int<lower=1, upper=3> n_latent;
  vector[n_total_visits] normalized_obs;
  array[n_patients + 1] int patient_visit_pos;
  array[n_total_visits] int visit_time_idx;
  real<lower=0> measure_sd;
  real log_lod;
  int<lower=0, upper=1> laplace_mode;
}

transformed data {
  int latent_dim = n_patients * n_latent;
  vector[latent_dim] theta_0 = rep_vector(0.0, latent_dim);
  real tolerance = 1e-6;
  int max_num_steps = 500;
  int hessian_block_size = n_latent;
  int solver = 3;                 // LU of I+KW: tolerant of indefinite Hessian
  int max_steps_line_search = 100;
  int allow_fallback = 1;
}

parameters {
  real tr_loc_pop;
  real<lower=0> tr_sd;
  real frac_logit_pop;
  real<lower=0> frac_sd;
  real init_logit_pop;
  real<lower=0> init_sd;

  vector[laplace_mode == 0 ? n_patients : 0] z_tr;
  vector[(laplace_mode == 0 && n_latent >= 2) ? n_patients : 0] z_frac;
  vector[(laplace_mode == 0 && n_latent >= 3) ? n_patients : 0] z_init;
}

model {
  tr_loc_pop     ~ normal(-3, 1);
  tr_sd          ~ normal(0, 1);
  frac_logit_pop ~ normal(0.5, 1);
  frac_sd        ~ normal(0, 1);
  init_logit_pop ~ normal(0, 1);
  init_sd        ~ normal(0, 1);

  if (laplace_mode == 0) {
    z_tr ~ std_normal();
    if (n_latent >= 2) z_frac ~ std_normal();
    if (n_latent >= 3) z_init ~ std_normal();

    for (i in 1:n_patients) {
      real ztr   = z_tr[i];
      real zfrac = n_latent >= 2 ? z_frac[i] : 0.0;
      real zinit = n_latent >= 3 ? z_init[i] : 0.0;

      real tr_loc     = tr_loc_pop     + tr_sd   * ztr;
      real frac_logit = frac_logit_pop + frac_sd * zfrac;
      real init_logit = init_logit_pop + init_sd * zinit;

      real log_dec_frac = log_inv_logit(frac_logit);
      real log_gro_frac = log1m_inv_logit(frac_logit);
      real dec_rate = exp(tr_loc + log_dec_frac);
      real gro_rate = exp(tr_loc + log_gro_frac);
      real init_log_dec = log_inv_logit(init_logit);
      real init_log_gro = log1m_inv_logit(init_logit);

      int v_start = patient_visit_pos[i];
      int v_end   = patient_visit_pos[i + 1] - 1;
      for (v in v_start:v_end) {
        real dt = visit_time_idx[v] - 1.0;
        real log_pred = log_pred_burden(dt, init_log_dec, init_log_gro,
                                        dec_rate, gro_rate);
        if (normalized_obs[v] > 0)
          target += normal_lpdf(log(normalized_obs[v]) | log_pred, measure_sd);
        else
          target += normal_lcdf(log_lod | log_pred, measure_sd);
      }
    }
  } else {
    target += laplace_marginal_tol(
      tumor_ll,
      (tr_loc_pop, tr_sd, frac_logit_pop, frac_sd,
       init_logit_pop, init_sd, measure_sd, log_lod,
       n_patients, n_latent, normalized_obs, patient_visit_pos, visit_time_idx),
      hessian_block_size,
      K_fn,
      (latent_dim, latent_dim),
      (theta_0, tolerance, max_num_steps, solver, max_steps_line_search, allow_fallback)
    );
  }
}
