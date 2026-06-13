// Phase 1: 3D bi-exponential TUMOR (SLD) Laplace experiment
//          using the built-in laplace_marginal_tol (Stan >= 2.39)
// ============================================================================
// GO/NO-GO gate for marginalizing patient-level latents of a *backgrounded*
// historical trial in sf-ssls-lfo.stan.
//
// This mirrors stan/experiments/laplace_psa_test.stan but uses the REAL tumor
// observation model from stan/modules/state_space/sf.stanfunctions:
//   sf_log_space_obs_lpdf — Gaussian on the log scale for observed SLD, and
//   left-censoring via normal_lcdf(log_lod | log_pred, sd) for below-LOD visits.
// (The docstring there mentions Student-t, but the code uses normal_lpdf /
//  normal_lcdf — we replicate the *code*, not the docstring.)
//
// Latents per patient: z_tr, z_frac, z_init (the bi-exponential NCP triple).
//   tr_loc      = tr_loc_pop      + tr_sd   * z_tr
//   frac_logit  = frac_logit_pop  + frac_sd * z_frac
//   init_logit  = init_logit_pop  + init_sd * z_init
//
// Two fitting modes controlled by `laplace_mode`:
//   0 = Full HMC (z params are sampled, standard NCP)
//   1 = Built-in Laplace marginalization (z params integrated out)
//
// For a CORRECT Laplace result the population posteriors from mode 0 and mode 1
// should agree within a few MCSE. Because the log-SLD-vs-z map runs through a
// nonlinear log_sum_exp (and a normal_lcdf censoring tail), this likelihood is
// NOT guaranteed log-concave in z — which is exactly why we gate on agreement.
// ============================================================================

functions {
  // Bi-exponential predicted log-burden at a single time offset dt (weeks since
  // baseline). Shared by the likelihood functor and the HMC path so the two
  // modes are byte-for-byte the same model.
  real log_pred_burden(real dt,
                       real init_log_dec, real init_log_gro,
                       real dec_rate, real gro_rate) {
    real state_dec = init_log_dec - dec_rate * dt;
    real state_gro = fmin(init_log_gro + gro_rate * dt, 500.0);
    return log_sum_exp(state_dec, state_gro);
  }

  // Log-likelihood functor for laplace_marginal_tol.
  // theta is the stacked NCP vector: [z_tr_1, z_frac_1, z_init_1, z_tr_2, ...].
  // All remaining (non-theta) args are passed as a tuple and marked `data`
  // where they are pure data (required for the higher-order autodiff path).
  real tumor_ll(vector theta,
                real tr_loc_pop, real tr_sd,
                real frac_logit_pop, real frac_sd,
                real init_logit_pop, real init_sd,
                real measure_sd, real log_lod,
                data int n_patients,
                data vector normalized_obs,
                data array[] int patient_visit_pos,
                data array[] int visit_time_idx) {

    int d = 3;  // NCP params per patient
    real lp = 0;

    for (i in 1:n_patients) {
      real z_tr   = theta[(i - 1) * d + 1];
      real z_frac = theta[(i - 1) * d + 2];
      real z_init = theta[(i - 1) * d + 3];

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
        if (normalized_obs[v] > 0) {
          lp += normal_lpdf(log(normalized_obs[v]) | log_pred, measure_sd);
        } else {
          // Below LOD: left-censored, matches sf_log_space_obs_lpdf.
          lp += normal_lcdf(log_lod | log_pred, measure_sd);
        }
      }
    }
    return lp;
  }

  // Prior covariance functor: identity matrix (standard normal NCP prior).
  // Args must be a tuple of >= 2 ints; (dim, dim) is a convenient dummy pair.
  matrix K_fn(int dim, int dummy) {
    return identity_matrix(dim);
  }
}

data {
  int<lower=1> n_patients;
  int<lower=1> n_total_visits;
  vector[n_total_visits] normalized_obs;        // <= 0 marks below-LOD (censored)
  array[n_patients + 1] int patient_visit_pos;
  array[n_total_visits] int visit_time_idx;     // 1-indexed week (dt = idx - 1)
  real<lower=0> measure_sd;
  real log_lod;                                 // log normalized limit of detection
  int<lower=0, upper=1> laplace_mode;
}

transformed data {
  int latent_dim = n_patients * 3;
  vector[latent_dim] theta_0 = rep_vector(0.0, latent_dim);
  real tolerance = 1e-6;
  int max_num_steps = 500;       // raised from 100: some hyperparam draws need more
  int hessian_block_size = 3;    // 3 NCP params per patient -> block diagonal
  // solver 3 = LU of (I + KW): the "last resort" that does NOT require the
  // log-likelihood Hessian W to be positive definite. The tumor likelihood is
  // non-log-concave in z (log_sum_exp + normal_lcdf tail), so solver 1
  // ("Hessian-root Cholesky") fails with "Hessian not positive definite".
  int solver = 3;
  int max_steps_line_search = 100;  // generous Wolfe line search for stability
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
  vector[laplace_mode == 0 ? n_patients : 0] z_frac;
  vector[laplace_mode == 0 ? n_patients : 0] z_init;
}

model {
  // Priors (match laplace_psa_test.stan so the two experiments are comparable)
  tr_loc_pop     ~ normal(-3, 1);
  tr_sd          ~ normal(0, 1);
  frac_logit_pop ~ normal(0.5, 1);
  frac_sd        ~ normal(0, 1);
  init_logit_pop ~ normal(0, 1);
  init_sd        ~ normal(0, 1);

  if (laplace_mode == 0) {
    // MODE 0: Full HMC over the NCP latents
    z_tr   ~ std_normal();
    z_frac ~ std_normal();
    z_init ~ std_normal();

    for (i in 1:n_patients) {
      real tr_loc     = tr_loc_pop     + tr_sd   * z_tr[i];
      real frac_logit = frac_logit_pop + frac_sd * z_frac[i];
      real init_logit = init_logit_pop + init_sd * z_init[i];

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
        if (normalized_obs[v] > 0) {
          target += normal_lpdf(log(normalized_obs[v]) | log_pred, measure_sd);
        } else {
          target += normal_lcdf(log_lod | log_pred, measure_sd);
        }
      }
    }
  } else {
    // MODE 1: Built-in Laplace marginalization of the NCP latents
    target += laplace_marginal_tol(
      tumor_ll,
      (tr_loc_pop, tr_sd, frac_logit_pop, frac_sd,
       init_logit_pop, init_sd, measure_sd, log_lod,
       n_patients, normalized_obs, patient_visit_pos, visit_time_idx),
      hessian_block_size,
      K_fn,
      (latent_dim, latent_dim),
      (theta_0, tolerance, max_num_steps, solver, max_steps_line_search, allow_fallback)
    );
  }
}
