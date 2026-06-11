// Phase 1d: JOINT SLD + multistate(0->1, 0->3) surrogate gate WITH CORRELATED
// PATIENT-LEVEL FRAILTY. Contract §9 step 2.
// ============================================================================
// Extends laplace_joint_test.stan (the validated d=2 burden-only gate) to the
// publication design: marginalize each background patient's latent vector
//   theta_i = (b1_i, b2_i, u01_i, u03_i)        d = 2 + n_frailty_slots = 4
// where
//   (b1,b2)   = surrogate quadratic burden slope/curvature (intercept pinned 0),
//   (u01,u03) = CORRELATED patient-level frailty intercepts on the 0->1 and 0->3
//               baseline log-hazards (the model's ms_corr_u patient-level block).
//
// Both hazards couple to burden via the (level, velocity) basis (the step-0
// production basis), feeding the SURROGATE quadratic:
//   level    g(w)  = b1*w + b2*w^2        (linear in (b1,b2))
//   velocity g'(w) = b1 + 2*b2*w          (linear in (b1,b2))
//   loghaz_T(w) = base_T + u_T
//               + cf_lvl_T * std_lvl(g(w)) + cf_vel_T * std_vel(g'(w))
// 0->1 is visit-gated (hazard only at observed visit weeks); 0->3 is continuous
// (every week). Frailty enters as a pure ADDITIVE intercept shift -> the survival
// term stays -sum exp(linear) => log-concave in all 4 latents (the property this
// gate verifies numerically before production).
//
// Prior covariance per patient: K_i = blockdiag(Sigma_beta[2x2], Sigma_u[2x2]),
//   Sigma_u = diag(s01,s03) * (L L') * diag(s01,s03)   (matches production
//   transformed_parameters.stan:37 ms_corr_u construction). Zero cross-block
//   (burden RE independent of frailty RE).
//
// Modes (laplace_mode): 0 = full HMC over (b1,b2,u01,u03); 1 = joint-Laplace.
// PASS: pop params AND coef_lvl/vel_01/03 AND sigma_01/03 AND rho(01,03) agree
// HMC-vs-Laplace; joint-Laplace converges clean on solver 1.
// ============================================================================

functions {
  real std_lvl(real lvl, data real median_lvl, data real iqr_lvl) {
    return (lvl - median_lvl) / iqr_lvl;
  }
  real std_vel(real vel, data real median_vel, data real iqr_vel) {
    return (vel - median_vel) / iqr_vel;
  }

  // Joint per-patient log-likelihood functor for laplace_marginal_tol.
  // WHITENED (NCP) latents: theta carries UNIT-scale latents; physical
  // parameters are recovered inside the functor via the population means and
  // scales. This keeps the prior covariance K near-identity (well-conditioned),
  // moving the burden vs frailty scale mismatch out of the Hessian the inner
  // Newton solver inverts. theta block per patient (d=4):
  //   [zb1, zb2, zu01, zu03]  with
  //     b1  = b1_pop + b1_sd * zb1,  b2 = b2_pop + b2_sd * zb2,
  //     (u01,u03) = diag(s01,s03) * L_frailty * (zu01, zu03)
  // Prior on theta is N(0, I) (burden) x N(0, LL') (frailty correlation only).
  real joint_frailty_ll(vector theta,
                        vector beta_pop,              // (b1_pop, b2_pop)
                        real measure_sd,
                        real base_01, real cf_lvl_01, real cf_vel_01,
                        real base_03, real cf_lvl_03, real cf_vel_03,
                        data real median_lvl, data real iqr_lvl,
                        data real median_vel, data real iqr_vel,
                        data int n_patients,
                        data int n_wk,
                        data vector sld_obs,
                        data array[] int sld_pos,
                        data array[] int sld_time,
                        data array[] int sld_is_visit,   // 1 if this SLD obs week is a 0->1 gating visit
                        data array[] int ms_event_wk_01, data array[] int ms_censored_01,
                        data array[] int ms_event_wk_03, data array[] int ms_censored_03,
                        data array[,] int visit_wk,      // [n_patients, n_wk] 1 if week is an observed visit (0->1 gating)
                        data int d) {
    real lp = 0;
    // Latents enter RAW (no internal scaling); the prior covariance K carries the
    // burden + frailty scales/correlation (see joint_frailty_K_fn). This matches
    // the validated d=2 toy and production surrogate_ll: laplace_marginal uses K
    // to regularize the mode-finding, so the scale MUST live in K, not the functor
    // (whitening to K=I left the burden directions near-flat => singular Hessian
    // => "not finite at initial theta"; the d=4 debugging root cause).
    for (i in 1:n_patients) {
      real b1  = beta_pop[1] + theta[(i - 1) * d + 1];
      real b2  = beta_pop[2] + theta[(i - 1) * d + 2];
      real u01 = theta[(i - 1) * d + 3];
      real u03 = theta[(i - 1) * d + 4];

      // --- SLD term ---
      int s_start = sld_pos[i];
      int s_end   = sld_pos[i + 1] - 1;
      for (v in s_start:s_end) {
        real t = sld_time[v];
        real mu = b1 * t + b2 * t * t;
        lp += normal_lpdf(sld_obs[v] | mu, measure_sd);
      }

      // --- 0->1 visit-gated hazard (only at observed visit weeks) ---
      int te1 = ms_censored_01[i] ? n_wk : ms_event_wk_01[i];
      real cum_haz_01 = 0;
      for (w in 1:te1) {
        if (visit_wk[i, w] == 1) {
          real lvl = b1 * w + b2 * w * w;
          real vel = b1 + 2 * b2 * w;
          real loghaz = base_01 + u01
                      + cf_lvl_01 * std_lvl(lvl, median_lvl, iqr_lvl)
                      + cf_vel_01 * std_vel(vel, median_vel, iqr_vel);
          cum_haz_01 += exp(loghaz);
          if (!ms_censored_01[i] && w == te1) lp += loghaz;
        }
      }
      lp += -cum_haz_01;

      // --- 0->3 continuous hazard (every week) ---
      int te3 = ms_censored_03[i] ? n_wk : ms_event_wk_03[i];
      real cum_haz_03 = 0;
      for (w in 1:te3) {
        real lvl = b1 * w + b2 * w * w;
        real vel = b1 + 2 * b2 * w;
        real loghaz = base_03 + u03
                    + cf_lvl_03 * std_lvl(lvl, median_lvl, iqr_lvl)
                    + cf_vel_03 * std_vel(vel, median_vel, iqr_vel);
        cum_haz_03 += exp(loghaz);
        if (!ms_censored_03[i] && w == te3) lp += loghaz;
      }
      lp += -cum_haz_03;
    }
    return lp;
  }

  // Per-patient prior covariance block, tiled block-diagonal (matches the
  // validated d=2 toy / production surrogate_K_fn pattern: the REAL scale lives
  // here, NOT in the functor — laplace_marginal uses K to regularize mode-finding).
  // K_i = blockdiag(Sigma_beta[2x2], Sigma_u[2x2]); zero burden<->frailty cross-block.
  matrix joint_frailty_K_fn(matrix Sigma_beta, matrix Sigma_u,
                            int n_patients, data int d) {
    matrix[n_patients * d, n_patients * d] K = rep_matrix(0, n_patients * d, n_patients * d);
    for (i in 1:n_patients) {
      int s = (i - 1) * d + 1;
      K[s:(s + 1), s:(s + 1)]             = Sigma_beta;   // burden block (rows 1-2)
      K[(s + 2):(s + 3), (s + 2):(s + 3)] = Sigma_u;      // frailty block (rows 3-4)
    }
    return K;
  }
}

data {
  int<lower=1> n_patients;
  int<lower=1> n_total_sld;
  int<lower=1> n_wk;
  vector[n_total_sld] sld_obs;
  array[n_patients + 1] int sld_pos;
  array[n_total_sld] int sld_time;
  array[n_total_sld] int sld_is_visit;
  array[n_patients] int ms_event_wk_01, ms_censored_01;
  array[n_patients] int ms_event_wk_03, ms_censored_03;
  array[n_patients, n_wk] int visit_wk;     // 0->1 visit-gating mask
  real median_lvl, iqr_lvl, median_vel, iqr_vel;
  real<lower=0> measure_sd;
  int<lower=0, upper=1> laplace_mode;
}

transformed data {
  int d = 4;                         // 2 burden + 2 frailty
  int latent_dim = n_patients * d;
  vector[latent_dim] theta_0 = rep_vector(0.0, latent_dim);
  real tolerance = 1e-8;
  int max_num_steps = 100;
  int hessian_block_size = d;        // full per-patient block (burden+frailty correlated within block)
  int solver = 1;                    // log-concave joint => solver 1
  int max_steps_line_search = 0;
  int allow_fallback = 1;
}

parameters {
  // Burden population
  real b1_pop;
  real b2_pop;
  real<lower=0> b1_sd;
  real<lower=0> b2_sd;
  // 0->1 hazard
  real base_01;
  real cf_lvl_01;
  real cf_vel_01;
  // 0->3 hazard
  real base_03;
  real cf_lvl_03;
  real cf_vel_03;
  // Correlated frailty scale + correlation (patient-level block)
  real<lower=0> s01;
  real<lower=0> s03;
  cholesky_factor_corr[2] L_frailty;
  // HMC latents (mode 0 only)
  vector[laplace_mode == 0 ? n_patients : 0] z_b1;
  vector[laplace_mode == 0 ? n_patients : 0] z_b2;
  matrix[2, laplace_mode == 0 ? n_patients : 0] z_frailty;   // std_normal, scaled by diag(s)*L
}

model {
  b1_pop  ~ normal(-0.04, 0.05);
  b2_pop  ~ normal(0.001, 0.005);
  b1_sd   ~ normal(0, 0.05);
  b2_sd   ~ normal(0, 0.005);
  base_01 ~ normal(-4, 1);
  base_03 ~ normal(-4, 1);
  cf_lvl_01 ~ normal(0, 1);
  cf_vel_01 ~ normal(0, 1);
  cf_lvl_03 ~ normal(0, 1);
  cf_vel_03 ~ normal(0, 1);
  s01 ~ normal(0, 1);
  s03 ~ normal(0, 1);
  L_frailty ~ lkj_corr_cholesky(2);

  if (laplace_mode == 0) {
    z_b1 ~ std_normal();
    z_b2 ~ std_normal();
    to_vector(z_frailty) ~ std_normal();
    // u_i = diag(s01,s03) * L_frailty * z_i  (matches ms_corr_u construction)
    matrix[2, n_patients] u = diag_pre_multiply([s01, s03]', L_frailty) * z_frailty;
    for (i in 1:n_patients) {
      real b1 = b1_pop + b1_sd * z_b1[i];
      real b2 = b2_pop + b2_sd * z_b2[i];
      real u01 = u[1, i];
      real u03 = u[2, i];

      int s_start = sld_pos[i];
      int s_end   = sld_pos[i + 1] - 1;
      for (v in s_start:s_end) {
        real t = sld_time[v];
        target += normal_lpdf(sld_obs[v] | b1 * t + b2 * t * t, measure_sd);
      }

      int te1 = ms_censored_01[i] ? n_wk : ms_event_wk_01[i];
      real cum_haz_01 = 0;
      for (w in 1:te1) {
        if (visit_wk[i, w] == 1) {
          real lvl = b1 * w + b2 * w * w;
          real vel = b1 + 2 * b2 * w;
          real loghaz = base_01 + u01
                      + cf_lvl_01 * std_lvl(lvl, median_lvl, iqr_lvl)
                      + cf_vel_01 * std_vel(vel, median_vel, iqr_vel);
          cum_haz_01 += exp(loghaz);
          if (ms_censored_01[i] == 0 && w == te1) target += loghaz;
        }
      }
      target += -cum_haz_01;

      int te3 = ms_censored_03[i] ? n_wk : ms_event_wk_03[i];
      real cum_haz_03 = 0;
      for (w in 1:te3) {
        real lvl = b1 * w + b2 * w * w;
        real vel = b1 + 2 * b2 * w;
        real loghaz = base_03 + u03
                    + cf_lvl_03 * std_lvl(lvl, median_lvl, iqr_lvl)
                    + cf_vel_03 * std_vel(vel, median_vel, iqr_vel);
        cum_haz_03 += exp(loghaz);
        if (ms_censored_03[i] == 0 && w == te3) target += loghaz;
      }
      target += -cum_haz_03;
    }
  } else {
    vector[2] beta_pop = [b1_pop, b2_pop]';
    matrix[2, 2] Sigma_beta = diag_matrix(square([b1_sd, b2_sd]') + 1e-10);
    // Frailty covariance: diag(s)*LL'*diag(s) (matches ms_corr_u construction).
    matrix[2, 2] Lu = diag_pre_multiply([s01, s03]', L_frailty);
    matrix[2, 2] Sigma_u = Lu * Lu' + diag_matrix(rep_vector(1e-10, 2));
    target += laplace_marginal_tol(
      joint_frailty_ll,
      (beta_pop, measure_sd, base_01, cf_lvl_01, cf_vel_01,
       base_03, cf_lvl_03, cf_vel_03,
       median_lvl, iqr_lvl, median_vel, iqr_vel,
       n_patients, n_wk, sld_obs, sld_pos, sld_time, sld_is_visit,
       ms_event_wk_01, ms_censored_01, ms_event_wk_03, ms_censored_03,
       visit_wk, d),
      hessian_block_size,
      joint_frailty_K_fn,
      (Sigma_beta, Sigma_u, n_patients, d),
      (theta_0, tolerance, max_num_steps, solver, max_steps_line_search, allow_fallback)
    );
  }
}

generated quantities {
  // recovered 01<->03 frailty correlation
  real rho_frailty = multiply_lower_tri_self_transpose(L_frailty)[1, 2];
}
