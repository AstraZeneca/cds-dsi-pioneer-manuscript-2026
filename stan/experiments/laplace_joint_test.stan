// Phase 1c: JOINT SLD + multistate(0->1) surrogate gate.
// ============================================================================
// Validates marginalizing background patients' tumor latents from the JOINT
// likelihood (SLD + piecewise-exponential 0->1 hazard coupled to the tumor
// burden), not SLD alone. Confirmed log-concave in the quadratic coefficients
// (R check: max Hessian eig -6188 << 0), because feeding the SURROGATE quadratic
// burden (linear in the coefficients) into the hazard keeps the survival term
// log-concave (event term linear; -integral-of-exp(linear) concave).
//
// Latents per patient: (b1, b2) = surrogate quadratic slope/curvature (intercept
// pinned at 0 by baseline normalization). They drive BOTH:
//   - SLD: log_burden(t) = b1*t + b2*t^2 + N(0, measure_sd)
//   - 0->1 hazard: loghaz(w) = base_01 + coef_01 * standardize(b1*w + b2*w^2)
//     piecewise-exponential over weekly grid; survival = -sum_w exp(loghaz(w)).
//
// Modes (laplace_mode): 0 = full HMC over (b1,b2); 1 = joint-Laplace marginalization.
// PASS: population params AND coef_01 (the burden->hazard coupling) agree HMC vs Laplace.
// ============================================================================

functions {
  // standardize log-burden: affine => linear in (b1,b2). data median/iqr.
  real standardize_lb(real lb, data real median_lb, data real iqr_lb) {
    return (lb - median_lb) / iqr_lb;
  }

  // Joint per-patient log-likelihood functor for laplace_marginal_tol.
  // theta = stacked [b1_1,b2_1, b1_2,b2_2, ...] (2 latents/patient, intercept pinned 0).
  real joint_ll(vector theta,
                vector beta_pop,              // (b1_pop, b2_pop) bridge mean (here length 2)
                real measure_sd,
                real base_01, real coef_01,
                data real median_lb, data real iqr_lb,
                data int n_patients,
                data int n_wk,                // weekly hazard grid length
                data vector sld_obs,          // log-SLD observations, stacked
                data array[] int sld_pos,     // per-patient position into sld_obs
                data array[] int sld_time,    // visit week (t) per sld obs
                data array[] int ms_event_wk, // 0->1 event week per patient (0 if censored at n_wk)
                data array[] int ms_censored) {
    int d = 2;
    real lp = 0;
    for (i in 1:n_patients) {
      real b1 = beta_pop[1] + theta[(i - 1) * d + 1];
      real b2 = beta_pop[2] + theta[(i - 1) * d + 2];

      // --- SLD term ---
      int s_start = sld_pos[i];
      int s_end   = sld_pos[i + 1] - 1;
      for (v in s_start:s_end) {
        real t = sld_time[v];
        real mu = b1 * t + b2 * t * t;
        lp += normal_lpdf(sld_obs[v] | mu, measure_sd);
      }

      // --- 0->1 piecewise-exponential hazard term (burden-coupled) ---
      int te = ms_censored[i] ? n_wk : ms_event_wk[i];
      real cum_haz = 0;
      for (w in 1:te) {
        real muw = b1 * w + b2 * w * w;
        real loghaz = base_01 + coef_01 * standardize_lb(muw, median_lb, iqr_lb);
        cum_haz += exp(loghaz);
        if (!ms_censored[i] && w == te) lp += loghaz;  // event contribution
      }
      lp += -cum_haz;  // integrated hazard (survival)
    }
    return lp;
  }

  // Identity-ish prior covariance: 2x2 per patient block. Here we pass the bridge
  // covariance directly (diagonal for the toy); tiled block-diagonal.
  matrix joint_K_fn(matrix Sigma_beta, int n_patients) {
    int d = 2;
    matrix[n_patients * d, n_patients * d] K = rep_matrix(0, n_patients * d, n_patients * d);
    for (i in 1:n_patients) {
      int s = (i - 1) * d + 1;
      K[s:(s + d - 1), s:(s + d - 1)] = Sigma_beta;
    }
    return K;
  }
}

data {
  int<lower=1> n_patients;
  int<lower=1> n_total_sld;
  int<lower=1> n_wk;                       // weekly hazard grid (e.g. 28)
  vector[n_total_sld] sld_obs;             // log-SLD
  array[n_patients + 1] int sld_pos;
  array[n_total_sld] int sld_time;
  array[n_patients] int ms_event_wk;       // 0->1 event week (1..n_wk), or n_wk if censored
  array[n_patients] int ms_censored;       // 1 if censored
  real median_lb;
  real iqr_lb;
  real<lower=0> measure_sd;
  int<lower=0, upper=1> laplace_mode;
}

transformed data {
  int latent_dim = n_patients * 2;
  vector[latent_dim] theta_0 = rep_vector(0.0, latent_dim);
  real tolerance = 1e-8;
  int max_num_steps = 100;
  int hessian_block_size = 2;
  int solver = 1;                 // log-concave joint => PD Hessian => solver 1
  int max_steps_line_search = 0;
  int allow_fallback = 1;
}

parameters {
  // Population burden coefficients (the bridge mean, sampled here for the toy)
  real b1_pop;
  real b2_pop;
  real<lower=0> b1_sd;
  real<lower=0> b2_sd;
  // MS 0->1 hazard population params (THE coupling is coef_01)
  real base_01;
  real coef_01;
  // HMC latents (mode 0 only)
  vector[laplace_mode == 0 ? n_patients : 0] z_b1;
  vector[laplace_mode == 0 ? n_patients : 0] z_b2;
}

model {
  b1_pop  ~ normal(-0.04, 0.05);
  b2_pop  ~ normal(0.001, 0.005);
  b1_sd   ~ normal(0, 0.05);
  b2_sd   ~ normal(0, 0.005);
  base_01 ~ normal(-4, 1);
  coef_01 ~ normal(0, 1);

  if (laplace_mode == 0) {
    z_b1 ~ std_normal();
    z_b2 ~ std_normal();
    for (i in 1:n_patients) {
      real b1 = b1_pop + b1_sd * z_b1[i];
      real b2 = b2_pop + b2_sd * z_b2[i];
      int s_start = sld_pos[i];
      int s_end   = sld_pos[i + 1] - 1;
      for (v in s_start:s_end) {
        real t = sld_time[v];
        real mu = b1 * t + b2 * t * t;
        target += normal_lpdf(sld_obs[v] | mu, measure_sd);
      }
      int te = ms_censored[i] ? n_wk : ms_event_wk[i];
      real cum_haz = 0;
      for (w in 1:te) {
        real muw = b1 * w + b2 * w * w;
        real loghaz = base_01 + coef_01 * standardize_lb(muw, median_lb, iqr_lb);
        cum_haz += exp(loghaz);
        if (ms_censored[i] == 0 && w == te) target += loghaz;
      }
      target += -cum_haz;
    }
  } else {
    vector[2] beta_pop = [b1_pop, b2_pop]';
    matrix[2, 2] Sigma_beta = diag_matrix(square([b1_sd, b2_sd]') + 1e-10);
    target += laplace_marginal_tol(
      joint_ll,
      (beta_pop, measure_sd, base_01, coef_01, median_lb, iqr_lb,
       n_patients, n_wk, sld_obs, sld_pos, sld_time, ms_event_wk, ms_censored),
      hessian_block_size,
      joint_K_fn,
      (Sigma_beta, n_patients),
      (theta_0, tolerance, max_num_steps, solver, max_steps_line_search, allow_fallback)
    );
  }
}
