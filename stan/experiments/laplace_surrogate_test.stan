// Phase 1: log-concave QUADRATIC surrogate for backgrounded-trial patients,
// validated against the full bi-exponential via the built-in laplace_marginal_tol.
// ============================================================================
// Backgrounded patients' normalized log-burden is modeled as a quadratic in t:
//     log_burden_i(t) = beta_i0 + beta_i1*t + beta_i2*t^2 + N(0, measure_sd^2)
// Latents enter LINEARLY => Gaussian marginal => Laplace is EXACT => solver 1
// (PD-Hessian Cholesky) is valid (the bi-exponential needed solver 3 to survive).
//
// REVISED after the first gate FAIL (see spec section 3a):
//   (1) The intercept beta_0 is PINNED at beta_pop[1] (= g(0) = 0 by baseline
//       normalization), NOT marginalized. g(0)=0 for all theta makes its bridge
//       variance structurally zero => Sigma_beta rank-deficient => solver 1 fails
//       to factor a singular prior. So we integrate out only (beta_1, beta_2):
//       hessian_block_size = 2, latent_dim = n_patients * 2.
//   (2) Sigma_beta is the GAUSS-HERMITE pushforward Cov[beta(theta)] over
//       theta ~ N(theta_pop, diag(sd^2)) (3 nodes/dim, 27 evals), NOT the
//       first-order J diag(sd^2) J'. First-order captured only ~88% of the true
//       variance (biasing tr_sd low); GH-3 captures ~99.9%.
//
// (beta_pop, Sigma_beta) are NOT free: they are the mechanistic-bridge image of
// the population rate params (tr_loc_pop, frac_logit_pop, init_logit_pop) and the
// per-patient SDs (tr_sd, frac_sd, init_sd):
//   beta_pop   = Vinv * [g(t0), g(t1), g(t2)]   (g = exact bi-exponential log-burden)
//   Sigma_beta = GH-3 pushforward of (beta_1, beta_2)
//
// Modes (laplace_mode): 0 = full-HMC bi-exponential reference, 1 = surrogate-Laplace.
// ============================================================================

functions {
  // Exact bi-exponential predicted normalized log-burden at time offset t (weeks).
  // Same rate construction as the full model (sf.stanfunctions / tumor module).
  real bi_exp_log_burden(real t,
                         real tr_loc, real frac_logit, real init_logit) {
    real log_dec_frac = log_inv_logit(frac_logit);
    real log_gro_frac = log1m_inv_logit(frac_logit);
    real dec_rate = exp(tr_loc + log_dec_frac);
    real gro_rate = exp(tr_loc + log_gro_frac);
    real init_log_dec = log_inv_logit(init_logit);
    real init_log_gro = log1m_inv_logit(init_logit);
    real state_dec = init_log_dec - dec_rate * t;
    real state_gro = fmin(init_log_gro + gro_rate * t, 500.0);
    return log_sum_exp(state_dec, state_gro);
  }

  // Quadratic coefficients (beta0, beta1, beta2) interpolating the exact
  // bi-exponential log-burden at the 3 fixed anchor times. Vinv is the constant
  // inverse Vandermonde passed in as data. (beta0 == g(t0) == 0 when t0 == 0.)
  vector surrogate_anchor_betas(real tr_loc_pop, real frac_logit_pop,
                                real init_logit_pop,
                                data matrix Vinv, data vector anchor_times) {
    vector[3] g;
    for (k in 1:3)
      g[k] = bi_exp_log_burden(anchor_times[k], tr_loc_pop,
                               frac_logit_pop, init_logit_pop);
    return Vinv * g;
  }

  // Bridge covariance of the MARGINALIZED coefficients (beta_1, beta_2), via a
  // Gauss-Hermite pushforward of beta(theta) under theta ~ N(theta_pop,diag(sd^2)).
  // gh_x / gh_w are the constant standard-normal GH nodes/weights (weights sum 1).
  // Returns a 2x2 covariance; a diagonal jitter keeps it PD if an sd -> 0.
  matrix surrogate_bridge_cov(real tr_loc_pop, real frac_logit_pop,
                              real init_logit_pop,
                              real tr_sd, real frac_sd, real init_sd,
                              data matrix Vinv, data vector anchor_times,
                              data vector gh_x, data vector gh_w,
                              data real jitter) {
    int m = num_elements(gh_x);
    int nq = m * m * m;
    array[nq] vector[2] beta_pts;
    vector[nq] w;
    vector[2] mean_b = rep_vector(0.0, 2);
    int q = 1;
    for (i in 1:m) {
      for (j in 1:m) {
        for (k in 1:m) {
          real th1 = tr_loc_pop     + tr_sd   * gh_x[i];
          real th2 = frac_logit_pop + frac_sd * gh_x[j];
          real th3 = init_logit_pop + init_sd * gh_x[k];
          vector[3] beta_full =
            surrogate_anchor_betas(th1, th2, th3, Vinv, anchor_times);
          beta_pts[q] = beta_full[2:3];        // (beta_1, beta_2) only
          w[q] = gh_w[i] * gh_w[j] * gh_w[k];
          mean_b += w[q] * beta_pts[q];
          q += 1;
        }
      }
    }
    matrix[2, 2] S = rep_matrix(0.0, 2, 2);
    for (r in 1:nq)
      S += w[r] * (beta_pts[r] - mean_b) * (beta_pts[r] - mean_b)';
    return S + diag_matrix(rep_vector(jitter, 2));
  }

  // laplace_marginal_tol functor. theta = stacked per-patient [db1, db2] latents
  // (deviations of the slope/curvature from beta_pop). The intercept is pinned at
  // beta_pop[1] (= 0); only beta_1, beta_2 are marginalized. K_fn returns the
  // block-diagonal Sigma_beta, so theta is in the natural (b1,b2) coordinate
  // centered at (beta_pop[2], beta_pop[3]).
  real surrogate_ll(vector theta,
                    vector beta_pop,
                    real measure_sd, real log_lod,
                    data int n_patients,
                    data vector normalized_obs,
                    data array[] int patient_visit_pos,
                    data array[] int visit_time,
                    data array[] int patient_of_visit) {
    int d = 2;                       // marginalized latents per patient: b1, b2
    real lp = 0;
    for (i in 1:n_patients) {
      real b0 = beta_pop[1];                                  // pinned (= 0)
      real b1 = beta_pop[2] + theta[(i - 1) * d + 1];
      real b2 = beta_pop[3] + theta[(i - 1) * d + 2];
      int v_start = patient_visit_pos[i];
      int v_end   = patient_visit_pos[i + 1] - 1;
      for (v in v_start:v_end) {
        real t = visit_time[v] - 1.0;
        real mu = b0 + b1 * t + b2 * t * t;
        if (normalized_obs[v] > 0)
          lp += normal_lpdf(log(normalized_obs[v]) | mu, measure_sd);
        else
          lp += normal_lcdf(log_lod | mu, measure_sd);
      }
    }
    return lp;
  }

  // Prior covariance functor: the bridged 2x2 Sigma_beta tiled block-diagonally,
  // one 2x2 block per patient (matches hessian_block_size = 2).
  matrix surrogate_K_fn(matrix Sigma_beta, int n_patients) {
    int d = 2;
    matrix[n_patients * d, n_patients * d] K =
      rep_matrix(0, n_patients * d, n_patients * d);
    for (i in 1:n_patients) {
      int s = (i - 1) * d + 1;
      K[s:(s + d - 1), s:(s + d - 1)] = Sigma_beta;
    }
    return K;
  }
}

data {
  int<lower=1> n_patients;
  int<lower=1> n_total_visits;
  vector[n_total_visits] normalized_obs;       // <= 0 marks below-LOD
  array[n_patients + 1] int patient_visit_pos;
  array[n_total_visits] int visit_time;        // 1-indexed week (t = idx - 1)
  array[n_total_visits] int patient_of_visit;  // patient id per visit row
  real<lower=0> measure_sd;
  real log_lod;
  vector[3] anchor_times;                      // fixed calendar anchors (weeks)
  int<lower=0, upper=1> laplace_mode;
}

transformed data {
  int latent_dim = n_patients * 2;             // marginalize (b1, b2) per patient
  vector[latent_dim] theta_0 = rep_vector(0.0, latent_dim);
  real tolerance = 1e-8;
  int max_num_steps = 100;        // log-concave => Newton converges fast
  int hessian_block_size = 2;     // 2 marginalized latents per patient
  int solver = 1;                 // PD-Hessian Cholesky: valid for log-concave
  int max_steps_line_search = 0;  // not needed when well-conditioned
  int allow_fallback = 1;
  real jitter = 1e-10;

  // Constant inverse Vandermonde for the 3 fixed anchors.
  matrix[3, 3] V;
  for (k in 1:3) {
    V[k, 1] = 1.0;
    V[k, 2] = anchor_times[k];
    V[k, 3] = anchor_times[k] * anchor_times[k];
  }
  matrix[3, 3] Vinv = inverse(V);

  // 3-point Gauss-Hermite nodes/weights for the standard normal (weights sum 1).
  vector[3] gh_x = [-sqrt(3.0), 0.0, sqrt(3.0)]';
  vector[3] gh_w = [1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0]';
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
  tr_loc_pop     ~ normal(-3, 1);
  tr_sd          ~ normal(0, 1);
  frac_logit_pop ~ normal(0.5, 1);
  frac_sd        ~ normal(0, 1);
  init_logit_pop ~ normal(0, 1);
  init_sd        ~ normal(0, 1);

  if (laplace_mode == 0) {
    // MODE 0: full-HMC bi-exponential reference (the true generative model)
    z_tr ~ std_normal();
    z_frac ~ std_normal();
    z_init ~ std_normal();
    for (i in 1:n_patients) {
      real tr_loc     = tr_loc_pop     + tr_sd   * z_tr[i];
      real frac_logit = frac_logit_pop + frac_sd * z_frac[i];
      real init_logit = init_logit_pop + init_sd * z_init[i];
      int v_start = patient_visit_pos[i];
      int v_end   = patient_visit_pos[i + 1] - 1;
      for (v in v_start:v_end) {
        real t = visit_time[v] - 1.0;
        real mu = bi_exp_log_burden(t, tr_loc, frac_logit, init_logit);
        if (normalized_obs[v] > 0)
          target += normal_lpdf(log(normalized_obs[v]) | mu, measure_sd);
        else
          target += normal_lcdf(log_lod | mu, measure_sd);
      }
    }
  } else {
    // MODE 1: surrogate-Laplace (the thing under test)
    vector[3] beta_pop = surrogate_anchor_betas(
      tr_loc_pop, frac_logit_pop, init_logit_pop, Vinv, anchor_times);
    matrix[2, 2] Sigma_beta = surrogate_bridge_cov(
      tr_loc_pop, frac_logit_pop, init_logit_pop,
      tr_sd, frac_sd, init_sd, Vinv, anchor_times, gh_x, gh_w, jitter);

    target += laplace_marginal_tol(
      surrogate_ll,
      (beta_pop, measure_sd, log_lod, n_patients,
       normalized_obs, patient_visit_pos, visit_time, patient_of_visit),
      hessian_block_size,
      surrogate_K_fn,
      (Sigma_beta, n_patients),
      (theta_0, tolerance, max_num_steps, solver, max_steps_line_search, allow_fallback)
    );
  }
}
