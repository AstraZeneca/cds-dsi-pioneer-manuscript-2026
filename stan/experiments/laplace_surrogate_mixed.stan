// Phase 1b: MIXED-COHORT gate for the backgrounded-trial surrogate.
// ============================================================================
// This is the test that matches PRODUCTION (spec section 5, revised): forecast
// patients use the full bi-exponential (they pin tr_sd/frac_sd/init_sd), while
// backgrounded patients use the log-concave quadratic surrogate (they contribute
// population MEANS + whatever SD info the quadratic carries). The all-backgrounded
// toy (laplace_surrogate_test.stan) was stricter than production: it forced the
// surrogate to identify all 3 population SDs alone, which 2 coefficients cannot do
// (tr_sd's footprint is the ~1250x-weaker curvature channel, below the noise floor).
//
// Patients are ordered FORECAST-FIRST: indices 1..n_forecast are forecast,
// (n_forecast+1)..n_patients are backgrounded.
//
// Modes (mixed_mode):
//   0 = reference: ALL patients full-HMC bi-exponential.
//   1 = mixed:     forecast = full-HMC bi-exponential; backgrounded = surrogate.
// Both contribute to the SAME population parameters. PASS = the two runs' population
// posteriors agree (the surrogate faithfully stands in for the backgrounded subset).
//
// Surrogate internals (validated separately): quadratic in t, intercept pinned at
// g(0)=0, marginalize (beta_1, beta_2), GH-3 pushforward covariance bridge,
// hessian_block_size = 2.
// ============================================================================

functions {
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

  // Quadratic coefficients interpolating g at the anchors, for a SINGLE theta.
  vector surrogate_anchor_betas(real tr_loc_pop, real frac_logit_pop,
                                real init_logit_pop,
                                data matrix Vinv, data vector anchor_times) {
    vector[3] g;
    for (k in 1:3)
      g[k] = bi_exp_log_burden(anchor_times[k], tr_loc_pop,
                               frac_logit_pop, init_logit_pop);
    return Vinv * g;
  }

  // GH-3 pushforward of beta(theta), theta ~ N(theta_pop, diag(sd^2)). Returns
  // BOTH the population mean E[beta] AND the marginalized-coef covariance
  // Cov[(beta_1,beta_2)] from the SAME quadrature:
  //   - mean E[beta] (NOT g(theta_pop)): the population mean log-burden is the
  //     average of trajectories, not the trajectory of the average patient. g is
  //     nonlinear so these differ by a Jensen gap that biases the location params
  //     when many backgrounded patients make them precise (the mixed-gate FAIL).
  //   - Cov[(beta_1,beta_2)]: the random-effect spread (intercept pinned at 0).
  // Tuple return: (mean_beta[3], cov_beta12[2,2]).
  tuple(vector, matrix) surrogate_bridge(real tr_loc_pop, real frac_logit_pop,
                                         real init_logit_pop,
                                         real tr_sd, real frac_sd, real init_sd,
                                         data matrix Vinv, data vector anchor_times,
                                         data vector gh_x, data vector gh_w,
                                         data real jitter) {
    int m = num_elements(gh_x);
    int nq = m * m * m;
    array[nq] vector[3] beta_pts;
    vector[nq] w;
    vector[3] mean_full = rep_vector(0.0, 3);
    int q = 1;
    for (i in 1:m) {
      for (j in 1:m) {
        for (k in 1:m) {
          real th1 = tr_loc_pop     + tr_sd   * gh_x[i];
          real th2 = frac_logit_pop + frac_sd * gh_x[j];
          real th3 = init_logit_pop + init_sd * gh_x[k];
          beta_pts[q] = surrogate_anchor_betas(th1, th2, th3, Vinv, anchor_times);
          w[q] = gh_w[i] * gh_w[j] * gh_w[k];
          mean_full += w[q] * beta_pts[q];
          q += 1;
        }
      }
    }
    vector[2] mean12 = mean_full[2:3];
    matrix[2, 2] S = rep_matrix(0.0, 2, 2);
    for (r in 1:nq)
      S += w[r] * (beta_pts[r][2:3] - mean12) * (beta_pts[r][2:3] - mean12)';
    return (mean_full, S + diag_matrix(rep_vector(jitter, 2)));
  }

  real surrogate_ll(vector theta,
                    vector beta_pop,
                    real measure_sd, data vector log_lod_per_visit,
                    data int n_bg,
                    data vector bg_obs,
                    data array[] int bg_pos,
                    data array[] int bg_time) {
    int d = 2;
    real lp = 0;
    for (i in 1:n_bg) {
      real b0 = beta_pop[1];
      real b1 = beta_pop[2] + theta[(i - 1) * d + 1];
      real b2 = beta_pop[3] + theta[(i - 1) * d + 2];
      int v_start = bg_pos[i];
      int v_end   = bg_pos[i + 1] - 1;
      for (v in v_start:v_end) {
        real t = bg_time[v] - 1.0;
        real mu = b0 + b1 * t + b2 * t * t;
        if (bg_obs[v] > 0)
          lp += normal_lpdf(log(bg_obs[v]) | mu, measure_sd);
        else
          lp += normal_lcdf(log_lod_per_visit[v] | mu, measure_sd);
      }
    }
    return lp;
  }

  matrix surrogate_K_fn(matrix Sigma_beta, int n_bg) {
    int d = 2;
    matrix[n_bg * d, n_bg * d] K = rep_matrix(0, n_bg * d, n_bg * d);
    for (i in 1:n_bg) {
      int s = (i - 1) * d + 1;
      K[s:(s + d - 1), s:(s + d - 1)] = Sigma_beta;
    }
    return K;
  }
}

data {
  int<lower=1> n_patients;
  int<lower=1, upper=n_patients> n_forecast;   // first n_forecast patients = forecast
  int<lower=1> n_total_visits;
  vector[n_total_visits] normalized_obs;       // <= 0 marks below-LOD
  array[n_patients + 1] int patient_visit_pos;
  array[n_total_visits] int visit_time;        // 1-indexed week (t = idx - 1)
  real<lower=0> measure_sd;
  real log_lod;
  vector[3] anchor_times;
  int<lower=0, upper=1> mixed_mode;
}

transformed data {
  int n_background = n_patients - n_forecast;

  // Compact background-only views (patients ordered forecast-first).
  int bg_v_start = patient_visit_pos[n_forecast + 1];
  int bg_v_end   = patient_visit_pos[n_patients + 1] - 1;
  int n_bg_visits = bg_v_end - bg_v_start + 1;
  vector[n_bg_visits] bg_obs = normalized_obs[bg_v_start:bg_v_end];
  array[n_bg_visits] int bg_time;
  for (v in 1:n_bg_visits) bg_time[v] = visit_time[bg_v_start + v - 1];
  // Constant per-visit LOD vector (toy uses a single global log_lod). Mirrors
  // the production per-visit signature exactly: a constant vector reproduces the
  // validated PASS.
  vector[n_bg_visits] bg_log_lod = rep_vector(log_lod, n_bg_visits);
  array[n_background + 1] int bg_pos;
  for (j in 1:(n_background + 1))
    bg_pos[j] = patient_visit_pos[n_forecast + j] - bg_v_start + 1;

  int latent_dim = n_background * 2;
  vector[latent_dim] theta_0 = rep_vector(0.0, latent_dim);
  real tolerance = 1e-8;
  int max_num_steps = 100;
  int hessian_block_size = 2;
  int solver = 1;
  int max_steps_line_search = 0;
  int allow_fallback = 1;
  real jitter = 1e-10;

  matrix[3, 3] V;
  for (k in 1:3) {
    V[k, 1] = 1.0;
    V[k, 2] = anchor_times[k];
    V[k, 3] = anchor_times[k] * anchor_times[k];
  }
  matrix[3, 3] Vinv = inverse(V);

  vector[3] gh_x = [-sqrt(3.0), 0.0, sqrt(3.0)]';
  vector[3] gh_w = [1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0]';

  // In mode 0 (reference) ALL patients are full-HMC; in mode 1 only forecast are.
  int n_hmc = mixed_mode == 0 ? n_patients : n_forecast;
}

parameters {
  real tr_loc_pop;
  real<lower=0> tr_sd;
  real frac_logit_pop;
  real<lower=0> frac_sd;
  real init_logit_pop;
  real<lower=0> init_sd;

  vector[n_hmc] z_tr;
  vector[n_hmc] z_frac;
  vector[n_hmc] z_init;
}

model {
  tr_loc_pop     ~ normal(-3, 1);
  tr_sd          ~ normal(0, 1);
  frac_logit_pop ~ normal(0.5, 1);
  frac_sd        ~ normal(0, 1);
  init_logit_pop ~ normal(0, 1);
  init_sd        ~ normal(0, 1);

  // Full-HMC bi-exponential for the first n_hmc patients (all of them in mode 0,
  // just the forecast cohort in mode 1).
  z_tr   ~ std_normal();
  z_frac ~ std_normal();
  z_init ~ std_normal();
  for (i in 1:n_hmc) {
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

  // Backgrounded cohort via the surrogate (mode 1 only).
  if (mixed_mode == 1 && n_background > 0) {
    vector[3] beta_pop;
    matrix[2, 2] Sigma_beta;
    (beta_pop, Sigma_beta) = surrogate_bridge(
      tr_loc_pop, frac_logit_pop, init_logit_pop,
      tr_sd, frac_sd, init_sd, Vinv, anchor_times, gh_x, gh_w, jitter);

    target += laplace_marginal_tol(
      surrogate_ll,
      (beta_pop, measure_sd, bg_log_lod, n_background, bg_obs, bg_pos, bg_time),
      hessian_block_size,
      surrogate_K_fn,
      (Sigma_beta, n_background),
      (theta_0, tolerance, max_num_steps, solver, max_steps_line_search, allow_fallback)
    );
  }
}
