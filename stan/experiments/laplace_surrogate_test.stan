// Phase 1: log-concave QUADRATIC surrogate for backgrounded-trial patients,
// validated against the full bi-exponential via the built-in laplace_marginal_tol.
// ============================================================================
// Backgrounded patients' normalized log-burden is modeled as a quadratic in t:
//     log_burden_i(t) = beta_i0 + beta_i1*t + beta_i2*t^2 + N(0, measure_sd^2)
// Latents beta_i enter LINEARLY => Gaussian marginal => Laplace is EXACT
// => solver 1 (PD-Hessian Cholesky) is valid (unlike the bi-exponential, which
//    needed solver 3 just to survive).
//
// (beta_pop, Sigma_beta) are NOT free: they are the mechanistic-bridge image of
// the population rate params (tr_loc_pop, frac_logit_pop, init_logit_pop) and
// the per-patient SDs (tr_sd, frac_sd, init_sd):
//   beta_pop = Vinv * [g(t0), g(t1), g(t2)]   (g = exact bi-exponential log-burden)
//   Sigma_beta = J * diag(sd^2) * J'          (J = d beta_pop / d theta_pop, autodiff)
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

  // Quadratic coefficients (beta0, beta1, beta2) that interpolate the exact
  // bi-exponential log-burden at the 3 fixed anchor times. Vinv is the constant
  // inverse Vandermonde passed in as data.
  vector surrogate_anchor_betas(real tr_loc_pop, real frac_logit_pop,
                                real init_logit_pop,
                                data matrix Vinv, data vector anchor_times) {
    vector[3] g;
    for (k in 1:3)
      g[k] = bi_exp_log_burden(anchor_times[k], tr_loc_pop,
                               frac_logit_pop, init_logit_pop);
    return Vinv * g;
  }

  // Analytic Jacobian of the exact bi-exponential log-burden g(t) w.r.t. the
  // population location params (tr_loc, frac_logit, init_logit), at a single t.
  // g = log_sum_exp(a, b), a = init_log_dec - dec_rate*t, b = init_log_gro + gro_rate*t.
  // dg = w_dec*da + w_gro*db with w_dec = softmax weight on the decay branch.
  // (Stan has NO callable autodiff Jacobian of a user function — the `jacobian`
  //  block/`jacobian +=` are for custom transforms, not a returnable matrix — so
  //  we differentiate the closed form directly. It is elementary here.)
  row_vector bi_exp_log_burden_grad(real t, real tr_loc, real frac_logit,
                                    real init_logit) {
    real p = inv_logit(frac_logit);          // dec fraction
    real q = inv_logit(init_logit);           // initial dec share
    real log_dec_frac = log_inv_logit(frac_logit);
    real log_gro_frac = log1m_inv_logit(frac_logit);
    real dec_rate = exp(tr_loc + log_dec_frac);
    real gro_rate = exp(tr_loc + log_gro_frac);
    real a = log_inv_logit(init_logit) - dec_rate * t;
    real b = fmin(log1m_inv_logit(init_logit) + gro_rate * t, 500.0);
    real m = fmax(a, b);
    real w_dec = exp(a - m) / (exp(a - m) + exp(b - m));
    real w_gro = 1 - w_dec;
    // partials of a, b w.r.t. (tr_loc, frac_logit, init_logit)
    real da_dtr  = -t * dec_rate;             real db_dtr  =  t * gro_rate;
    real da_dfr  = -t * dec_rate * (1 - p);   real db_dfr  = -t * gro_rate * p;
    real da_din  = 1 - q;                     real db_din  = -q;
    row_vector[3] g_grad;
    g_grad[1] = w_dec * da_dtr + w_gro * db_dtr;
    g_grad[2] = w_dec * da_dfr + w_gro * db_dfr;
    g_grad[3] = w_dec * da_din + w_gro * db_din;
    return g_grad;
  }

  // Bridge covariance Sigma_beta = J diag(sd^2) J', J = d beta_pop / d theta_pop.
  // beta_pop = Vinv * g(anchors), so J = Vinv * [grad g(t_k)]_k (a 3x3 stack of
  // the per-anchor gradients). A diagonal jitter keeps K PD if an sd -> 0.
  matrix surrogate_bridge_cov(real tr_loc_pop, real frac_logit_pop,
                              real init_logit_pop,
                              real tr_sd, real frac_sd, real init_sd,
                              data matrix Vinv, data vector anchor_times,
                              data real jitter) {
    matrix[3, 3] g_jac;   // row k = d g(t_k) / d theta
    for (k in 1:3)
      g_jac[k] = bi_exp_log_burden_grad(anchor_times[k], tr_loc_pop,
                                        frac_logit_pop, init_logit_pop);
    matrix[3, 3] J = Vinv * g_jac;
    matrix[3, 3] D = diag_matrix(square([tr_sd, frac_sd, init_sd]'));
    return J * D * J' + diag_matrix(rep_vector(jitter, 3));
  }

  // laplace_marginal_tol functor. theta = stacked per-patient [b0,b1,b2] latents
  // in the WHITENED coordinate (K = identity, mean 0); we map to the bridged
  // distribution inside via the Cholesky of Sigma_beta passed through phi.
  // Here we use the simpler route: K_fn returns Sigma_beta directly, so theta is
  // in the natural beta coordinate centered at beta_pop. So the functor receives
  // beta_pop and adds (theta_i - 0) usage: latents are deviations from beta_pop.
  real surrogate_ll(vector theta,
                    vector beta_pop,
                    real measure_sd, real log_lod,
                    data int n_patients,
                    data vector normalized_obs,
                    data array[] int patient_visit_pos,
                    data array[] int visit_time,
                    data array[] int patient_of_visit) {
    int d = 3;
    real lp = 0;
    for (i in 1:n_patients) {
      real b0 = beta_pop[1] + theta[(i - 1) * d + 1];
      real b1 = beta_pop[2] + theta[(i - 1) * d + 2];
      real b2 = beta_pop[3] + theta[(i - 1) * d + 3];
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

  // Prior covariance functor for laplace_marginal_tol: the bridged Sigma_beta,
  // block-replicated per patient is handled by hessian_block_size=3 + this K
  // returning the per-block covariance tiled. We return a full block-diagonal.
  matrix surrogate_K_fn(matrix Sigma_beta, int n_patients) {
    int d = 3;
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
  int latent_dim = n_patients * 3;
  vector[latent_dim] theta_0 = rep_vector(0.0, latent_dim);
  real tolerance = 1e-8;
  int max_num_steps = 100;        // log-concave => Newton converges fast
  int hessian_block_size = 3;
  int solver = 1;                 // PD-Hessian Cholesky: valid for log-concave
  int max_steps_line_search = 0;  // not needed when well-conditioned
  int allow_fallback = 1;
  real jitter = 1e-8;

  // Constant inverse Vandermonde for the 3 fixed anchors.
  matrix[3, 3] V;
  for (k in 1:3) {
    V[k, 1] = 1.0;
    V[k, 2] = anchor_times[k];
    V[k, 3] = anchor_times[k] * anchor_times[k];
  }
  matrix[3, 3] Vinv = inverse(V);
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
    matrix[3, 3] Sigma_beta = surrogate_bridge_cov(
      tr_loc_pop, frac_logit_pop, init_logit_pop,
      tr_sd, frac_sd, init_sd, Vinv, anchor_times, jitter);

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
