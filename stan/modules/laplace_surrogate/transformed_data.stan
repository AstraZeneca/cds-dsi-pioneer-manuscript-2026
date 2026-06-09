// --- Laplace surrogate: constants + validity guard --------------------------
// Constant inverse Vandermonde for the 3 fixed anchors (first anchor must be 0).
matrix[3, 3] surrogate_V;
for (k in 1:3) {
  surrogate_V[k, 1] = 1.0;
  surrogate_V[k, 2] = surrogate_anchor_times[k];
  surrogate_V[k, 3] = surrogate_anchor_times[k] * surrogate_anchor_times[k];
}
matrix[3, 3] surrogate_Vinv = inverse(surrogate_V);

// 3-point Gauss-Hermite nodes/weights for the standard normal (weights sum 1).
vector[3] surrogate_gh_x = [-sqrt(3.0), 0.0, sqrt(3.0)]';
vector[3] surrogate_gh_w = [1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0]';

// Laplace inner-solver knobs (linear-Gaussian surrogate => solver 1; may fall
// back to solver 2 on the collinear (t,t^2) basis — harmless, validated).
real surrogate_tolerance = 1e-8;
int surrogate_max_num_steps = 100;
int surrogate_hessian_block_size = 2;   // marginalize (beta_1, beta_2)
int surrogate_solver = 1;
int surrogate_max_steps_line_search = 0;
int surrogate_allow_fallback = 1;
real surrogate_jitter = 1e-10;

// GUARD: the bridge SD assumes the simple summed-variance patient SD, which
// holds for frac/init always and for tr UNLESS the tr SD sub-hierarchy is
// active (then the per-patient SD varies and the single bridge is wrong).
if (enable_background_surrogate == 1 && n_subhier_active_tr_intercept > 0)
  fatal_error("enable_background_surrogate=1 is incompatible with an active tr ",
              "SD sub-hierarchy (n_subhier_active_tr_intercept=",
              n_subhier_active_tr_intercept, "): the surrogate bridge assumes a ",
              "single patient-level tr SD. Disable the surrogate or the tr SD ",
              "sub-hierarchy.");

// First anchor must be 0 (the pinned intercept depends on g(0)=0).
if (enable_background_surrogate == 1 && surrogate_anchor_times[1] != 0.0)
  fatal_error("surrogate_anchor_times[1] must be 0 (baseline); got ",
              surrogate_anchor_times[1]);

vector[enable_background_surrogate == 1 ? n_background_patients * 2 : 0]
  surrogate_theta_0 = rep_vector(0.0,
    enable_background_surrogate == 1 ? n_background_patients * 2 : 0);

// Background-only compact views, assembled from transformed data only (so they
// satisfy the data-only argument qualifiers of surrogate_ll). Per-patient LOD
// offset (burden normalized to each patient's baseline => LOD shifts by
// log_baseline_sld[p]).
array[enable_background_surrogate == 1 ? n_background_patients + 1 : 0] int surrogate_bg_pos;
int surrogate_n_bg_visits = 0;
if (enable_background_surrogate == 1 && n_background_patients > 0) {
  surrogate_bg_pos[1] = 1;
  for (j in 1:n_background_patients) {
    int p = background_patient_idx[j];
    int vs, ve;
    (vs, ve) = get_pos(patient_visit_pos, p);
    surrogate_bg_pos[j + 1] = surrogate_bg_pos[j] + (ve - vs + 1);
  }
  surrogate_n_bg_visits = surrogate_bg_pos[n_background_patients + 1] - 1;
}
vector[surrogate_n_bg_visits] surrogate_bg_obs;
array[surrogate_n_bg_visits] int surrogate_bg_time;
vector[surrogate_n_bg_visits] surrogate_bg_log_lod;
if (enable_background_surrogate == 1 && n_background_patients > 0) {
  int w = 1;
  for (j in 1:n_background_patients) {
    int p = background_patient_idx[j];
    int vs, ve;
    (vs, ve) = get_pos(patient_visit_pos, p);
    for (v in vs:ve) {
      surrogate_bg_obs[w]     = normalized_sld[v];
      surrogate_bg_time[w]    = t_patient_visit_idx[v];
      surrogate_bg_log_lod[w] = log_lod - log_baseline_sld[p];
      w += 1;
    }
  }
}
