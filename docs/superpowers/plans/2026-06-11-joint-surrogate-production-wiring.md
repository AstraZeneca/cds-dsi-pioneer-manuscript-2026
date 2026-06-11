# Joint SLD+MS Laplace Surrogate — Production Wiring Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Wire the validated joint SLD+multistate Laplace surrogate (with correlated patient-level frailty) into the publication model `stan/tumor/sf-ssm-log-space.stan`, so ~419 background-trial patients contribute to population/coupling parameters via marginalized latents instead of full trajectory simulation.

**Architecture:** Extend the existing `stan/modules/laplace_surrogate/` module in place from d=2 (SLD-only) to d=2+n_frailty (joint SLD+MS). The functor gains visit-gated 0→1 + continuous 0→3 hazard terms whose burden feature is built with the SAME `standardize_log_burden()` helper the forecast path uses (coordinate-frame identity is mandatory — `tv_coef` is shared model-wide). A d=2 Stan overload preserves the LFO consumer (`sf-ssls-lfo.stan`) numerically. Per-patient prior covariance is `K = blockdiag(Σ_β, Σ_u)`.

**Tech Stack:** Stan (cmdstan 2.39, `laplace_marginal_tol`), R (cmdstanr driver for the feature-scale check), `stanc` for compile gates.

**Spec:** `docs/superpowers/specs/2026-06-11-joint-surrogate-production-wiring-design.md` (reviewed `wf_922c4cd6`).
**Math template:** `stan/experiments/laplace_joint_frailty_test.stan` (validated d=4 gate functor `joint_frailty_ll`).

---

## File Structure

| File | Responsibility | Action |
|---|---|---|
| `stan/modules/laplace_surrogate/surrogate.stanfunctions` | The joint functor (parametric) + d=2 overload + `surrogate_K_fn` (parametric) + d=2 overload | Modify |
| `stan/modules/laplace_surrogate/transformed_data.stan` | `surrogate_d` / `frailty_slots` detection, guards, bg-compact MS arrays (event/censor/visit mask, baseline_week) | Modify |
| `stan/modules/laplace_surrogate/likelihood.stan` | Assemble static baselines + Σ_u + K; branch d=2 (overload) vs joint | Modify |
| `stan/tumor/sf-ssm-log-space.stan` | Host model: 4 includes + line-124 flag fix | Modify |
| `r/experiments/test_surrogate_feature_scale.R` | Feature-scale equivalence check (surrogate vs forecast frame) | Create |

**Interdependency note:** Tasks 1–4 edit `#include` fragments + the host. Stan fragments do NOT compile in isolation, and the functor signature + its single call site must change together. The integration compile gate is **Task 5** — it validates Tasks 1–4 as a unit. Do not expect a green compile before Task 5.

---

## Task 1: Detect `d` / frailty slots, build bg MS arrays + guards (`transformed_data.stan`)

**Files:**
- Modify: `stan/modules/laplace_surrogate/transformed_data.stan` (append after the existing line 76 SLD-view block)

Today the module's `transformed_data.stan` ends at line 76 with the SLD-compact views (`surrogate_bg_obs/time/log_lod`) and a hardcoded `surrogate_hessian_block_size = 2` (line 19). This task adds the d-parametric detection, the guards, and the bg-compact MS arrays the joint functor needs.

- [ ] **Step 1: Replace the hardcoded `surrogate_hessian_block_size = 2` and add d-detection**

In `stan/modules/laplace_surrogate/transformed_data.stan`, change line 19 from:

```stan
int surrogate_hessian_block_size = 2;   // marginalize (beta_1, beta_2)
```

to:

```stan
// d-parametric latent dim: 2 burden REs + one frailty latent per frailty slot.
// Detected from the MS intercept-mode config (patient level = last hierarchy col).
int surrogate_patient_lv = n_levels;
// Coupled slots = transitions whose hazard reads the burden coupling. Publication
// = {01, 03}; both read ms_time_varying_covar_01 (transformed_parameters.stan:176,780).
// We detect frailty slots as coupled slots with a patient-level RE-NCP intercept.
int surrogate_n_frailty_slots = 0;
array[2] int surrogate_frailty_slot = {0, 0};   // ascending slot order; 0 = unused
if (enable_background_surrogate == 1) {
  // MS_SLOT_01 / MS_SLOT_03 are in scope (modules/multistate/transformed_data.stan).
  if (enable_ms_01 && ms_level_intercept_mode[MS_SLOT_01, surrogate_patient_lv] == 2) {
    surrogate_n_frailty_slots += 1;
    surrogate_frailty_slot[surrogate_n_frailty_slots] = MS_SLOT_01;
  }
  if (enable_ms_03 && ms_level_intercept_mode[MS_SLOT_03, surrogate_patient_lv] == 2) {
    surrogate_n_frailty_slots += 1;
    surrogate_frailty_slot[surrogate_n_frailty_slots] = MS_SLOT_03;
  }
}
int surrogate_d = 2 + surrogate_n_frailty_slots;
int surrogate_hessian_block_size = surrogate_d;
```

- [ ] **Step 2: Add the guards (after the d-detection block)**

Append immediately after Step 1's block. (`surrogate_frailty_block` is declared
here at file scope so `likelihood.stan` can read it in Task 3.)

```stan
int surrogate_frailty_block = 0;   // patient-level MS correlation block index (set by guard 4)
if (enable_background_surrogate == 1) {
  // Guard 6: split must be at a trial (non-patient) level — the baseline gather
  // depends on it (a patient-level split falls into the ->1 fallback).
  if (!(forecast_split_level > 0 && forecast_split_level < n_levels))
    fatal_error("enable_background_surrogate=1 requires 0 < forecast_split_level < ",
                "n_levels (got ", forecast_split_level, "); the surrogate baseline ",
                "gather is only valid for a trial-level split.");

  // Guards 1 & 2: every coupled/frailty slot's patient-level baseline must be a
  // marginalizable Gaussian intercept (mode 0 or 2) with NO patient-level GP.
  array[2] int surrogate_coupled_slot = {MS_SLOT_01, MS_SLOT_03};
  for (s_idx in 1:2) {
    int slot = surrogate_coupled_slot[s_idx];
    int slot_on = (slot == MS_SLOT_01) ? enable_ms_01 : enable_ms_03;
    if (slot_on) {
      if (enable_ms_level_gp[slot, surrogate_patient_lv] != 0)
        fatal_error("enable_background_surrogate=1 forbids a patient-level GP on ",
                    "coupled MS slot ", slot, " (not a marginalizable scalar).");
      int mode = ms_level_intercept_mode[slot, surrogate_patient_lv];
      if (!(mode == 0 || mode == 2))
        fatal_error("enable_background_surrogate=1: coupled MS slot ", slot,
                    " patient-level intercept mode must be 0 (off) or 2 (RE-NCP); ",
                    "got ", mode, ". Patient-level baseline beyond marginalized ",
                    "Gaussian frailty is unsupported by the surrogate.");
    }
  }

  // Guard 7: the joint functor accesses tv_coef_01[2] / tv_coef_03[2] (level +
  // velocity). time_varying_coef_01 is sized 1 when
  // (enable_ms_visit_gated_01 && !enable_ms_visit_gated_latent_01), else
  // n_time_varying_covar (parameters.stan:44). The surrogate requires the 2-feature
  // (level, velocity) basis => assert n_time_varying_covar == 2 AND not the size-1
  // observed-gating case. Publication: enable_ms_visit_gated_latent_01 = TRUE => size 2 (OK).
  if (n_time_varying_covar != 2)
    fatal_error("enable_background_surrogate=1 requires the 2-feature (level, ",
                "velocity) coupling basis (n_time_varying_covar=2); got ",
                n_time_varying_covar, ".");
  if (enable_ms_visit_gated_01 && !enable_ms_visit_gated_latent_01)
    fatal_error("enable_background_surrogate=1 is incompatible with OBSERVED ",
                "visit-gated 0->1 (enable_ms_visit_gated_latent_01=0), which sizes ",
                "time_varying_coef_01 to 1; the surrogate needs the latent 2-feature basis.");
  // Guard 7b: the functor unconditionally reads tv_coef_01[1..2] / tv_coef_03[1..2]
  // whenever a slot is a frailty slot. time_varying_coef_01 is size 0 unless
  // enable_ms_pop_time_varying_cov (parameters.stan:43); time_varying_coef_03 size 0
  // unless enable_ms_pop_time_varying_cov && enable_ms_03_time_varying_cov
  // (parameters.stan:149). Require coupling ON for every detected frailty slot
  // (frailty_slots subset of coupled-TV slots, per spec §2).
  if (surrogate_n_frailty_slots > 0 && !enable_ms_pop_time_varying_cov)
    fatal_error("enable_background_surrogate=1 with frailty slots requires ",
                "enable_ms_pop_time_varying_cov=1 (else time_varying_coef_01 is size 0 ",
                "and the functor reads tv_coef_01[2] out of bounds).");
  for (m in 1:surrogate_n_frailty_slots)
    if (surrogate_frailty_slot[m] == MS_SLOT_03 && !enable_ms_03_time_varying_cov)
      fatal_error("enable_background_surrogate=1 with 0->3 frailty requires ",
                  "enable_ms_03_time_varying_cov=1 (else time_varying_coef_03 is size 0).");

  // Guard 4: the patient-level correlation block's member slots must equal the
  // detected frailty_slots in the same (ascending) order, so Sigma_u rows align
  // with the theta frailty entries. Also record the block index for likelihood.stan.
  if (surrogate_n_frailty_slots >= 2) {
    for (b in 1:n_ms_corr_blocks) {
      if (ms_corr_block_level[b] == surrogate_patient_lv) {
        surrogate_frailty_block = b;
        for (m in 1:surrogate_n_frailty_slots)
          if (ms_corr_block_member_slot[b, m] != surrogate_frailty_slot[m])
            fatal_error("enable_background_surrogate=1: patient-level correlation ",
                        "block member slot order does not match detected frailty ",
                        "slots; Sigma_u would misalign with the latent vector.");
      }
    }
    if (surrogate_frailty_block == 0)
      fatal_error("enable_background_surrogate=1: expected a patient-level MS ",
                  "correlation block for ", surrogate_n_frailty_slots,
                  " frailty slots, found none.");
  }
}
```

- [ ] **Step 3: Resize `surrogate_theta_0` to d-parametric length**

Change the existing `surrogate_theta_0` declaration (lines 40-42) from:

```stan
vector[enable_background_surrogate == 1 ? n_background_patients * 2 : 0]
  surrogate_theta_0 = rep_vector(0.0,
    enable_background_surrogate == 1 ? n_background_patients * 2 : 0);
```

to:

```stan
vector[enable_background_surrogate == 1 ? n_background_patients * surrogate_d : 0]
  surrogate_theta_0 = rep_vector(0.0,
    enable_background_surrogate == 1 ? n_background_patients * surrogate_d : 0);
```

- [ ] **Step 4: Build the bg-compact MS arrays (append at end of file)**

Append after the existing SLD-view block (current line 76). These mirror the gate's `visit_wk` / `ms_event_wk_01/03` / `ms_censored_01/03` (`r/experiments/test_laplace_joint_frailty.R:47-95`) but are derived from the unified production MS data fields. `n_wk` for the surrogate is `max_all_t`.

```stan
// --- bg-compact MS arrays for the joint surrogate hazard terms ---------------
// Built only when frailty/coupling is active; otherwise zero-sized (d=2 path).
int surrogate_ms_active = (enable_background_surrogate == 1 && surrogate_n_frailty_slots > 0) ? 1 : 0;
int surrogate_n_wk = surrogate_ms_active ? max_all_t : 0;

// Per-bg event week + censor flag for 0->1 and 0->3, and the per-bg baseline week
// (calendar week of the patient's baseline visit) used to anchor g(tau).
array[surrogate_ms_active ? n_background_patients : 0] int surrogate_bg_event_wk_01;
array[surrogate_ms_active ? n_background_patients : 0] int surrogate_bg_censored_01;
array[surrogate_ms_active ? n_background_patients : 0] int surrogate_bg_event_wk_03;
array[surrogate_ms_active ? n_background_patients : 0] int surrogate_bg_censored_03;
array[surrogate_ms_active ? n_background_patients : 0] int surrogate_bg_baseline_week;
// 0->1 visit-gating mask: [bg, week] == 1 at observed visit weeks.
array[surrogate_ms_active ? n_background_patients : 0,
      surrogate_ms_active ? max_all_t : 0] int surrogate_bg_visit_wk_01
  = rep_array(0, surrogate_ms_active ? n_background_patients : 0,
                 surrogate_ms_active ? max_all_t : 0);

if (surrogate_ms_active) {
  for (j in 1:n_background_patients) {
    int p = background_patient_idx[j];
    // 0->1 event/censor: ms_censored_01 is a raw data field.
    surrogate_bg_censored_01[j] = ms_censored_01[p];
    surrogate_bg_event_wk_01[j] = ms_censored_01[p] ? max_all_t : ms_time_01[p];
    // 0->3 event/censor: derived. ms_final_state == 3 => dropout event at ms_time_03;
    // otherwise censored at ms_time_03 (= patient_max_t).
    surrogate_bg_censored_03[j] = (ms_final_state[p] == 3) ? 0 : 1;
    surrogate_bg_event_wk_03[j] = ms_time_03[p];
    // Baseline week = calendar week of the patient's BASELINE visit = the LAST
    // screening visit, NOT the first visit. Matches _full_model_transformed_data.stan:20-21
    // (baseline_visit_idx = curr_patient_visit_pos + n_patient_screening_visits[i] - 1),
    // which is the anchor t_patient_visit_idx uses for the SLD term. Using the first
    // visit here would mis-anchor g(tau) by n_screening-1 weeks and break the
    // shared-frame identity between the SLD and hazard terms. n_patient_screening_visits
    // is in scope from _visit_transformed_data.stan:16 (included before this module).
    int vs, ve;
    (vs, ve) = get_pos(patient_visit_pos, p);
    surrogate_bg_baseline_week[j] = t_patient_visits[vs + n_patient_screening_visits[p] - 1];
    // Visit-gating mask: mark each observed visit's calendar week.
    for (v in vs:ve) {
      int wk = t_patient_visits[v];
      if (wk >= 1 && wk <= max_all_t) surrogate_bg_visit_wk_01[j, wk] = 1;
    }
  }
}
```

- [ ] **Step 5: Syntax sanity (deferred to Task 5 compile gate)**

This fragment cannot compile standalone. Record that it will be validated by the Task 5 `stanc` gate. Do a manual read-through: confirm every symbol referenced (`n_levels`, `MS_SLOT_01/03`, `enable_ms_01/03`, `ms_level_intercept_mode`, `enable_ms_level_gp`, `n_ms_corr_blocks`, `ms_corr_block_level`, `ms_corr_block_member_slot`, `forecast_split_level`, `background_patient_idx`, `ms_censored_01`, `ms_time_01`, `ms_final_state`, `ms_time_03`, `patient_visit_pos`, `t_patient_visits`, `max_all_t`, `get_pos`) is in scope at the module's include point in `sf-ssm-log-space.stan` (after `modules/multistate/transformed_data.stan`).

- [ ] **Step 6: Commit**

```bash
git add stan/modules/laplace_surrogate/transformed_data.stan
git commit -m "feat(laplace): d-parametric detection + guards + bg MS arrays in surrogate transformed_data"
```

---

## Task 2: Parametric joint functor + d=2 overload (`surrogate.stanfunctions`)

**Files:**
- Modify: `stan/modules/laplace_surrogate/surrogate.stanfunctions` (replace `surrogate_ll` at lines 81-106 and `surrogate_K_fn` at lines 108-116)

The current `surrogate_ll` (lines 81-106) is SLD-only with `int d = 2`. Replace it with a parametric joint functor that adds the visit-gated 0→1 + continuous 0→3 hazard terms, PLUS a d=2 overload that forwards with degenerate args so the LFO call site is unchanged.

- [ ] **Step 1: Replace `surrogate_ll` (lines 81-106) with the parametric joint functor**

```stan
// JOINT SLD + multistate(0->1 visit-gated, 0->3 continuous) functor with
// correlated patient-level frailty. Ported from the validated gate
// (laplace_joint_frailty_test.stan:51-121). Latents enter RAW; scale lives in K.
// Burden feature is built via standardize_log_burden() — the SAME helper the
// forecast path uses (_ms_burden_tv_covar.stan:64-65) — so the shared tv_coef
// applies an identical-frame coupling to background and forecast patients.
// theta block per patient (length d): [b1, b2, (u01), (u03)] (frailty optional).
real surrogate_ll(vector theta,
                  vector beta_pop,                     // [b0, b1, b2] (b0 pinned)
                  real measure_sd, data vector log_lod_per_visit,
                  data int n_bg,
                  data vector bg_obs,
                  data array[] int bg_pos,
                  data array[] int bg_time,
                  // --- joint MS extension (degenerate/empty for the d=2 path) ---
                  data int d,
                  data int n_frailty_slots,
                  data int n_wk,
                  vector base01_static_flat,           // [n_bg * n_wk], row-major per bg patient
                  vector base03_static_flat,
                  vector log_baseline_burden_bg,       // [n_bg]
                  vector tv_coef_01,                   // [2] (level, velocity); empty if 0->1 not coupled
                  vector tv_coef_03,                   // [2]; empty if 0->3 not coupled
                  data real median_log_burden_obs, data real iqr_log_burden_obs,
                  data real median_velocity_obs, data real iqr_velocity_obs,
                  data array[] int bg_event_wk_01, data array[] int bg_censored_01,
                  data array[] int bg_event_wk_03, data array[] int bg_censored_03,
                  data array[] int bg_baseline_week,
                  data array[,] int bg_visit_wk_01) {
  real b0 = beta_pop[1];
  real lp = 0;
  for (i in 1:n_bg) {
    real b1 = beta_pop[2] + theta[(i - 1) * d + 1];
    real b2 = beta_pop[3] + theta[(i - 1) * d + 2];

    // --- SLD term (anchored time t = bg_time[v] - 1; baseline -> t=0 -> g(0)=b0) ---
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

    if (n_frailty_slots > 0) {
      // Frailty latents follow the 2 burden latents, ascending slot order [u01,u03].
      real u01 = theta[(i - 1) * d + 3];
      real u03 = (n_frailty_slots >= 2) ? theta[(i - 1) * d + 4] : 0.0;
      int bw = bg_baseline_week[i];

      // --- 0->1 visit-gated hazard ---
      int te1 = bg_censored_01[i] ? n_wk : bg_event_wk_01[i];
      real cum_haz_01 = 0;
      for (w in 1:te1) {
        if (bg_visit_wk_01[i, w] == 1) {
          real tau = w - bw;                            // anchored time (same frame as g fit)
          real g_norm = b0 + b1 * tau + b2 * tau * tau; // NORMALIZED log-burden
          real lvl_feat = standardize_log_burden(g_norm, log_baseline_burden_bg[i],
                                                  median_log_burden_obs, iqr_log_burden_obs);
          real vel_feat = standardize_velocity(b1 + 2 * b2 * tau,
                                                median_velocity_obs, iqr_velocity_obs);
          real loghaz = base01_static_flat[(i - 1) * n_wk + w] + u01
                      + tv_coef_01[1] * lvl_feat + tv_coef_01[2] * vel_feat;
          cum_haz_01 += exp(loghaz);
          if (!bg_censored_01[i] && w == te1) lp += loghaz;
        }
      }
      lp += -cum_haz_01;

      // --- 0->3 continuous hazard (every week) ---
      int te3 = bg_censored_03[i] ? n_wk : bg_event_wk_03[i];
      real cum_haz_03 = 0;
      for (w in 1:te3) {
        real tau = w - bw;
        real g_norm = b0 + b1 * tau + b2 * tau * tau;
        real lvl_feat = standardize_log_burden(g_norm, log_baseline_burden_bg[i],
                                                median_log_burden_obs, iqr_log_burden_obs);
        real vel_feat = standardize_velocity(b1 + 2 * b2 * tau,
                                              median_velocity_obs, iqr_velocity_obs);
        real loghaz = base03_static_flat[(i - 1) * n_wk + w] + u03
                    + tv_coef_03[1] * lvl_feat + tv_coef_03[2] * vel_feat;
        cum_haz_03 += exp(loghaz);
        if (!bg_censored_03[i] && w == te3) lp += loghaz;
      }
      lp += -cum_haz_03;
    }
  }
  return lp;
}

// d=2 SLD-only OVERLOAD: the original signature LFO #includes. Forwards to the
// parametric body with degenerate (empty/zero) MS inputs => numerically the
// original quadratic-surrogate SLD likelihood. n_frailty_slots=0 => hazard terms
// are skipped entirely.
real surrogate_ll(vector theta,
                  vector beta_pop,
                  real measure_sd, data vector log_lod_per_visit,
                  data int n_bg,
                  data vector bg_obs,
                  data array[] int bg_pos,
                  data array[] int bg_time) {
  array[0] int empty_i;
  array[0, 0] int empty_mask;
  return surrogate_ll(theta, beta_pop, measure_sd, log_lod_per_visit, n_bg,
                      bg_obs, bg_pos, bg_time,
                      2, 0, 0,
                      rep_vector(0.0, 0), rep_vector(0.0, 0), rep_vector(0.0, 0),
                      rep_vector(0.0, 0), rep_vector(0.0, 0),
                      1.0, 1.0, 0.0, 1.0,
                      empty_i, empty_i, empty_i, empty_i, empty_i, empty_mask);
}
```

- [ ] **Step 2: Replace `surrogate_K_fn` (lines 108-116) with parametric + overload**

```stan
// Parametric K tiler: per-patient block = blockdiag(Sigma_beta[2x2], Sigma_u).
// Sigma_u is (d-2)x(d-2); for d=2 it is empty and K is the pure burden tiling.
matrix surrogate_K_fn(matrix Sigma_beta, matrix Sigma_u, int n_bg, data int d) {
  matrix[n_bg * d, n_bg * d] K = rep_matrix(0, n_bg * d, n_bg * d);
  int nf = d - 2;
  for (i in 1:n_bg) {
    int s = (i - 1) * d + 1;
    K[s:(s + 1), s:(s + 1)] = Sigma_beta;                       // burden block
    if (nf > 0)
      K[(s + 2):(s + d - 1), (s + 2):(s + d - 1)] = Sigma_u;    // frailty block
  }
  return K;
}

// d=2 OVERLOAD (LFO call site): no frailty block.
matrix surrogate_K_fn(matrix Sigma_beta, int n_bg) {
  return surrogate_K_fn(Sigma_beta, rep_matrix(0.0, 0, 0), n_bg, 2);
}
```

- [ ] **Step 3: Manual read-through for log-concavity + symbol scope**

Confirm: (a) `loghaz` is affine in `(b1,b2,u01,u03)` — `g_norm` and `b1+2·b2·τ` are linear, `standardize_*` are affine, `base*_static_flat` and `u*` are additive constants/latents; (b) `standardize_log_burden` / `standardize_velocity` are defined in `_burden.stanfunctions` (included before this file in the host's `functions` block); (c) the flattened index `(i-1)*n_wk + w` matches the row-major flattening the model block will use in Task 3.

- [ ] **Step 4: Commit**

```bash
git add stan/modules/laplace_surrogate/surrogate.stanfunctions
git commit -m "feat(laplace): parametric joint surrogate_ll + surrogate_K_fn with d=2 overloads"
```

---

## Task 3: Assemble baselines + Σ_u + K and branch the call (`likelihood.stan`)

**Files:**
- Modify: `stan/modules/laplace_surrogate/likelihood.stan` (extend the `target +=` block at lines 32-41)

The current `likelihood.stan` computes `bg_beta_pop` / `bg_Sigma` via `surrogate_bridge` then calls `laplace_marginal_tol` with the d=2 tuples. This task assembles the static baselines, Σ_u, and the per-bg `log_baseline_burden`, then branches: joint call when `surrogate_n_frailty_slots > 0`, else the existing d=2 overload call.

- [ ] **Step 1: Insert the joint-input assembly before the `target +=` call**

In `stan/modules/laplace_surrogate/likelihood.stan`, after the `surrogate_bridge(...)` call (current line 30) and before `target += laplace_marginal_tol(` (current line 32), insert:

```stan
  if (surrogate_n_frailty_slots > 0) {
    // --- per-bg static baselines (theta-independent), row-major flattened ---
    vector[n_background_patients * surrogate_n_wk] base01_static_flat;
    vector[n_background_patients * surrogate_n_wk] base03_static_flat;
    vector[n_background_patients] log_baseline_burden_bg;
    // TI-cov linear predictors in QR space on the BACKGROUND rows (mirror
    // transformed_parameters.stan:216,787). Publication: pop-level term only
    // (enable_ms_level_cov both FALSE), so no level-slope path. GUARD with the
    // SAME conditions the forecast path uses — time_invariant_coef_qr_01/03 are
    // vector[0] when their flags are off (parameters.stan:47,154), so an
    // unconditional multiply would be a (n_bg × n_covar)·vector[0] mismatch.
    vector[n_background_patients] linpred_bg_01 = zeros_vector(n_background_patients);
    vector[n_background_patients] linpred_bg_03 = zeros_vector(n_background_patients);
    if (enable_ms_pop_time_invariant_cov && n_time_invariant_covar > 0)
      linpred_bg_01 = Q_covar_design_matrix[background_patient_idx, ] * time_invariant_coef_qr_01;
    if (enable_ms_pop_time_invariant_cov && enable_ms_03_time_invariant_cov && n_time_invariant_covar > 0)
      linpred_bg_03 = Q_covar_design_matrix[background_patient_idx, ] * time_invariant_coef_qr_03;
    for (j in 1:n_background_patients) {
      int p = background_patient_idx[j];
      log_baseline_burden_bg[j] = log_baseline_burden[p];
      row_vector[surrogate_n_wk] b01 = log_pop_lambda_01;
      row_vector[surrogate_n_wk] b03 = log_pop_lambda_03;
      for (lv in 1:n_levels) {
        if (lv < n_levels) {  // trial-level baseline residual (patient level is frailty, in theta)
          if (enable_ms_01)
            b01 += log_level_lambda_01_residual[patient_ms_baseline_flat_idx_slot[MS_SLOT_01, p, lv]];
          if (enable_ms_03)
            b03 += log_level_lambda_03_residual[patient_ms_baseline_flat_idx_slot[MS_SLOT_03, p, lv]];
        }
      }
      b01 += linpred_bg_01[j];
      b03 += linpred_bg_03[j];
      for (w in 1:surrogate_n_wk) {
        base01_static_flat[(j - 1) * surrogate_n_wk + w] = b01[w];
        base03_static_flat[(j - 1) * surrogate_n_wk + w] = b03[w];
      }
    }

    // --- Sigma_u: diag(sigma) * (L L') * diag(sigma) (matches ms_corr_u, tp:37) ---
    vector[surrogate_n_frailty_slots] sigma_u;
    sigma_u[1] = (surrogate_frailty_slot[1] == MS_SLOT_01)
               ? log_lambda_gp_01_level_intercept_sd[surrogate_patient_lv]
               : log_lambda_gp_03_level_intercept_sd[surrogate_patient_lv];
    if (surrogate_n_frailty_slots >= 2)
      sigma_u[2] = log_lambda_gp_03_level_intercept_sd[surrogate_patient_lv];
    matrix[surrogate_n_frailty_slots, surrogate_n_frailty_slots] Lu =
      diag_pre_multiply(sigma_u, L_ms_intercept_corr[surrogate_frailty_block]);
    matrix[surrogate_n_frailty_slots, surrogate_n_frailty_slots] Sigma_u =
      Lu * Lu' + diag_matrix(rep_vector(surrogate_jitter, surrogate_n_frailty_slots));

    target += laplace_marginal_tol(
      surrogate_ll,
      (bg_beta_pop, measure_sd_sld, surrogate_bg_log_lod, n_background_patients,
       surrogate_bg_obs, surrogate_bg_pos, surrogate_bg_time,
       surrogate_d, surrogate_n_frailty_slots, surrogate_n_wk,
       base01_static_flat, base03_static_flat, log_baseline_burden_bg,
       time_varying_coef_01, time_varying_coef_03,
       median_log_burden_obs, iqr_log_burden_obs, median_velocity_obs, iqr_velocity_obs,
       surrogate_bg_event_wk_01, surrogate_bg_censored_01,
       surrogate_bg_event_wk_03, surrogate_bg_censored_03,
       surrogate_bg_baseline_week, surrogate_bg_visit_wk_01),
      surrogate_hessian_block_size,
      surrogate_K_fn,
      (bg_Sigma, Sigma_u, n_background_patients, surrogate_d),
      (surrogate_theta_0, surrogate_tolerance, surrogate_max_num_steps,
       surrogate_solver, surrogate_max_steps_line_search, surrogate_allow_fallback)
    );
  } else {
```

- [ ] **Step 2: Close the branch around the existing d=2 call**

The existing `target += laplace_marginal_tol(...)` (lines 32-41) becomes the `else` body. Indent it under the `else {` opened in Step 1 and close the brace. The existing call uses the d=2 overloads (`surrogate_ll` 8-arg, `surrogate_K_fn` 2-arg), so it is unchanged except for indentation and the closing `}`:

```stan
    target += laplace_marginal_tol(
      surrogate_ll,
      (bg_beta_pop, measure_sd_sld, surrogate_bg_log_lod, n_background_patients,
       surrogate_bg_obs, surrogate_bg_pos, surrogate_bg_time),
      surrogate_hessian_block_size,
      surrogate_K_fn,
      (bg_Sigma, n_background_patients),
      (surrogate_theta_0, surrogate_tolerance, surrogate_max_num_steps,
       surrogate_solver, surrogate_max_steps_line_search, surrogate_allow_fallback)
    );
  }
```

(Note: `surrogate_hessian_block_size == surrogate_d == 2` in this branch, so the d=2 overload's K is correct.)

- [ ] **Step 3: Note on `surrogate_frailty_block`**

`surrogate_frailty_block` is already declared and set in Task 1 (file-scope int,
assigned by guard 4). `likelihood.stan` reads it only inside the
`surrogate_n_frailty_slots > 0` branch. Caveat for a *single*-frailty-slot config
(not the publication, which has 2): `L_ms_intercept_corr` is
`array[n_ms_corr_blocks] cholesky_factor_corr[ms_corr_dim]` and a correlation block
only forms when `ms_corr_dim >= 2` — so a 1-slot frailty case would have
`surrogate_frailty_block == 0` and the `Lu` line would need a 1×1 identity fallback.
The publication is 2-slot, so this path is exercised correctly; flag the 1-slot case
as a compile-time `fatal_error` if it ever arises rather than silently mis-indexing.

- [ ] **Step 4: Manual read-through for symbol scope**

Confirm in scope at the `likelihood.stan` include point (model block, after `modules/multistate/transformed_parameters.stan`): `log_pop_lambda_01/03`, `log_level_lambda_01/03_residual`, `patient_ms_baseline_flat_idx_slot`, `Q_covar_design_matrix`, `time_invariant_coef_qr_01/03`, `time_varying_coef_01/03`, `log_baseline_burden`, `log_lambda_gp_01/03_level_intercept_sd`, `L_ms_intercept_corr`, `median_log_burden_obs`, `iqr_log_burden_obs`, `median_velocity_obs`, `iqr_velocity_obs`, plus all `surrogate_*` transformed-data symbols from Task 1.

- [ ] **Step 5: Commit**

```bash
git add stan/modules/laplace_surrogate/likelihood.stan stan/modules/laplace_surrogate/transformed_data.stan
git commit -m "feat(laplace): assemble static baselines + Sigma_u + branch joint vs d=2 surrogate call"
```

---

## Task 4: Wire the module into `sf-ssm-log-space.stan` + fix line-124 flag

**Files:**
- Modify: `stan/tumor/sf-ssm-log-space.stan` (functions/data/transformed-data/model blocks)

- [ ] **Step 1: Add the functions include**

In the `functions` block (after line 13 `#include "modules/tumor/tumor.stanfunctions"`), add:

```stan
  #include "modules/laplace_surrogate/surrogate.stanfunctions"
```

- [ ] **Step 2: Add the data includes**

In the `data` block, after line 33 (`#include "modules/init/flags.stan"`), add:

```stan
  #include "modules/laplace_surrogate/flags.stan"
  #include "modules/laplace_surrogate/data.stan"
```

- [ ] **Step 3: Add the transformed-data include**

In the `transformed data` block, after line 53 (`#include "modules/multistate/transformed_data.stan"`), add:

```stan
  #include "modules/laplace_surrogate/transformed_data.stan"
```

(Must come AFTER multistate transformed_data — needs `MS_SLOT_*`, `ms_level_intercept_mode`, `n_ms_corr_blocks`; AFTER `_full_model_transformed_data.stan` at line 46 — needs `background_patient_idx`; AFTER `_tumor_observed_covar_transformed_data.stan` would be ideal for the velocity constants but those are only read in `likelihood.stan`, not here. The constants are at line 54 which is after — fine.)

- [ ] **Step 4: Add the model-block likelihood include**

The model-block brace structure is: line 124 `0` (the flag arg), 125 `);`, 126 `}` (close `profile("multistate loglik")`), 127 `}` (close `if (fit_multistate_data)`), 128 `}` (close `profile("loglik")`), 129 `}` (close `model`). Place the include **between line 127 and 128** — inside `profile("loglik")`, after the `if (fit_multistate_data)` block closes, so it runs whenever the model block runs and the surrogate's own `if (enable_background_surrogate == 1 && n_background_patients > 0)` guard (already in `likelihood.stan:5`) gates it:

```stan
    }   // <- line 127: closes if (fit_multistate_data)

    #include "modules/laplace_surrogate/likelihood.stan"
  }     // <- closes profile("loglik")
```

(Must run after `modules/multistate/transformed_parameters.stan` — it does; transformed parameters are computed before the model block.)

- [ ] **Step 5: Fix the hardcoded `0` at the multistate visit-gating argument**

The last `multistate(...)` argument (currently line 124, the `0` before the closing `);`) must become the real flag. Change:

```stan
          log_cond_surv_03,
          log_cond_surv_32,
          0
        );
```

to:

```stan
          log_cond_surv_03,
          log_cond_surv_32,
          enable_ms_visit_gated_01
        );
```

- [ ] **Step 6: Commit**

```bash
git add stan/tumor/sf-ssm-log-space.stan
git commit -m "feat(laplace): wire joint surrogate into sf-ssm-log-space + fix enable_ms_visit_gated_01"
```

---

## Task 5: Integration compile gate (`stanc` on `sf-ssm-log-space.stan`)

**Files:** none (verification only)

- [ ] **Step 1: Run the syntax check**

Run:
```bash
~/.cmdstan/cmdstan-2.39.0/bin/stanc --include-paths=stan --include-paths=stan/tumor stan/tumor/sf-ssm-log-space.stan
```

Expected: exits 0, no output (clean parse). This validates Tasks 1–4 as a unit.

- [ ] **Step 2: Fix any compile errors, then re-run**

If errors appear, they will cite `file:line` in the module fragments. Common classes: tuple-arg type mismatch (the parameter tuple in `laplace_marginal_tol` must match the functor signature positionally and by type), out-of-scope symbol (an include ordering problem — revisit Task 4 placement), or a `data`-qualifier mismatch (an argument the functor marks `data` but the call passes a parameter-derived value). Fix in the relevant module file, re-commit, re-run until clean.

- [ ] **Step 3: Commit any fixes**

```bash
git add -A stan/
git commit -m "fix(laplace): resolve sf-ssm joint surrogate compile errors"
```

(Skip if Step 1 was clean.)

---

## Task 6: LFO d=2 numeric-equivalence gate

**Files:** none (verification only) — proves the shared-module edit preserved the LFO consumer.

- [ ] **Step 1: Compile-check `sf-ssls-lfo.stan`**

Run:
```bash
~/.cmdstan/cmdstan-2.39.0/bin/stanc --include-paths=stan --include-paths=stan/tumor stan/tumor/sf-ssls-lfo.stan
```

Expected: exits 0. This confirms the d=2 overloads resolve at the LFO call site (which passes the 8-arg `surrogate_ll` and 2-arg `surrogate_K_fn`).

- [ ] **Step 2: Numeric check — LFO log-density unchanged**

The overload forwards to the parametric body with `n_frailty_slots=0`, so the hazard terms are skipped and the result must equal the pre-change SLD-only likelihood. Verify by comparing the marginalized log-density on a fixture before/after.

Run (uses the existing LFO experiment driver if present, else a minimal harness):
```bash
git stash   # revert to pre-change module
# compile + evaluate LFO log-density at a fixed parameter point on a small fixture -> save lp_before
git stash pop
# recompile + evaluate at the SAME point -> lp_after
```

Concretely, if `r/experiments/test_laplace_surrogate_mixed.R` (or the LFO gate driver) exists, run it at fixed seed before and after and diff the reported `lp__` for the surrogate term. Expected: identical to < 1e-6 (only the call path changed, not the math).

- [ ] **Step 3: Record the result**

If lp matches: note PASS in the commit message. If it diverges: the overload forwarding introduced a difference (most likely the degenerate-arg defaults leaking into the SLD term) — debug before proceeding. Do NOT proceed to the production fit on a failed LFO gate.

```bash
git commit --allow-empty -m "test(laplace): LFO d=2 numeric-equivalence gate PASS (lp matches pre-change)"
```

---

## Task 7: Forecast-vs-surrogate feature-scale check

**Files:**
- Create: `r/experiments/test_surrogate_feature_scale.R`

This is the spec's §5-step-3 gate and the one validation the standalone gate deferred: confirm the surrogate's level/velocity features land in the SAME numeric frame as the forecast features, since `tv_coef` is shared and forecast-calibrated.

- [ ] **Step 1: Write the check script**

Create `r/experiments/test_surrogate_feature_scale.R`:

```r
# Forecast-vs-surrogate feature-SCALE check (spec 2026-06-11 §5 step 3).
# For matched (patient, week) pairs, the surrogate level feature
#   standardize_log_burden(b0 + b1*tau + b2*tau^2, log_baseline_burden[p], med, iqr)
# must equal the forecast feature
#   (log_baseline_burden[p] + fmin(log_burden_normalized, 10) - med) / iqr
# to within the uncapped-quadratic vs capped-bi-exponential approximation.
library(tidyverse)
library(cmdstanr)

# Pull the fitted model's standardization constants + a few background patients'
# (b0,b1,b2) surrogate coefficients and their observed log-burden trajectory.
# (Load from the most recent surrogate fit store; parameterize the path.)
store <- file.path("/mnt/data/analysis-results", Sys.getenv("DOMINO_STARTING_USERNAME"),
                   "publication", Sys.getenv("TAR_RUN", "laplace"), "_targets")
stopifnot(dir.exists(store))

# median/iqr_log_burden_obs and median/iqr_velocity_obs are in the stan data.
# Reconstruct g(tau) for sampled background patients from posterior (b0,b1,b2)
# means, evaluate both features on the patient's observed visit weeks, and
# compare. Report max abs diff in level-feature and velocity-feature units.
standardize_log_burden <- function(g_norm, log_base, med, iqr) {
  (log_base + pmin(g_norm, 10)) - med
} # divided by iqr below

# ... (load constants + per-patient coefficients from the store; compute both
#     feature vectors; tabulate diffs) ...

cat(sprintf("Level feature: max |surrogate - forecast| = %.4f (std units)\n", max_lvl_diff))
cat(sprintf("Velocity feature: max |surrogate - forecast| = %.4f (std units)\n", max_vel_diff))
stopifnot(max_lvl_diff < 0.5)  # threshold: well within one IQR
```

(The store-loading details depend on which surrogate fit is available; parameterize `TAR_RUN` and the patient sample. The check is meaningful only against a fit where the surrogate ran — so this task may run AFTER an initial short production fit, or against the gate's synthetic trajectories as a dry run.)

- [ ] **Step 2: Run the check (dry-run on synthetic if no fit yet)**

Run:
```bash
Rscript r/experiments/test_surrogate_feature_scale.R
```

Expected: max level-feature diff < 0.5 std units (the only material difference is the `fmin(·,10)` cap, which binds only in the extreme growth tail). A larger diff means the frame match is wrong — revisit the §3 functor.

- [ ] **Step 3: Commit**

```bash
git add r/experiments/test_surrogate_feature_scale.R
git commit -m "test(laplace): forecast-vs-surrogate feature-scale equivalence check"
```

---

## Task 8: Update model specification docs

**Files:**
- Modify: `quarto/sclc/website/documentation/multistate-specification.qmd` (or the publication equivalent)

Per the r-guidelines "model specification is the blueprint" rule: when Stan changes, the spec page must reflect it.

- [ ] **Step 1: Add a subsection describing the background-patient marginalization**

Document: the forecast/background split, the joint SLD+MS surrogate, the d=2+n_frailty latent layout, the block-diagonal K, and that background patients' hazard coupling uses the same (level, velocity) feature frame as forecast patients. Reference the spec doc.

- [ ] **Step 2: Commit**

```bash
git add quarto/sclc/website/documentation/multistate-specification.qmd
git commit -m "docs(spec): document background-patient Laplace marginalization in multistate spec"
```

---

## Final verification

- [ ] All 8 tasks committed.
- [ ] `stanc` clean on both `sf-ssm-log-space.stan` (Task 5) and `sf-ssls-lfo.stan` (Task 6).
- [ ] LFO numeric-equivalence gate PASS (Task 6).
- [ ] Feature-scale check PASS (Task 7).
- [ ] After all gates: hand to an adversarial-review workflow on the IMPLEMENTATION (spec §5 step 4) before the production fit.
