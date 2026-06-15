# Joint SLD+MS Laplace Surrogate — Production Wiring Design (2026-06-11)

Status: **DESIGN — awaiting human review.**

This design specifies how to wire the **validated joint SLD+multistate Laplace
surrogate** (with correlated patient-level frailty) into the publication model
`stan/tumor/sf-ssm-log-space.stan`. It is the implementation-side companion to the
reviewed **index contract** (`2026-06-09-joint-surrogate-index-contract.md`), which
remains the source of truth for index spaces, invariants, and the worked numeric
example. Where this design and the contract overlap, the contract governs; this
document adds the concrete file-level plan and resolves the implementation choices
the contract left open.

Prerequisites already satisfied (do not re-litigate):
- d=4 joint+frailty marginalization validated standalone
  (`laplace_joint_frailty_test.stan`, gate **conditional PASS** 2026-06-11, commit
  `f64d4105`). Correctness signal: min SD ratio 0.96 (Laplace not overconfident),
  0 divergences both modes; the borderline 5.59 diff/MCSE is a noisy-HMC-reference
  artifact, not a Laplace error.
- SLD-only (d=2) surrogate integrated as `stan/modules/laplace_surrogate/`, live in
  `sf-ssls-lfo.stan`, flag-gated behind `enable_background_surrogate`.
- Forecast/background split present in `sf-ssm-log-space.stan`: the tumor SLD loop
  already gates over `forecast_patient_idx` (model block :91-92) and the
  `multistate(...)` call is already sliced by `forecast_patient_idx` (:107-125),
  via commits `98ac9041` / `1952630a`. **Two of the contract's three §1 BLOCKER
  edits are therefore already done.**
- (level, velocity) coupling basis landed model-wide (`enable_ms_velocity_basis`,
  `n_time_varying_covar = 2`); publication config confirmed.
- cmdstan 2.39 toolchain pinned.

---

## 0. Decision log (this session)

1. **Extend the module in place** (d-parametric), NOT a parallel functor. Rationale:
   the publication model compiles to identical C++ either way (same d=4 joint
   functor → same `laplace_marginal_tol` → same block-diagonal Hessian), so there is
   **no computational difference**; the source-org choice only affects whether the
   LFO model shares the bytes. Extend-in-place gives a single source-of-truth `d`
   (contract §3) and a smaller review surface.
2. **One parametric functor + a d=2 overload for LFO** (review `wf_922c4cd6`
   BLOCKER 2/3). The shared `surrogate_ll` / `surrogate_K_fn` / `likelihood.stan`
   call site is `#include`d by both models and IS rewritten; a Stan overload with the
   original d=2 signature forwards to the parametric body with degenerate args so
   LFO's numeric result is unchanged. "Byte-identical generated code" is NOT claimed
   (runtime guards transpile differently); **numerically-identical d=2 reduction** is,
   verified by a numeric gate (§5 step 2).
3. **Both 0→1 and 0→3 couple to burden** (confirmed: `time_varying_coef_03` reads the
   same `ms_time_varying_covar_01` level/velocity features as 0→1,
   `transformed_parameters.stan:780`). The functor carries coupling on both.
4. **Surrogate feature MUST match the forecast feature's coordinate frame** (review
   `wf_922c4cd6` BLOCKER 1). `tv_coef` is shared model-wide and forecast-calibrated,
   so the surrogate level feature is built with the SAME `standardize_log_burden(...)`
   helper, including `b0`, the per-patient `log_baseline_burden[p]` shift, and the
   anchored time `τ = w − baseline_week[p]`. See §3.
5. **No fresh synthetic marginalization gate** — the d=4 marginalization is already
   validated. The production-specific risk is the coordinate-frame match (decision 4),
   validated empirically at wiring time (§5 step 3).

---

## 1. File scope & include wiring

The change is confined to two files plus config.

### A. `stan/modules/laplace_surrogate/` — extend in place
- **`surrogate.stanfunctions`**: generalize the joint functor and `surrogate_K_fn`
  to `d = 2 + n_frailty_slots`, adding the MS-hazard coupling + frailty terms ported
  from the validated gate functor `joint_frailty_ll`
  (`laplace_joint_frailty_test.stan:51-121`). `surrogate_bridge` (Σ_β GH-3
  pushforward) is unchanged.
  - **Functor: ONE parametric implementation + a thin LFO overload** (decision
    2026-06-11; review `wf_922c4cd6` BLOCKER 2/3). The full `surrogate_ll` takes the
    extended argument set (static baselines, Σ_u-derived inputs are in K, `tv_coef`,
    `log_baseline_burden`, standardization constants, `d`, MS event/visit arrays).
    Provide a **Stan function overload** `surrogate_ll(theta, beta_pop, measure_sd,
    log_lod_per_visit, n_bg, bg_obs, bg_pos, bg_time)` — the EXISTING d=2 signature —
    that simply forwards to the parametric body with degenerate
    (`[0]`-sized / zero-length) hazard inputs and `d=2`. LFO calls the overload
    unchanged; `sf-ssm` calls the full signature. Stan resolves overloads by arity +
    type, so the LFO `#include`d call site needs **no edit** and its numeric result is
    the d=2 SLD-only likelihood. (`surrogate_K_fn` similarly: keep a 2-arg
    `(Σ_β, n_bg)` overload forwarding to the parametric `(Σ_β, Σ_u, n_patients, d)`
    with a `0×0` Σ_u and `d=2`.)
- **`transformed_data.stan`**: compute `surrogate_d` / `n_frailty_slots` /
  `frailty_slots` / `coupled_tv_slots` from the MS config; build the bg-compact MS
  arrays (event/censor weeks per transition, 0→1 visit-gating mask, per-bg
  `baseline_week`) alongside the existing SLD views; resize `surrogate_theta_0` and
  `surrogate_hessian_block_size` to `surrogate_d`; add the §5 guards.
- **`likelihood.stan`**: assemble the param-derived static baselines
  (`base01_static`, `base03_static`), `Σ_u`, and `K_i = blockdiag(Σ_β, Σ_u)`; gather
  per-bg `log_baseline_burden`; call `laplace_marginal_tol` with the extended
  functor/data/param tuples. **This single call site is `#include`d by BOTH models**
  — so the LFO branch (`n_frailty_slots==0`) must select the overload / degenerate
  assembly. The cleanest form: an `if (n_frailty_slots > 0) { full call } else { d=2
  overload call }` so LFO's transpiled path is the original d=2 call.
- **`data.stan` / `flags.stan`**: unchanged (config flows through
  `enable_background_surrogate` + the existing MS flags).

### B. `stan/tumor/sf-ssm-log-space.stan` — 4 includes + 1 fix
- Add to the `functions` block: `#include "modules/laplace_surrogate/surrogate.stanfunctions"`.
- Add to the `data` block: `modules/laplace_surrogate/flags.stan`,
  `modules/laplace_surrogate/data.stan`.
- Add to the `transformed data` block: `modules/laplace_surrogate/transformed_data.stan`
  (AFTER `modules/multistate/transformed_data.stan` — needs `MS_SLOT_*`,
  `ms_level_intercept_mode`, `background_patient_idx`).
- Add to the `model` block: `modules/laplace_surrogate/likelihood.stan` (AFTER
  `modules/multistate/transformed_parameters.stan` has run — needs
  `log_level_lambda_01/03_residual`, `ms_corr_u` construction inputs).
- **Fix the hardcoded `0`** at the last `multistate(...)` argument (model block
  :124) → the real `enable_ms_visit_gated_01` flag. This is a **pre-existing
  forecast-path config bug, independent of the surrogate** (review `wf_922c4cd6`
  MINOR 4): with the hardcoded `0`, FORECAST patients are routed through the
  continuous/interval-censored 0→1 path instead of the visit-gated path that
  `publication_targets.R:319` (`enable_ms_visit_gated_01 = TRUE`) requests. The
  dense path covers forecast patients and the surrogate covers background patients
  — DISJOINT sets — so this fix is NOT about aligning the dense path with the
  surrogate; it is about the forecast likelihood honoring its own config flag. It
  lands alongside the surrogate work because both touch the same model block. (The
  line-297 GQ comment "no PSA covariate" is stale pioneer-PSA framing, irrelevant
  to the SLD likelihood — leave the GQ `burden_enable_ms_visit_gated_01 = 0` alone;
  it governs simulation, not the likelihood.)

Mirror the include placement used by `sf-ssls-lfo.stan` (lines 13/37/38/61/120).

### LFO safety
`surrogate.stanfunctions` and `likelihood.stan` are shared with `sf-ssls-lfo.stan`,
which `#include`s them verbatim — Stan has no optional args, so the SAME call site
runs in both models. The plan therefore does NOT claim "LFO is not modified": the
shared call site IS rewritten. What is guaranteed is a **numerically-identical d=2
reduction** for LFO, achieved by (a) the d=2 `surrogate_ll` / `surrogate_K_fn`
overloads that forward to the parametric body with degenerate hazard inputs, and
(b) the `n_frailty_slots==0` branch in `likelihood.stan` selecting the original d=2
call. Runtime guards transpile to C++ that differs from today's — so the guarantee
is byte-identical *numeric behavior*, not byte-identical *generated code*. Enforced
by §5 step 2: recompile `sf-ssls-lfo.stan` AND a numeric d=2-equivalence check (run
the LFO surrogate path on a fixture before/after; log-density must match to
tolerance).

---

## 2. Latent layout & `d` detection

Per-background-patient latent block (length `d`), stacked burden-first then frailty
in ascending MS-slot order:

```
θ_i = [ b1_i, b2_i,            ← burden RE (always): SLD μ + both hazards' coupling
        u_{s1,i}, …, u_{sF,i} ]  ← one frailty latent per frailty slot, ascending slot order
```

```
d = 2 + n_frailty_slots
```

Detection (module `transformed_data.stan`, data-int shaped):

```
patient_lv     = n_levels                                  // patient level = last hierarchy column
coupled_tv_slots = { slots whose hazard reads the burden coupling }   // {01,03} publication
frailty_slots  = { slot ∈ coupled : ms_level_intercept_mode[slot, patient_lv] == 2 }
n_frailty_slots = |frailty_slots|
```

- **Publication**: `ms_level_intercept_mode[MS_SLOT_01, patient_lv] == 2` and
  `[MS_SLOT_03, patient_lv] == 2` (`publication_targets.R:408-409`) ⇒
  `frailty_slots = {01,03}`, `n_frailty_slots = 2`, `d = 4`.
- **LFO**: neither set ⇒ `frailty_slots = ∅`, `d = 2` — the **numerically-identical
  d=2 reduction** (via the d=2 overload; §1.A). Functor collapses to the SLD term only.
- `coupled_tv_slots` (drives functor body) and `frailty_slots` (drives latent dim +
  K) are detected and asserted **separately** (contract §2): they coincide as
  {01,03} for the publication but must not be assumed equal in general.

`b0 = beta_pop[1]` stays **pinned** (intercept fixed at the population value, not a
latent) — consistent with both the validated gate and today's d=2 `surrogate_ll`.

**Single source of truth for `d`**: the literal `int d = 2` currently hardcoded in
`surrogate_ll` (`:88`), `surrogate_K_fn` (`:109`), and the `transformed_data` setup
is replaced by one `surrogate_d`, passed to the functor as a `data int` (the gate
already does this), with `assert surrogate_hessian_block_size == surrogate_d`.

---

## 3. Extended `surrogate_ll` functor body

Per background patient `i` (bg-compact). Unpack:

```
b0  = beta_pop[1]                          // pinned intercept
b1  = beta_pop[2] + θ[(i-1)d + 1]
b2  = beta_pop[3] + θ[(i-1)d + 2]
u01 = θ[(i-1)d + 3]   (if 01 ∈ frailty_slots)
u03 = θ[(i-1)d + 4]   (if 03 ∈ frailty_slots)
```

**Term 1 — SLD (unchanged from d=2 path):**
```
for each visit v of patient i:
    t = bg_time[v] − 1.0                       // existing anchored offset (baseline → t=0 → g(0)=b0)
    μ = b0 + b1·t + b2·t²
    obs > 0 :  lp += normal_lpdf(log(obs) | μ, measure_sd)
    obs ≤ 0 :  lp += normal_lcdf(log_lod_per_visit[v] | μ, measure_sd)   // LOD-censored
```

**Term 2 — 0→1 hazard (visit-gated; covariate = surrogate burden):**
```
te01 = ms_censored_01[i] ? n_wk : ms_event_wk_01[i]
cum = 0
for w in 1:te01:                                          // w = calendar week (matches base_T_static column)
    if visit_wk_01[i,w] == 1:                             // only at observed visit weeks
        τ      = w − baseline_week[p]                      // anchored time (baseline → 0); SAME frame g was fit in
        g_norm = b0 + b1·τ + b2·τ²                         // NORMALIZED log-burden (b0 PINNED, included)
        lvl_feat = standardize_log_burden(g_norm, log_baseline_burden[p],
                                          median_log_burden_obs, iqr_log_burden_obs)
        vel_feat = standardize_velocity(b1 + 2·b2·τ, median_velocity_obs, iqr_velocity_obs)
        loghaz = base01_static[i,w] + u01
               + tv_coef_01[1]·lvl_feat + tv_coef_01[2]·vel_feat
        cum += exp(loghaz)
        if (!ms_censored_01[i] && w == te01) lp += loghaz
lp += −cum
```

**Term 3 — 0→3 hazard (continuous; every week):** identical shape, NO `visit_wk`
gate, runs `w in 1:te03`, uses `base03_static[i,w]`, `u03`, `tv_coef_03[1..2]`, the
SAME `lvl_feat`/`vel_feat` (0→3 reads `ms_time_varying_covar_01`,
`transformed_parameters.stan:780`).

### Coordinate-frame match (contract §4 pin — RESOLVED by matching the forecast feature exactly)
**This is the load-bearing correctness requirement** (adversarial review
`wf_922c4cd6` BLOCKER 1). `tv_coef_01/03` are shared model-wide and calibrated on
the FORECAST patients' feature, so the surrogate's feature MUST live in the
identical coordinate frame, built by the SAME helper:

- **Forecast feature** (`_ms_burden_tv_covar.stan:64-65`):
  `log_burden_absolute = log_baseline_burden[p] + fmin(log_burden_normalized, 10)`;
  `feature = (log_burden_absolute − median_log_burden_obs) / iqr_log_burden_obs`
  = `standardize_log_burden(log_burden_normalized, log_baseline_burden[p],
  median_log_burden_obs, iqr_log_burden_obs)` (`_burden.stanfunctions:34-39`).
- **Surrogate feature** MUST call the **same `standardize_log_burden(...)` helper**
  with `g_norm = b0 + b1·τ + b2·τ²` as the normalized-burden argument and the
  patient's `log_baseline_burden[p]`. Three corrections vs. the naive
  `lvl = b1·w + b2·w²` the review caught:
  1. **Include `b0`** (the pinned intercept) — `g` is `b0+b1·τ+b2·τ²`, not `b1·τ+b2·τ²`.
  2. **Add `log_baseline_burden[p]`** — lifts normalized→absolute; a per-patient
     additive shift a global median/iqr cannot absorb. The helper does this.
  3. **Anchor time correctly**: the surrogate quadratic is fit in baseline-anchored
     time (`t = bg_time[v]−1`, baseline→0, `surrogate.stanfunctions:97`). The hazard
     loops calendar week `w` (to index `base_T_static[i,w]`), so evaluate `g` at
     `τ = w − baseline_week[p]` to land in the same anchored frame. (The forecast
     trajectory column `w` is the burden at calendar week `w` via
     `states_start_col = 2 − first_visit`; matching `base_T_static`'s calendar-week
     column keeps all three — baseline, hazard week, residual column — consistent.)
- **The `fmin(·,10)` cap is NOT applied to the surrogate** — the uncapped analytic
  quadratic is the deliberate log-concave surrogate (contract O1). Applying `fmin`
  would break log-concavity. The standardize helper's internal `fmin` operates on
  the normalized argument; passing the uncapped quadratic means the cap only binds
  in the extreme tail, an accepted approximation.
- **New functor arguments** (add to contract §4 table): `log_baseline_burden[p]`
  (per-bg-patient, param-derived → parameter tuple), `median/iqr_log_burden_obs`,
  `median/iqr_velocity_obs`, `baseline_week[p]` (data) — all already exist in the
  model (`_tumor_observed_covar` transformed data + `_full_model_transformed_data`).

### Properties preserved
- **Log-concavity.** `g_norm` and `b1+2·b2·τ` are affine in `(b1,b2)`;
  `standardize_*` are affine maps (uncapped in the operative range);
  `base_T_static[i,w]` and `u` are additive. So `loghaz` is affine in `(b1,b2,u)`
  and `−Σexp(loghaz)` stays concave in all `d` latents — the property the gate
  verified. Solver 1.
- **Numerically-identical d=2 reduction.** When `n_frailty_slots = 0` AND no coupled
  slots (LFO), Terms 2-3 contribute nothing and the `u` latents do not exist; the
  functor reduces to Term 1, **numerically** identical to today's d=2 SLD-only
  result (NOT byte-identical generated code — runtime guards transpile to different
  C++; see §1 LFO safety). Guard the hazard terms behind coupled-slot presence.

---

## 4. Model-block assembly (`likelihood.stan`)

Computed once per leapfrog step, OUTSIDE the Laplace integral (the Newton solver
re-evaluates the functor many times per step; these must not be inside it).

### (a) Static baseline vectors
`base01_static`, `base03_static` sized `[n_background_patients, max_all_t]`:

```
// TI-cov linear predictors in QR space, evaluated on the BACKGROUND rows
// (mirrors the forecast linpred at transformed_parameters.stan:218,375 but with
//  background_patient_idx; population-level term only — publication sets
//  enable_ms_level_cov = c(trial=FALSE, patient=FALSE), so NO level slopes):
linpred_bg_01 = Q_covar_design_matrix[background_patient_idx, :] * time_invariant_coef_qr_01   // [n_bg]
linpred_bg_03 = Q_covar_design_matrix[background_patient_idx, :] * time_invariant_coef_qr_03   // [n_bg]

for i in 1:n_background_patients:
    p = background_patient_idx[i]
    base01_static[i] = log_pop_lambda_01                                   // population row (row_vector over weeks)
    base03_static[i] = log_pop_lambda_03
    for lv in enabled baseline levels with lv < n_levels:                  // trial-level baseline residual
        base01_static[i] += log_level_lambda_01_residual[ patient_ms_baseline_flat_idx_slot[MS_SLOT_01, p, lv] ]
        base03_static[i] += log_level_lambda_03_residual[ patient_ms_baseline_flat_idx_slot[MS_SLOT_03, p, lv] ]
    base01_static[i] += linpred_bg_01[i]                                   // TI-cov offset, broadcast over weeks
    base03_static[i] += linpred_bg_03[i]
```

Two mirrored sources: the **baseline residual gather** copies the forecast path at
`transformed_parameters.stan:164`; the **TI-cov linear predictor** copies the
forecast path at `:218` / `:375`, both keyed by `background_patient_idx` instead of
`forecast_patient_idx`, accumulated into a plain matrix instead of `log_cond_surv`.
TI-cov uses the QR-space coefficients `time_invariant_coef_qr_01/03` against
`Q_covar_design_matrix` (NOT a raw `coef · covar` product). Level-slope TI terms are
absent under the publication config (`enable_ms_level_cov` both FALSE); if a future
config enables them, the gather must extend to the `ms_scaled_level_slope_*` path
(`:223-254`) — assert or handle.

**Background-safety** (contract §1 BLOCKER resolution): for `lv < n_levels`,
`patient_ms_baseline_flat_idx_slot[slot, p, lv]` returns the **trial-level group
row** (`get_global_group_idx`), correct for every patient. The `→1` fallback
(`transformed_data.stan:391`) fires only at the **patient level**
(`lv == n_levels`), which is never gathered — the patient-level intercept is the
frailty `u`, which lives in θ. No forecast-sized array is touched. Guard 6 protects
the invariant.

Passed into the functor through the **parameter tuple** (parameter-derived ⇒ cannot
be `data`-qualified like the SLD views).

### (b) Σ_u frailty covariance (n_frailty_slots × n_frailty_slots)
Same construction the model already uses for `ms_corr_u`
(`transformed_parameters.stan:37`) — no new hyperparameters:

```
σ[m] = log_lambda_gp_<slot_m>_level_intercept_sd[patient_lv]   for slot_m ∈ frailty_slots
     // publication: σ = [ log_lambda_gp_01_level_intercept_sd[patient_lv],
     //                     log_lambda_gp_03_level_intercept_sd[patient_lv] ]
L    = L_ms_intercept_corr[b]                                  // patient-level corr block b
                                                               //   (ms_corr_block_level[b] == patient_lv)
Σ_u  = diag(σ) · (L Lᵀ) · diag(σ) + jitter·I
```

### (c) Σ_β burden covariance
Unchanged: `surrogate_bridge(...)` GH-3 pushforward, as the d=2 path computes it.

### (d) K_i = blockdiag(Σ_β, Σ_u)
`surrogate_K_fn` extended to tile the `d×d` block per patient (today tiles 2×2).
Zero burden↔frailty cross-block — **prior-exact** (contract O2): burden REs use
independent `std_normal` NCP priors, frailty is a separate LKJ+Gaussian block; the
block-diagonal K introduces zero prior approximation error.

---

## 5. Guards & validation

### Guards (module `transformed_data.stan`, fire only when `enable_background_surrogate == 1`)

| # | Assert | Prevents |
|---|---|---|
| 1 | `enable_ms_level_gp[slot, patient_lv] == 0` for slot ∈ coupled ∪ frailty | time-varying patient residual (not a marginalizable scalar) |
| 2 | `ms_level_intercept_mode[slot, patient_lv] ∈ {0, 2}` for coupled slots | FE/other patient intercept the surrogate can't integrate out |
| 3 | `n_frailty_slots == count(mode==2 on coupled slots)` | latent dim silently drifting from the model's frailty |
| 4 | `ms_corr_block_member_slot[b,·] == frailty_slots` (same order, block `b` with `ms_corr_block_level[b]==patient_lv`) | Σ_u rows misaligned with θ frailty entries |
| 5 | retained: tr SD sub-hier (`:28-33`), `surrogate_anchor_times[1]==0` (`:36-38`) | (existing) |
| 6 | `0 < forecast_split_level < n_levels` | baseline gather falling into the `→1` patient-level fallback |

Guards 3 & 4 are load-bearing: the latent vector, Σ_u, and the static gather are
three independent representations of the same frailty structure; any disagreement
would silently produce a wrong posterior. The asserts make drift a `fatal_error` at
data-prep time.

### Validation strategy (implementation-plan gates, in order)
1. **`stanc` syntax check** of `sf-ssm-log-space.stan` (new joint path).
2. **LFO d=2 equivalence** — `stanc` recompile of `sf-ssls-lfo.stan` (the shared
   module still compiles for the d=2 consumer) **AND a numeric check**: run the LFO
   surrogate path on a fixture before/after the change; the marginalized log-density
   must match to tolerance. The d=2 overload + `n_frailty_slots==0` branch guarantee
   numeric (not generated-code) identity — so a compile pass alone is insufficient.
3. **Forecast-vs-surrogate feature-scale check** (contract O1 / §9 step 2 leftover):
   confirm the surrogate level feature, built via `standardize_log_burden(b0+b1·τ+b2·τ²,
   log_baseline_burden[p], median_log_burden_obs, iqr_log_burden_obs)` (§3), is
   numerically in the SAME frame as the forecast feature
   (`_ms_burden_tv_covar.stan:64-65`) for matched (patient, week) pairs — since
   `tv_coef` is shared and forecast-calibrated. This is the §3 BLOCKER-1 fix; verify
   it empirically, not just by inspection. The only material frame difference left is
   the surrogate's uncapped quadratic vs. the forecast's `fmin(·,10)` cap (an accepted
   log-concave approximation, binds only in the extreme tail).
4. **Adversarial-review workflow** on the finished implementation (multi-agent;
   mirrors the contract's `wf_9c08ae27`): index correctness, log-concavity
   preservation, double-counting (surrogate replaces — not supplements — the dense
   likelihood for background patients), OOB safety.
5. **Production fit** — only after 1-4 pass.

No fresh synthetic gate: the d=4 marginalization is validated (commit `f64d4105`).
The production-specific risk is step 3, not the marginalization itself.

---

## 6. Out of scope
- The (level, velocity) coupling-basis core model change (landed; its own spec).
- Any change to shared multistate array *shapes* (`log_cond_surv_*` stay
  `[n_forecast_patients, max_all_t]`; 8 models `#include` them).
- The LFO model's behavior (must stay **numerically identical** in its d=2 path; the
  shared call site IS rewritten — see §1 LFO safety — but the d=2 overload +
  `n_frailty_slots==0` branch preserve the result, verified by the §5 step-2 numeric
  gate).
- Re-validating the marginalization math (done).
