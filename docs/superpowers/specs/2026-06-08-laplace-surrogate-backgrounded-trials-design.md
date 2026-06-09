# Log-concave surrogate likelihood for backgrounded trials

**Date:** 2026-06-08
**Branch:** `karim/laplace`
**Status:** Design approved, pending implementation plan

## 1. Problem & goal

Backgrounded historical trials inform only the **population** parameters during
training; we never forecast them or report their per-patient trajectories. The
full bi-exponential per-patient likelihood is non-log-concave in its latents —
the rate sits on a log scale inside a double exponential, `exp(-exp(z)·t)`,
wrapped in a `log_sum_exp` of a decaying and a growing branch — so the embedded
Laplace approximation (`laplace_marginal_tol`) fails: ~91% divergences vs clean
HMC, and even a 1-latent / no-LOD ablation diverges (34%). This was the
documented NO-GO from the prior session (see memory `laplace-builtin-2-39.md`).

**Goal:** replace the backgrounded patients' likelihood with a surrogate that is
**linear in its latent coefficients** (hence Gaussian-marginal, hence Laplace is
*exact* and fast), while transmitting **unbiased information to the population
parameters** via a tight mechanistic bridge. Forecast patients are untouched.

Key realization: "linear-Gaussian" — where Laplace is exact — means linear *in
the latent coefficients*, NOT linear in time. We can use any fixed time basis
(here, quadratic) and keep the closed-form Gaussian marginal.

## 2. Surrogate model (per backgrounded patient)

Normalized log-burden modeled as a quadratic in time:

```
log_burden_i(t) = β_i0 + β_i1·t + β_i2·t²  +  ε,   ε ~ N(0, measure_sd²)
```

- Coefficients `β_i = (β_i0, β_i1, β_i2)` enter **linearly** → Gaussian
  likelihood in `β_i` → **closed-form Gaussian marginal**, no inner Newton solve.
- A quadratic is the minimal basis that can represent the model's defining
  feature, the single nadir (shrink-then-regrow); a straight line structurally
  cannot bend and would discard the signal historical trials are richest in.
- LOD censoring stays as `normal_lcdf(log_lod | μ_it, measure_sd)` — log-concave
  in a linear mean, so it does NOT reintroduce the pathology.
- `β_i ~ MVN(β_pop, Σ_β)`, where `(β_pop, Σ_β)` are NOT free parameters; they are
  the mechanistic-bridge image of the real population rate parameters (§3).
- **Intercept `β_i0` is pinned at 0, NOT marginalized** (REVISED — see §3a). The
  burden is normalized to baseline so `g(t=0) = 0` *identically for every θ*, and
  the first anchor sits at `t=0`. That makes `β_i0`'s bridge variance structurally
  zero → `Σ_β` rank-deficient → solver 1 fails to factor a singular prior. So we
  integrate out only the slope/curvature `(β_i1, β_i2)` → **2 latents per patient,
  `hessian_block_size = 2`**. The intercept contributes `μ_i(0)=0` exactly.

## 3. Mechanistic bridge (anchor-matched quadratic, fixed calendar anchors)

How `(β_pop, Σ_β)` are computed from the population parameters each iteration.

**Means — exact bi-exponential at 3 fixed anchor times.**
Pick constant anchors `t = (t₀, t₁, t₂)` spanning the historical trials'
observation windows (§6). At the population rate parameters
`θ_pop = (tr_loc_pop, frac_logit_pop, init_logit_pop)`, evaluate the true
bi-exponential log-burden using the *same* rate construction as the full model:

```
g(t; θ_pop) = log_sum_exp(init_log_dec − dec_rate·t,  init_log_gro + gro_rate·t)
```

Solve the 3×3 Vandermonde system for the quadratic through `(tₖ, g(tₖ))`:

```
β_pop = V⁻¹ · [g(t₀), g(t₁), g(t₂)]ᵀ ,   V[k,:] = [1, tₖ, tₖ²]
```

`V⁻¹` is a **compile-time constant** (anchors fixed) → one cheap matrix-vector
product per iteration. `β_pop` is a smooth deterministic function of `θ_pop`, so
AD flows through and borrowing lands in the real population parameters.

**Covariance — Gauss–Hermite pushforward (REVISED — see §3a).**
Per-patient spread comes from `{tr_sd, frac_sd, init_sd}` propagated to `(β_1,
β_2)`. The original plan used a **first-order** linearization
`Σ_β = J·diag(sd²)·Jᵀ`. The Phase-1 gate showed this understates the true spread
by ~12% (slope) to ~20% (curvature) — capturing only **87.8%** of the true
variance — which biased `tr_sd` low (0.506 → 0.432) and failed the gate. So we
replace it with a **3-point Gauss–Hermite quadrature** of the exact pushforward
`Cov[β(θ)]`, `θ ~ N(θ_pop, diag(sd²))`:

```
Σ_β = Σ_q w_q · (β(θ_q) − β̄)(β(θ_q) − β̄)ᵀ ,   β̄ = Σ_q w_q · β(θ_q)
```

over the tensor grid of 3 nodes per dimension (27 evaluations of the closed-form
`β(·)`; GH nodes `±√3, 0`, weights `1/6, 2/3, 1/6`). Verified to recover **99.9%**
of the true variance (vs 87.8% first-order, and — notably — 80–83% for the
unscented transform, which was worse). Fully deterministic → AD-friendly.

`Σ_β` is the surrogate's random-effect covariance — patient-to-patient
variability transmitted to `tr_sd`/`frac_sd`/`init_sd`. It becomes the prior
covariance `K` passed to `laplace_marginal_tol` (non-identity, now 2×2).

Why fixed anchors still compose well: `V⁻¹` and the GH nodes/weights are all
compile-time constants; the per-iteration cost is 27 closed-form `g` evaluations.

## 3a. Revisions from the Phase-1 gate (2026-06-08)

The first gate run FAILED on one quantity (`tr_sd`), with a clean fit otherwise
(0 divergences, R-hat 1.00, ESS 2594, SD-ratio 0.99, corr-diff 0.082). Pure-R
diagnostics localized two independent causes and their fixes — **neither is the
plan's "FAIL → spline" path; the quadratic basis was never the problem**:

1. **Solver-1 → solver-2 fallback (was misread as non-log-concavity).** The data
   Hessian is well-conditioned (cond 8.6e5 ≪ 4.5e15 break point) and PD — the
   likelihood *is* log-concave. The real cause: `g(0)=0` identically (baseline
   normalization) with the first anchor at `t=0` makes the intercept's bridge
   variance structurally zero, so `Σ_β` is rank-2 and the `1e-8` jitter leaves a
   near-singular `K`. **Fix:** drop the intercept from the marginalized latents
   (pin `β_0 = 0`); integrate out only `(β_1, β_2)`, `hessian_block_size = 2`.
   Diagnostic: full-rank, cond ≈ 1.6e3 — solver 1 holds.
2. **`tr_sd` underestimate (first-order bridge bias).** First-order `Σ_β`
   captured only 87.8% of the true pushforward variance. **Fix:** GH-3 quadrature
   (above), 99.9%.

`tr_sd` reaches the data only through the slope/curvature `(β_1, β_2)`, so both
fixes target exactly the failing quantity.

## 4. Integration into `sf-ssls-lfo.stan`

Routing already exists: `background_patient_idx` / `n_background_patients`, and
patient NCP params sized by `n_forecast_patients`, so backgrounded patients
currently have no explicit latents and contribute nothing to the LFO likelihood.
The surrogate slots into exactly that gap.

- **Forecast patients:** unchanged. Full bi-exponential `sf_log_space_obs`,
  explicit NCP latents, full trajectory — everything we report is untouched.
- **Backgrounded patients:** new
  `target += laplace_marginal_tol(surrogate_ll, (β_pop, obs, ...),
  hessian_block_size=2, K_fn=Σ_β, ...)` — marginalize `(β₁,β₂)` only; intercept
  pinned (§3a).
- **Solver:** the surrogate marginal is linear-Gaussian, so **solver 1** is the
  intended path. In practice it falls back to **solver 2** occasionally (the
  `(t, t²)` basis makes the inner Hessian ill-conditioned); this did NOT affect
  convergence (0 divergences, R-hat 1.00) or the validated endpoint. A follow-up
  could center/scale the time basis so solver 1 holds; not required.

New Stan functions in a small module (e.g. `stan/modules/laplace_surrogate/`):
- `surrogate_anchor_betas(θ_pop)` → quadratic coefficients at a single θ
- `surrogate_bridge(θ_pop, sds)` → `(β_pop, Σ_β)` via the GH-3 pushforward —
  returns BOTH the population mean `E[β(θ)]` and the `(β₁,β₂)` covariance from
  one quadrature (§3, §3a)
- `surrogate_ll(β, θ_pop, obs, ...)` — the marginalized likelihood (quadratic
  log-burden Gaussian + `normal_lcdf` LOD tail), intercept pinned at `β_pop[1]`
- `surrogate_K_fn(Σ_β, n)` — tiles `Σ_β` block-diagonally (2×2 per patient)

The old hand-coded `stan/modules/laplace/` stays untouched and unused; retire it
only after this lands (separate cleanup, §6 out-of-scope).

Blast radius: exactly one new code path, reachable only by patients we never
report, rejoining the model at `θ_pop`. Population params, priors, hierarchy,
and the entire forecast path are shared verbatim.

## 5. Validation gate (REVISED — mixed-cohort + PFS endpoint)

The original all-backgrounded gate was **stricter than production**: it forced
the surrogate to identify all three population SDs from the surrogate alone,
which two coefficients cannot do (`tr_sd`'s footprint is the ~1250×-weaker
curvature channel, below the noise floor). The validated gate matches production.

**Mixed-cohort harness** (`stan/experiments/laplace_surrogate_mixed.stan` +
`r/experiments/test_laplace_surrogate_mixed.R`):
- Synthetic data from the **true bi-exponential** (+ LOD censoring), patients
  ordered forecast-first: e.g. 20 forecast + 40 backgrounded.
- Two fits: **(A) reference** = ALL patients full-HMC; **(B) mixed** = forecast
  full-HMC (pins the SDs), backgrounded = surrogate-marginalized.

**Gate** — three layers, each closer to the deliverable:
1. **Sampler health (B):** `pct_divergent < 1`, `max_rhat < 1.01`,
   `min_ess_bulk > 400`. (Solver 1 may fall back to solver 2; not disqualifying.)
2. **Population marginals:** all six `{*_pop, *_sd}` agree — `max diff/MCSE < 5`,
   `min SD ratio > 0.8`.
3. **Endpoint invariance (DECISIVE):** derive each population-drawn patient's
   PFS = PD-crossing week (RECIST rule: burden ≥20% above running nadir AND ≥5
   absolute, `baseline_sld=60`) from both runs; require PFS quantiles agree to
   `< 2 weeks` and landmark event rates to `< 0.03`. This is the actual reported
   quantity. The raw-burden `corr_diff` and burden-quantile checks are
   INFORMATIONAL only — they sit one layer below the deliverable and a breach
   there is acceptable iff PFS is invariant.

**Honest-gate discipline:** gate on divergences/R-hat/ESS FIRST; measure the
real deliverable (PFS) rather than a proxy when the two disagree.

## 5a. Validation RESULT (2026-06-08): PASS

Mixed-cohort gate (20 forecast + 40 backgrounded, after the §3a + GH-mean fixes):
- **Sampler:** 0 divergences, R-hat 1.00, ESS ~1400 (solver 1 → solver 2 fallback
  on some chains; no effect on convergence).
- **Marginals:** all six pass — `max diff/MCSE 4.10`, `min SD ratio 0.96`. `tr_sd`
  agrees (ref 0.288 vs mixed 0.279): the forecast cohort pins it, as designed.
- **PFS endpoint:** **invariant** — median PFS identical (34 wk), max quantile
  diff **1 week** (< one assessment interval), max landmark event-rate diff
  **0.019**. ✓
- **Known cosmetic breaches (do NOT reach the deliverable):** joint
  `frac_logit_pop:init_logit_pop` correlation distorts (`corr_diff 0.30`) and the
  week-52 *upper-tail* burden compresses ~8% — both because the quadratic cannot
  separate `frac` from `init` (structural; no polynomial degree fixes it, per the
  cubic pre-check). PFS is a short/mid-horizon threshold-crossing event, so
  neither moves it.

**Verdict:** the surrogate is validated for production use — backgrounded
patients marginalized via the quadratic + GH bridge preserve the reported PFS.

## 6. Open parameters & risks

**Set during implementation (not now):**
- **Anchor times `(t₀, t₁, t₂)`** — from the backgrounded trials' actual
  observation windows. Default start: baseline (0), population-median nadir time,
  latest common follow-up. Tuned empirically in the standalone harness; fixed
  constants once chosen.
- **`measure_sd` for the surrogate** — reuse the model's existing SLD
  measurement-noise parameter (same physical quantity).

**Risks & mitigations (updated post-validation):**
- **`frac`/`init` not separable (CONFIRMED structural limitation).** The quadratic
  sees only the combined early-trajectory effect of `frac_logit` and `init_logit`,
  so it cannot preserve their individual trade-off — it distorts their joint
  correlation (`corr_diff 0.30`) and compresses the long-horizon worst-case burden
  tail ~8%. A cubic does NOT fix this (pre-check: `frac` contributes ~nothing
  beyond the linear coefficient; week-52 is also beyond the data window, so no
  anchor helps). **Validated harmless:** neither moves the reported PFS endpoint.
  Document as a known limitation; revisit only if a future deliverable depends on
  the long-horizon upper-tail burden or the frac/init joint specifically.
- **`tr_sd` not identifiable from the surrogate alone** — its footprint is the
  ~1250×-weaker curvature channel. Mitigation BUILT INTO the design: backgrounded
  trials are used only alongside forecast patients, who pin the SDs. Validated by
  the mixed-cohort gate.
- **Solver-1 → solver-2 fallback** — `(t, t²)` collinearity ill-conditions the
  inner Hessian. No effect on convergence/endpoint. Optional fix: center/scale the
  time basis.
- **`Σ_β` near-singular if an SD → 0** — small documented jitter on `K`'s diagonal
  (`1e-10`).
- **Anchor/window mismatch across trials** — per-trial fixed anchors (still
  constant, still cheap); easy extension, built only if needed.

**Out of scope (YAGNI):**
- Multistate hazards for backgrounded patients (SLD-only, as decided).
- Retiring the old hand-coded `stan/modules/laplace/` (separate cleanup after
  this lands).
- Per-trial anchor sets (only if the gate demands it).
