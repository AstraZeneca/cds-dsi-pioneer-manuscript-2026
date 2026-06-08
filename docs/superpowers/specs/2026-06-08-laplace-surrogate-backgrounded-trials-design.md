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
- Normalization fixes `β_i0`'s population mean at 0 (burden normalized to
  baseline), but `β_i0` is kept as a free per-patient latent to absorb baseline
  measurement noise → **3 latents per patient**, matching `hessian_block_size=3`.

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

**Covariance — Jacobian propagation.**
Per-patient spread comes from `{tr_sd, frac_sd, init_sd}`. Propagate through the
same anchor map via its Jacobian `J = ∂β_pop/∂θ`:

```
Σ_β = J · diag(tr_sd², frac_sd², init_sd²) · Jᵀ
```

`Σ_β` is the surrogate's random-effect covariance — the **linearized image** of
the true latent covariance — so patient-to-patient variability the trial sees is
transmitted to `tr_sd`/`frac_sd`/`init_sd`, not just the mean trajectory. `Σ_β`
becomes the prior covariance `K` passed to `laplace_marginal_tol` (non-identity).

- **Jacobian source:** **analytic** `∂g/∂θ`. NOTE — the original intent was
  Stan's built-in autodiff, but Stan has *no callable autodiff Jacobian of a user
  function* (the `jacobian` block / `jacobian +=` are for custom-transform log-det
  adjustments, not a returnable matrix). Since `g = log_sum_exp(a, b)` is closed
  form, `∂g = w_dec·∂a + w_gro·∂b` (softmax weights) is elementary, and
  `J = Vinv · [∂g(t_k)/∂θ]_k`. A finite-difference unit test guards the
  hand-derivation against typos.
- **`Σ_β` rank:** `J·diag(·)·Jᵀ` with `J` 3×3 is generically full-rank → `K` is
  PD as Laplace wants. Degrades gracefully toward rank-2 as any SD → 0.

Why fixed anchors compose well with the tight bridge: `V⁻¹` is constant and `J`
is a small closed-form Jacobian. Nadir-tracking anchors would make `V` itself
depend on θ and `J` gain terms through the moving anchor times — the fragile AD
that sank the old hand-coded solver. Both earlier decisions (tight bridge, fixed
anchors) are what keep this section tractable.

## 4. Integration into `sf-ssls-lfo.stan`

Routing already exists: `background_patient_idx` / `n_background_patients`, and
patient NCP params sized by `n_forecast_patients`, so backgrounded patients
currently have no explicit latents and contribute nothing to the LFO likelihood.
The surrogate slots into exactly that gap.

- **Forecast patients:** unchanged. Full bi-exponential `sf_log_space_obs`,
  explicit NCP latents, full trajectory — everything we report is untouched.
- **Backgrounded patients:** new
  `target += laplace_marginal_tol(surrogate_ll, (θ_pop args + obs),
  hessian_block_size=3, K_fn=Σ_β, ...)`.
- **Solver:** because the surrogate marginal is genuinely log-concave
  (linear-Gaussian), use **solver 1** (Cholesky of the PD Hessian) — the fast,
  well-conditioned path — NOT the solver-3 fallback the bi-exponential needed
  merely to survive. If solver 1 misbehaves, the surrogate isn't log-concave as
  designed → immediate red flag.

New Stan functions in a small module (e.g. `stan/modules/laplace_surrogate/`):
- `surrogate_anchor_betas(θ_pop)` → `β_pop`
- `surrogate_bridge_cov(θ_pop, sds)` → `Σ_β` (via `jacobian`)
- `surrogate_ll(β, θ_pop, obs, ...)` — the marginalized likelihood (quadratic
  log-burden Gaussian + `normal_lcdf` LOD tail)

The old hand-coded `stan/modules/laplace/` stays untouched and unused; retire it
only after this lands (separate cleanup, §6 out-of-scope).

Blast radius: exactly one new code path, reachable only by patients we never
report, rejoining the model at `θ_pop`. Population params, priors, hierarchy,
and the entire forecast path are shared verbatim.

## 5. Validation gate

Same structure as the prior NO-GO experiment, retargeted at the population
posterior (the only thing backgrounded trials affect). Run **before** any
integration into the real model.

**Standalone harness** (`stan/experiments/laplace_surrogate_test.stan` +
`r/experiments/test_laplace_surrogate.R`):
- Synthetic data: a few backgrounded trials' worth of patients simulated from
  the **true bi-exponential** (test the surrogate against the real generative
  process, not its own assumptions), plus censored-below-LOD visits.
- Two fits: **(A)** full-HMC everywhere (reference); **(B)**
  surrogate-marginalized for the backgrounded patients.

**Gate** — surrogate fit clean FIRST, then population agreement:
1. **Sampler health (B):** solver 1; require `pct_divergent < 1`,
   `max_rhat < 1.01`, `min_ess_bulk > 400`. Solver 1 working *is* the design
   claim — if it doesn't, the surrogate isn't log-concave as intended.
2. **Population agreement, joint:** compare
   `{tr_loc_pop, frac_logit_pop, init_logit_pop, tr_sd, frac_sd, init_sd}`
   between A and B. Means within a few MCSE; SDs not understated (B mustn't be
   falsely confident); and the **correlation structure** of the joint population
   posterior preserved. The joint check subsumes target-forecast invariance: the
   target forecast is a nonlinear function of the full joint population
   posterior, so checking the joint catches correlation/tail drift a
   marginal-only gate would miss.

**Honest-gate discipline (carried over):** gate on divergences/R-hat/ESS FIRST;
never let a low-ESS fit fake a PASS through inflated MCSE.

**Optional, not a gate:** target-trial forecast predictive identical A vs B —
belt-and-suspenders only. We do NOT gate on backgrounded-trial predictive (that
would test the very thing we deliberately approximate).

## 6. Open parameters & risks

**Set during implementation (not now):**
- **Anchor times `(t₀, t₁, t₂)`** — from the backgrounded trials' actual
  observation windows. Default start: baseline (0), population-median nadir time,
  latest common follow-up. Tuned empirically in the standalone harness; fixed
  constants once chosen.
- **`measure_sd` for the surrogate** — reuse the model's existing SLD
  measurement-noise parameter (same physical quantity).

**Risks & mitigations:**
- **Quadratic too stiff for long post-nadir regrowth** — a single parabola can't
  track a sharp nadir + sustained exponential regrowth over a very long window.
  The gate will expose this as population bias; pre-agreed fallback is a
  cubic/4-knot spline (still linear-in-coefficients → Laplace stays exact; only
  latent count and one derivation grow). We do NOT silently accept a biased fit.
- **`Σ_β` near-singular if an SD → 0** — degrades to rank-2; solver 1 may
  complain. Mitigation: small documented jitter on `K`'s diagonal.
- **Anchor/window mismatch across trials** — if backgrounded trials have very
  different follow-up lengths, one fixed anchor set may not suit all. Mitigation:
  per-trial fixed anchors (still constant, still cheap) — easy extension, built
  only if needed.

**Out of scope (YAGNI):**
- Multistate hazards for backgrounded patients (SLD-only, as decided).
- Retiring the old hand-coded `stan/modules/laplace/` (separate cleanup after
  this lands).
- Per-trial anchor sets (only if the gate demands it).
