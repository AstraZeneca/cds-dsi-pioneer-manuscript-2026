# Laplace Surrogate — Performance Investigation (2026-06-12/13)

Measured investigation into why the **marginalized** publication fit
(`tumor_ssls_res_posterior_full`, `enable_surrogate=TRUE`) is slow, and whether it
can be made "scalable + faster" than the explicit all-patient fit. Side/exploratory
work — the MCMC (explicit) version carries production. **Conclusion: as formulated,
the surrogate is slower than the explicit fit. The only lever that touches the
actual (d=4) production fit is within-chain `reduce_sum` parallelism over the 419
block-diagonal patient solves (~4–5× wall-clock on the Laplace component, no math
change), and even that is worth building only if marginalization improves mixing —
which is unmeasured and predicted to fail. Run the ESS/sec gate first.**

> **2026-06-13 correction.** An earlier version of this doc named an *analytic
> linear-Gaussian marginal for the d=2 burden path* as "the only surviving lever."
> That is wrong for the fit that ships. The publication config is **d=4**
> (`enable_frailty=TRUE`, `publication_targets.R:253` → `surrogate_d=4`,
> `transformed_data.stan:38`), where the burden latents `b1,b2` are coupled through
> `standardize_log_burden`/`standardize_velocity` into `exp(loghaz)` per-week hazard
> loops *inside the same per-patient Laplace block* (`surrogate.stanfunctions:141–174`).
> That joint marginal is non-Gaussian — **no closed form exists**, regardless of how
> few patients are censored. The analytic marginal (and its aliases: analytic-Hessian
> hook, sufficient-stats/clustering, EP/VI alternatives) only applies to a d=2,
> SLD-only config that is **not** production. See "Ranked plan" below.

## Setup
- Split: 78 forecast (lilly_cxcr4) fit explicitly + 419 background (amgen_darbe)
  marginalized via `laplace_marginal_tol`. `surrogate_d=2` (burden) or 4 (+frailty).
- Reference: surrogate-OFF all-497 explicit fit (job #1945) **succeeded in 9.4 h**
  (33,817 s; 300 warmup + 500 sampling × 4 chains). This is the baseline to beat.

## Measured root cause
The surrogate `laplace_marginal_tol` call costs **~1.4 s per gradient evaluation**
and is **~75% of total runtime** (Stan `profile("surrogate_laplace")`: 458 s of 609 s;
1.59 s/call over 288 calls). A real fit needs ~36 gradient calls/iter → ~76 s/iter →
~17 h/chain — **slower than the 9.4 h explicit fit it was meant to accelerate.**

## What was ruled out (each by measurement, not reasoning)
| Hypothesis | Test | Verdict |
|---|---|---|
| Bad init / Pathfinder seed | seeded from surrogate-OFF all-497 PF | stall persists — not init |
| Patient-level frailty `Sigma_u` | `enable_frailty=FALSE` → d=2 | stall persists — frailty irrelevant |
| Burden bridge `Sigma_beta` singular | `[BRIDGE-DIAG]` eigenvalue print | never fired — `Sigma_beta` healthy |
| `max_treedepth` | clean td=6 8-iter run | 76 s/iter — treedepth only scales call count |
| Newton solver iterating to cap | `surrogate_max_num_steps` sweep {100,20,5} | **flat: 1.388/1.373/1.339 s/call** — solve converges in ≤5 steps |

The `block_matrix_sqrt` Schur failures and forecast-SLD `-nan` warnings are
**recoverable** (a run completed despite them) — not the cost.

## Why it's structural
The cost is the marginalization linear algebra itself: per leapfrog, for each of 419
patients, an inner Newton solve + 2×2 factorization + log-det, **reverse-mode
autodiffed**. It is block-diagonal/per-patient (Stan already exploits
`hessian_block_size=d`), so it's **O(n_bg), already optimal** — batching/restructuring
gains nothing. Laplace's per-gradient cost is intrinsically ~5–10× an explicit
likelihood eval (forward+reverse once). Laplace pays off only when latents are a
**sampling liability** (high-dim-per-unit, or funnels that wreck explicit HMC). Ours
are 419 units × **2 well-identified, NCP'd burden latents** — Laplace's worst regime:
full inner-solve tax, ~no mixing benefit, because explicit HMC handles 2 latents fine.

## Why the analytic marginal does NOT help production (the d=4 wall)
The tempting lever is a closed-form marginal. Background SLD censoring is low —
**6.3% of visits below-LOD, 10.3% of patients with ≥1 censored visit** — so ~90% of
background patients are fully observed, and *for a d=2 SLD-only patient* the model is
exact Bayesian linear regression: the marginal over `(b1,b2)` is one `multi_normal_lpdf`
with covariance `X·Sigma_beta·Xᵀ + σ²I`, no Newton solve. **But production is d=4, and
the closed form does not exist there.** With `enable_frailty=TRUE` the same `b1,b2`
flow through `standardize_log_burden`/`standardize_velocity` into two `exp(loghaz)`
per-week hazard loops, alongside the frailty latents `u01,u03`, **all in the same
per-patient Laplace block** (`surrogate.stanfunctions:141–174`). `exp` of a
quadratic-in-latent has no finite sufficient statistic; you cannot peel the SLD piece
into a `multi_normal` while `(b1,b2)` remain coupled into the non-Gaussian hazard. So:
- The analytic marginal gives **exactly 1.0× on the d=4 fit that ships.** It only
  helps a d=2 SLD-only ablation, and even there is Amdahl-capped to ~parity end-to-end
  (~2.0–2.5× on the surrogate component → ~6–7 h/chain vs the 9.4 h explicit fit).
- Its aliases all collapse to the same dead end: an **analytic-Hessian hook** (Stan
  2.39 `laplace_options` exposes no derivative callback — verified in
  `laplace_marginal_density_estimator.hpp`); **sufficient-stats / design-clustering**
  (each bg patient carries a distinct static baseline + event week, so even identical
  visit grids don't share a block); **EP/VI/INLA** (research-grade, would only
  *approximate* what is already exact at d=2).
- **Hoisting the prior** (`Sigma_beta`, `beta_pop`) to `transformed_data` is also out:
  `surrogate_bridge` rebuilds them every leapfrog from the live sampled
  (`tr_loc_pop`, frac/init SDs), so they are θ-dependent. Only the θ-independent
  design sufficient stats (`X'X`, `X'y`, inverse-Vandermonde) are precomputable, which
  does not remove the per-call cost.

## Warm-starting the inner solve is also dead
Seeding `theta_0` from the previous leapfrog step cannot help: the measured
`max_num_steps` sweep already proves the solve converges in ≤5 Newton steps, so the
iteration loop is only ~3.5% of the 1.4 s — the other ~96% is fixed Hessian
assembly + factorization + log-det + the reverse-mode/third-order autodiff sweep,
which a better initial guess does not touch. It is also **infeasible** in Stan 2.39
without forking C++: `surrogate_theta_0` is a `transformed_data` constant passed
unchanged each call, and there is no cross-gradient-call persistence of `theta*` in
`laplace_marginal_density_estimator.hpp`.

## Open question gating all of it
Whether marginalizing these latents **improves mixing at all** is UNMEASURED. We
compared wall-clock/iter, not effective-samples/sec. If the explicit fit mixes fine
(it converged in 9.4 h), a faster marginal might match but not beat it. **Decisive next
step if resumed: measure ESS/sec, explicit vs surrogate, at n_bg=419.** If
marginalization doesn't raise ESS/iter, no marginal (analytic or iterative) is worth
building.

## Ranked plan (from the 2026-06-13 multi-agent investigation)
8 efficiency angles investigated then adversarially refuted against the measured
findings above; **only 2 survived.** Sequence:

1. **ESS/sec gate — DO FIRST (low effort, no new Stan code).** The one measurement
   nobody has taken; it gates everything. Add a temporary publication tribble pair
   identical except `enable_surrogate = c(TRUE, FALSE)`, both forced to
   `warmstart=FALSE`, `metric_files=NULL`, same seed, short budget (≈150 warmup + 300
   sampling × 4). For each fit compute `posterior::ess_bulk` on the **shared**
   population vector only — `tr_loc_pop`, `frac_logit_loc_pop`, `init_logit_loc_pop`,
   `measure_sd_sld`, plus `matches('^(time_invariant_coef_qr|time_varying_coef|log_pop_lambda)')`
   — and divide by **sampling-only** seconds (`fit$time()$chains`, exclude warmup).
   **Decision rule:** the surrogate must clear the measured ~8× per-iter wall handicap
   (76 s/iter vs ~10.6 s/iter), i.e. surrogate ESS/iter ≥ ~8× explicit on the binding
   (min) shared param, with no shared-param `Rhat > 1.01` or divergences in either fit.
   `<2×` = clear fail (STOP, ship explicit); `>8×` = clear pass; in-between =
   inconclusive (likely not worth the build cost). **Caveat:** a positive result is
   partly the surrogate's smoother *approximate* (quadratic + GH-bridge) geometry, not
   pure latent removal — so a marginal positive is not a green light by itself.
   *Predicted outcome: NEGATIVE* (419 × 2 well-identified NCP'd latents are Laplace's
   worst regime; explicit HMC already converged in 9.4 h). Feasibility caveat: even a
   "short" surrogate fit is ~9.5 h/chain of cold-start warmup at 76 s/iter.

2. **`reduce_sum` within-chain parallelism — only if (1) passes.** The *only* lever
   that attacks the d=4 production fit without changing the math. The 419 per-patient
   solves are block-diagonal and independent, so splitting the monolithic
   `laplace_marginal_tol` into per-chunk calls inside a `partial_sum` functor divides
   the irreducible O(n_bg) work across cores. Contradicts no measured finding — "O(n_bg)
   already optimal" is about *FLOP count*; this attacks *wall-clock*, never measured.
   The machinery exists (`multistate.stanfunctions:399` already runs `reduce_sum`;
   `stan_threads=TRUE` wired in `util.R`). Realistic ~4–5× on the 75% Laplace component
   (sub-linear TBB scaling, per-chunk overhead, only ~8–12 threads/chain since 4 chains
   share 48 cores) → ~17 h/chain toward ~6.4–7.5 h, i.e. parity-to-modest-win over 9.4 h.
   **First step (de-risk before any production refactor):** a minimal Stan smoke test
   mirroring `stan/experiments/reduce_sum_toy.stan` that calls `laplace_marginal_tol`
   inside a `reduce_sum` `partial_sum` on a 4-patient d=2 synthetic problem; compile
   with `cpp_options=list(stan_threads=TRUE)`; confirm it (a) compiles and (b) matches
   the monolithic target to numerical tolerance. The central unknown is whether the
   2.39 embedded-Laplace nested AD tape composes inside `reduce_sum`'s per-worker
   `nested_rev_autodiff`.

3. **Free diagnostics + warmup-metric harvest — fold into runs you make anyway.**
   Density-invariant, zero-risk, low effort, but capped at parity ("match, never
   beat"): the td=6 run already cost 76 s/iter at low treedepth, so calls/iter is near
   its floor. On the rank-1 short fit, set `save_metric=TRUE` and inspect the CmdStan
   `n_leapfrog__`/`treedepth__` columns to *measure* actual calls/iter and whether
   treedepth hits the cap (metric could help) or sits low (lever dead). Reuse the run;
   do not launch a dedicated job.

**Dead ends (do not pursue):** analytic d=2 marginal and all its aliases (d=4
non-Gaussian — see above); analytic-Hessian injection (no 2.39 API hook); inner-solve
warm-start (≤3.5% of cost + infeasible without forking C++); precomputing the
θ-dependent bridge; any "batching" to beat O(n_bg) (patients are genuinely independent
— parallelism divides the work but cannot make FLOPs sub-linear).

## Tunable knobs added during this work (kept — useful, harmless when surrogate off)
- `enable_frailty` (publication tribble) — gates the Phase-2 patient-level frailty.
- `surrogate_max_num_steps_in` / `surrogate_tolerance_in` / `surrogate_solver_in`
  (data, default sentinel 0 = Stan builtin) — solver knobs tunable without recompile.
- `profile("surrogate_laplace")` block in `likelihood.stan` — for re-measuring cost.

See also memory: `laplace-marginalized-stall`, `laplace-surrogate-toggle`.
