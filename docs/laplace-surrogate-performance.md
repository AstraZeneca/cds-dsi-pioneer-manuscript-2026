# Laplace Surrogate — Performance Investigation (2026-06-12/13)

Measured investigation into why the **marginalized** publication fit
(`tumor_ssls_res_posterior_full`, `enable_surrogate=TRUE`) is slow, and whether it
can be made "scalable + faster" than the explicit all-patient fit. Side/exploratory
work — the MCMC (explicit) version carries production. **Conclusion: as formulated,
the surrogate is slower than the explicit fit; the only real lever left is an
analytic linear-Gaussian marginal for the d=2 burden path, and even that is only
worth pursuing if marginalization improves mixing (unmeasured).**

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

## The only surviving lever: analytic linear-Gaussian marginal (d=2)
Background SLD censoring is low — **6.3% of visits below-LOD, 10.3% of patients with
≥1 censored visit** — so ~90% of background patients are fully observed. For a
fully-observed patient the d=2 model is **exact Bayesian linear regression**: the
marginal over `(b1,b2)` is closed-form, one `multi_normal_lpdf` with covariance
`X·Sigma_beta·Xᵀ + σ²I`. No Newton solve, no `block_matrix_sqrt`. We pay the iterative
tax only because `laplace_marginal_tol` is general-purpose and doesn't know the d=2
core is conjugate. A bespoke analytic marginal would cut the ~1.4 s/call constant
(~5–10×). Caveats: (a) only helps d=2, not the d=4 frailty/hazard path (`exp(loghaz)`
is genuinely non-Gaussian); (b) the ~10% censored patients need the iterative path or
a small approximation; (c) **still O(n_bg) and still autodiffed** — wins the constant,
not the scaling.

## Open question gating all of it
Whether marginalizing these latents **improves mixing at all** is UNMEASURED. We
compared wall-clock/iter, not effective-samples/sec. If the explicit fit mixes fine
(it converged in 9.4 h), a faster marginal might match but not beat it. **Decisive next
step if resumed: measure ESS/sec, explicit vs surrogate, at n_bg=419.** If
marginalization doesn't raise ESS/iter, no marginal (analytic or iterative) is worth
building.

## Tunable knobs added during this work (kept — useful, harmless when surrogate off)
- `enable_frailty` (publication tribble) — gates the Phase-2 patient-level frailty.
- `surrogate_max_num_steps_in` / `surrogate_tolerance_in` / `surrogate_solver_in`
  (data, default sentinel 0 = Stan builtin) — solver knobs tunable without recompile.
- `profile("surrogate_laplace")` block in `likelihood.stan` — for re-measuring cost.

See also memory: `laplace-marginalized-stall`, `laplace-surrogate-toggle`.
