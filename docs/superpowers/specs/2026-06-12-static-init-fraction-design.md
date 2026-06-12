# Static Initial-Fraction Compartment — Design Spec

**Branch:** `karim/process-noise`
**Date:** 2026-06-12
**Status:** Approved design, pending implementation plan

## Motivation

The two-component log-space Stein-Fojo state-space model (`stan/tumor/sf-ssm-log-space.stan`)
decomposes baseline tumor burden into a *decreasing* compartment (treatment effect) and a
*growing* compartment (progression/resistance). The decreasing compartment always decays to
zero, so **any** burden that persists at long horizons must be attributed to the exponential
growth compartment. A patient with genuinely stable disease has no structural way to be
represented as flat — the model is *forced* to assign them a positive growth rate, which then
extrapolates to unbounded ("exploding") SLD at long forecast horizons.

This is a known, documented property of the biexponential family. The pharmacometric
literature (Stein 2008; Bruno 2023; Benzekry 2014; CPT:PSP 2025, `10.1002/psp4.70095`)
confirms the regrowth term `e^(kg·t)` carries a constant, unbounded rate with no carrying
capacity or rate attenuation, and that long-horizon extrapolation of such models produces
outlier predictions. The standard deterministic fixes (Gompertz growth-rate decay,
logistic carrying capacity, Simeoni exponential→linear switch) are documented, but **no
state-space / AR(1) / mean-reverting formulation of a time-varying growth rate exists in the
verified literature** — and no "static compartment" exists in the biexponential family. This
extension is therefore novel.

### Goal of this branch (Step 1 of 2)

Add a third **static** (zero-growth) compartment to the initial-state split. This relieves the
structural forcing that makes *stable-disease patients* explode, by letting persistent burden
sit in a flat compartment instead of being coerced into the growth compartment.

**Scope boundary:** This does *not* bound the growth rate. Genuine progressors (positive growth
fraction) will still grow exponentially. Whether residual progressor explosion remains a
practical problem is to be **measured** after this change; growth-rate decay (Step 2) is a
deferred follow-up, not part of this branch.

## Mathematical Parameterization

Baseline burden `B₀` splits three ways via a **nested conditional-logit** parameterization:

```
π_decrease = inv_logit(init_logit_loc)              # unchanged meaning: decrease vs. rest
rest       = 1 - π_decrease
π_static   = rest · inv_logit(init_logit_static)    # NEW: static vs. growth among the rest
π_growth   = rest · (1 - inv_logit(init_logit_static))
# π_decrease + π_static + π_growth = 1
```

Initial log-states:

```
init_log_decrease = log(π_decrease) + log(B₀)
init_log_static   = log(π_static)   + log(B₀)      # NEW — constant in time (rate ≡ 0)
init_log_growth   = log(π_growth)   + log(B₀)
```

SLD at any time:

```
log_SLD(t) = log_sum_exp(state_decrease(t), init_log_static, state_growth(t))
```

The static term carries no time argument because its rate is zero — its state is constant.

### Why nested conditional logits (not a flat 3-category softmax)

The nested form keeps `init_logit_loc` with its **exact current meaning** (decrease vs.
everything else). All existing fitted priors, covariate effects, and the N-level hierarchy on
`init_logit_loc` remain interpretable and unchanged. The 2-component model is recovered exactly
when `init_logit_static → -∞` (π_static → 0). A flat softmax would re-mix the existing
parameter's meaning and force re-derivation of its priors. The nested form adds exactly one new
orthogonal linear predictor.

### Hierarchy scope for the static fraction (decided)

`init_logit_static` gets a **population intercept only** — one global
`init_logit_static_loc_pop` scalar, no covariates and no random effects. Rationale:
maximally identifiable, easiest to validate recoverability. Hierarchy (covariates, trial-level
RE) can be added later if the data supports it.

## Code Surface

The model has a **single combine funnel**: every observation, replication, and forecast computes
SLD through `calc_log_burden_mean` (`stan/modules/state_space/sf.stanfunctions:187`). The two
components are never combined anywhere else. This makes the change isolated and flag-gated.

### A. Combine function — `calc_log_burden_mean` (arity-polymorphic)

`stan/modules/state_space/sf.stanfunctions:187`. Today asserts `cols == 2` and does a 2-arg
`log_sum_exp`. Make it handle 2 OR 3 columns, discriminating on column count:

```stan
// 2 cols → decrease/growth (static absent); 3 cols → decrease/static/growth.
// Stan's log_sum_exp accepts a row_vector, so the 3-way case reduces each row:
//   for (v in 1:rows(patient_states)) out[v] = log_sum_exp(patient_states[v]);
// The 2-way path keeps the existing vectorized 2-arg form for performance.
```

The exact reduction (per-row loop vs. a vectorized helper) is an implementation detail for the
plan; what the spec fixes is that the function discriminates on `cols(patient_states)` and the
2-column branch is byte-for-byte the current code. This is not a backward-compat alias — it is
one function genuinely handling both shapes. PSA / pioneer and any model with the static
compartment disabled keep passing a 2-column matrix and hit the identical old path.

### B. State array — `states_full_grid` (2 → 3 slices)

`stan/modules/state_space/transformed_parameters.stan:59`. Change
`array[2] matrix[...]` to `array[enable_static_init ? 3 : 2] matrix[...]`. The static slice
(index 3) is filled with `init_log_static_patient` broadcast across all time columns — constant,
no rate, no cumsum. The four existing process-noise branches are **untouched** for components 1
and 2; the constant third slice is appended once, outside the branching. The per-patient
`states` matrix used by the observation model gains a third column when the flag is on.

### C. `init` module (8-step "Adding Module Parameters" checklist)

1. `stan/modules/init/flags.stan` — `int<lower=0,upper=1> enable_static_init;`
2. `stan/modules/init/parameters.stan` — `real init_logit_static_loc_pop;` (population scalar only)
3. `stan/modules/init/transformed_parameters.stan` — compute the nested logits →
   `init_log_decrease_patient`, `init_log_static_patient`, `init_log_growth_patient`. When the
   flag is off, `init_log_static_patient` is absent / size 0 and the existing two-way computation
   is unchanged.
4. `stan/modules/init/hyperparams.stan` + `stan/modules/init/priors.stan` — prior on
   `init_logit_static_loc_pop` (weakly informative; default should put π_static modest, e.g.
   centered so static is a minority of the non-decreasing burden).
5. `stan/modules/init/generated_quantities.stan` — expose π_static (and optionally the three
   fractions) for diagnostics.

### D. R-side wiring

- `r/priors.R` — default prior for `init_logit_static_loc_pop`.
- `r/initializers.R` — init value for the new parameter.
- Targets pipeline (`targets/publication_targets.R` / `pub_targets.R`) — set
  `enable_static_init` flag in the stan-data assembly for the relevant variant(s).
- Stan-data assembly (`prepare_tumor_stan_data` path) — pass the flag through; default `0`.

### E. Tests

- Existing testthat fixtures that assemble `init` stan-data must default `enable_static_init = 0`
  so the current suite is unaffected.
- New tests covering the 3-way path: simplex sums to 1, recovers 2-way when flag off, static
  state is constant in time.

## Validation Strategy

Identifiability is the primary risk: over short follow-up, *static*, *slow-decrease*, and
*slow-growth* trajectories look alike. The validation plan is substantive, not an afterthought.

1. **Recoverability simulation (gating test).** Simulate data from a known 3-way split (e.g.
   π_static = 0.3), fit, confirm the posterior recovers the static fraction without collapsing it
   into the growth/decrease compartments (prior-predictive + simulation-based calibration). If the
   static fraction is not recoverable at realistic follow-up lengths, we learn that *before*
   spending an MCMC run on real data. Non-negotiable because there is no literature precedent for
   this structure.
2. **Backward-compatibility test.** With `enable_static_init = 0`, results must match current
   `main` (2-column path untouched). Existing testthat suite stays green.
3. **Explosion-reduction measurement (the actual goal).** Fit with the static compartment on;
   compare long-horizon forecast SLD distributions against the 2-component baseline. Success =
   fewer patients with runaway extrapolated SLD, and stable-disease patients forecasting flat
   instead of being forced into growth. This step decides whether Step 2 (rate decay) is needed.
4. **Syntax/compile gate.**
   `stanc --include-paths=stan --include-paths=stan/tumor stan/tumor/sf-ssm-log-space.stan`
   must pass before any fit.

## Out of Scope (deferred)

- Growth-rate attenuation / decay over time (Step 2) — only if Step 3 measurement shows residual
  progressor explosion is a practical problem.
- Covariates or random effects on the static fraction — add only if recoverability holds and the
  data supports heterogeneity.
- Any change to PSA/pioneer models — they keep the 2-column path via the flag default.
