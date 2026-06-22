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

### A. Combine function — `calc_log_burden_mean` (extra scalar argument)

`stan/modules/state_space/sf.stanfunctions:187`. The static compartment is **constant in
log-space** (rate ≡ 0), so it never needs to be a propagating state column — it is just one
extra constant term inside the existing `log_sum_exp`. Add an optional per-patient scalar
`static_log_level` (= `log(π_static)`, in the same log-proportion units as the state columns,
before `log(baseline)` is added), defaulting to `negative_infinity()`:

```stan
// Old (preserved exactly via overload with static_log_level = negative_infinity()):
//   log_sum_exp(state[,1], state[,2]) + log(baseline)
// New (3-way): add the constant static level elementwise inside the log_sum_exp.
vector calc_log_burden_mean(matrix patient_states, real baseline, real static_log_level) {
  assert_equal(cols(patient_states), 2);
  vector[rows(patient_states)] lse2 =
    to_vector(log_sum_exp(patient_states[, 1], patient_states[, 2]));
  return log_sum_exp(lse2, rep_vector(static_log_level, rows(patient_states))) + log(baseline);
}
// Backward-compat overload — byte-for-byte the current behavior:
vector calc_log_burden_mean(matrix patient_states, real baseline) {
  return calc_log_burden_mean(patient_states, baseline, negative_infinity());
}
```

Because `exp(-∞) = 0`, `log_sum_exp(lse2, -∞) = lse2` identically — so the 2-component model is
recovered **exactly** when `static_log_level = negative_infinity()`. The state array stays
`array[2]` everywhere; propagation, Kalman, and forecast helpers are untouched. PSA / pioneer
and any model with the static compartment disabled call the 2-arg overload and are unaffected.

### B. Threading the static scalar to the call sites

The per-patient static log-level `init_log_static_patient[p]` (= `log(π_static)`) is passed to
`calc_log_burden_mean` at each call site, traveling the same path as the already-in-scope
`baseline_obs_value` / `baseline_obs_per_patient[p]`. Call sites (all in
`stan/modules/state_space/sf.stanfunctions` unless noted):

- `:870`, `:889` — replication / forecast RNG (uses `baseline_obs_value`)
- `:1009`, `:1013` — mean log obs (uses `baseline_obs_value`)
- `:1128`, `:1135`, `:1142` — per-patient replication / forecast (uses `baseline_obs_per_patient[p]`)
- `stan/tumor/sf-ssls-lfo.stan:213,220,227` — LFO model call sites

When `enable_static_init = 0`, `init_log_static_patient` has size 0 and call sites pass
`negative_infinity()` (or call the 2-arg overload) — exact 2-way recovery. The
`states_full_grid` array and per-patient `patient_states` matrices stay 2-column throughout.

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
