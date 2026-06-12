# Gompertz Growth-Rate Decay — Design Spec

**Branch:** `karim/process-noise`
**Date:** 2026-06-12
**Status:** Approved design, pending implementation plan

## Motivation

The two-component log-space Stein-Fojo state-space model
(`stan/tumor/sf-ssm-log-space.stan`) accumulates the *growth* compartment linearly in
log-space:

```
state_growth(t) = init_log_growth + growth_rate · t
```

so total burden `SLD(t) ≈ exp(state_growth(t)) ∝ e^{growth_rate · t}` — a **constant,
unbounded specific growth rate**. At long forecast horizons this compounds without limit
("forever-exploding SLD"), producing biologically implausible trajectories. This is the
documented core defect of the biexponential family (Stein 2008; Bruno 2023; CPT:PSP 2025,
`10.1002/psp4.70095`): the regrowth term carries no carrying capacity and no rate
attenuation.

The pharmacometric standard fix is **Gompertz growth-rate decay**: the specific growth rate
itself decays exponentially over time, `g(t) = growth_rate · e^{−κ·t}`, so the log-burden
growth arm bends concave-down to a **finite asymptote** rather than rising linearly forever.

### Relationship to the prior (static-compartment) attempt

The earlier branch work added a third *static* (zero-growth) initial compartment
(`enable_static_init`, default OFF). A recoverability simulation
(`r/process_noise/recoverability_sim.R`) **ruled it out**: the static fraction is not
identifiable from the SLD likelihood, because `state_decrease(t) = init − rate·t` aliases a
static compartment as `rate → 0` — there is no temporal signature separating the two.

Gompertz decay does **not** have that problem. The decay parameter κ leaves a distinct,
datable signature: **curvature in the log-SLD growth arm**. κ = 0 produces a straight line;
κ > 0 produces concave-down deceleration. Unlike the static fraction (whose likelihood
gradient was *exactly zero* — structurally aliased by "decrease-rate→0"), κ has a **nonzero
gradient** `∂state_g/∂κ ≈ −growth_rate·t²/2`: it is **weakly identified, not structurally
non-identified**, so the data can move it off the prior where growth curvature is observed.
See "Identifiability" below for the honest precondition. This is the live path; the static
machinery stays committed but OFF.

## Goal

Replace the unbounded linear-in-time growth accumulation with Gompertz decay: the growth
rate attenuates as `growth_rate · e^{−κ_i·t}`, where κ_i is a **per-patient decay rate**
driven by a linear predictor (population intercept now; baseline-covariate slopes and
higher hierarchy levels available later). The whole mechanism is gated by a single master flag
`enable_gr_decay` (propensity-style): default `0` recovers today's model **exactly**, and this
branch sets it to `1` (pop-intercept-only) in `publication_targets.R` to activate the
mechanism. Other projects keep the safe `0` default.

This subsumes and replaces the legacy `growth_lag` / `growth_transition_rate` "slow-the-start"
sigmoid (`get_growth_lag_factor`), which only ever survived in `stan/legacy/` and is not wired
into the live engine. We are not reviving it; the user has confirmed the lag approach never
worked. Gompertz attenuates the *finish* (long-horizon plateau), which is the actual problem.

## Architecture

### The time warp

Define **warped growth-time** as the integral of the decaying specific growth rate:

```
φ_i(t) = ∫₀ᵗ e^{−κ_i·s} ds = (1 − e^{−κ_i·t}) / κ_i
```

Substitute `φ_i(t)` for bare `t` wherever the *growth* component accumulates. Properties:

- **κ_i → 0  ⟹  φ_i(t) → t**   — recovers today's linear model exactly (continuous limit).
- **t → ∞   ⟹  φ_i(t) → 1/κ_i** — log-growth-arm plateaus at `init_log_growth + growth_rate/κ_i`.

The **decrease** compartment is untouched — it already decays to zero, so it never explodes.

### Backward compatibility

`enable_gr_decay = 0` (default) takes a code path that uses the literal `growth_rate · t`
(no `exp`, no warp, zero added cost) — bit-identical to the current model.

The mechanism is an **embedded time-warp**, not an additive term, so the static
compartment's `negative_infinity()`-into-`log_sum_exp` analogy does **not** transfer (that was
additive: `log_sum_exp(x, −∞) = x`). Here the off-path is a literal branch guard at each growth
accumulation site:

```stan
// at every growth-accumulation site:
real warped = enable_gr_decay ? growth_warp(t, kappa_i) : t;   // off-path: warped == t exactly
state_g = init_g + growth_rate * warped;
```

`growth_warp()` is a single shared helper (see "Numerical evaluation" below); the `? :` guard
makes the off-path a true no-op with no `exp` call. The backward-compat gate is the full
existing testthat suite staying green with `enable_gr_decay = 0` defaulted.

### New module: `stan/modules/gr_decay/`

κ borrows the **N-level linear-model plumbing** from `frac` (the closest analog — both shape
the growth process: NCP/CP level modes, QR-space coefficients, flat-index gather), but uses
`tr`'s **log link** (`exp`, not `inv_logit`) and — critically — a **propensity-style
master-flag gating discipline**, not the always-on core pattern of `tr`/`frac`.

**Two distinct module patterns exist in this codebase; gr_decay follows the optional one.**
`tr`/`frac`/`init` are *core* modules — always present, every parameter unconditionally
declared. `propensity`/`multistate` are *optional* modules — a master flag
(`enable_propensity_weighting`, `enable_ms_*`) gates **every** parameter, transform, and prior
so an off model samples and stores nothing. Because Gompertz decay is opt-in, `gr_decay` is an
**optional** module: `enable_gr_decay` is the master gate on *all* of its parameters. The
hierarchy *shape* is copied from `frac`; the *gating* is copied from `propensity`.

Per the 8-step "Adding Module Parameters" checklist in `.claude/rules/stan-guidelines.md`,
the module contains:

| File | Contents (frac-shaped plumbing, tr log-link, propensity-style gating) |
|---|---|
| `flags.stan` | `enable_gr_decay` (master); `enable_pop_cov_gr_decay`; `enable_level_intercept_gr_decay[n_levels]`; `enable_level_cov_gr_decay[n_levels]` |
| `hyperparams.stan` | `gr_decay_log_loc_pop_mean`/`_sd`; per-level intercept/slope SD hyperpriors; QR-coef hyperparams. These are **data**, so per existing-module convention (cf. `init/hyperparams.stan`) they are **always declared with flag-independent shape, unused when off** — no sampling cost, no guard. |
| `parameters.stan` | `gr_decay_log_loc_pop`; QR-space `gr_decay_coef_qr_pop`; unified level intercept/slope SD + raw/cp draws — **every declaration sized `enable_gr_decay ? … : 0`** (outer gate), with `enable_pop_cov_gr_decay ? n_covar : 0` etc. as the inner gates, exactly as `propensity` gates on its master flag |
| `transformed_parameters.stan` | assembles `gr_decay_log_loc_patient` (pop intercept + QR pop covariate effects + level intercepts + level slopes), then `gr_decay_patient = exp(gr_decay_log_loc_patient)`. **Sizing decision (chosen): always `vector[n_forecast_patients]`**, computed inside `if (enable_gr_decay)` and **defaulting to `zeros_vector` when off** — so κ_i = exp(0) = 1 is *never* used (the branch guard picks `t` not the warp when off); the zeros default just keeps the vector well-shaped for the warp sites without a separate conditional-size at every downstream use. This is the propensity pattern (always-shaped, value-inert) and is simpler than multistate's size-0 + guard-every-use. |
| `priors.stan` | population-intercept prior + (guarded) covariate/level priors, identical loop structure to `frac/priors.stan`, all wrapped in `if (enable_gr_decay)` |
| `transformed_data.stan` | flat-index precomputation parallel to `frac` (`n_enabled_groups_*`, `*_flat_idx`, cp/ncp split counts). **Full build-out** (not a stub): the pop-only config makes most counts 0, but the machinery is compiled so enabling levels later is a flag flip. Wrapped so it is inert when `enable_gr_decay = 0`. |
| `generated_quantities.stan` | expose `gr_decay_kappa_pop = exp(gr_decay_log_loc_pop)` and the implied plateau offset `exp(tr_loc_pop + frac_log_growth_pop) / kappa_pop` for interpretability |

**Link function:** κ is modeled on the **log scale** (`gr_decay_log_loc_*`) so κ_i = exp(linpred) > 0 is guaranteed for every patient. The NCP/CP/QR/flat-index plumbing is identical across `tr`/`frac`/`init`/`gr_decay` — **only the inverse link differs** (`exp` for `gr_decay`/`tr`, `inv_logit` for `frac`).

**Initial hierarchy configuration (this branch):** intercept-only, population level only.
`enable_gr_decay = 1` to activate the mechanism, but `enable_pop_cov_gr_decay = 0`,
`enable_level_intercept_gr_decay = rep(LEVEL_MODE_NONE, n_levels)`,
`enable_level_cov_gr_decay = rep(0, n_levels)`. So the linear-predictor assembly reduces to a
scalar broadcast:

```stan
// pop-only reduction (all level/cov terms are zero-length / zeros):
vector[n_forecast_patients] gr_decay_linpred_pop =
  enable_pop_cov_gr_decay ? (Q_covar_design_matrix[forecast_patient_idx, :] * gr_decay_coef_qr_pop)
                          : zeros_vector(n_forecast_patients);
vector[n_forecast_patients] gr_decay_log_loc_patient = gr_decay_log_loc_pop   // scalar broadcast
  + gr_decay_linpred_pop                // zeros here
  + gr_decay_linpred_level_intercepts   // zeros here
  + gr_decay_linpred_level_slopes;      // zeros here
```

i.e. `gr_decay_log_loc_patient == gr_decay_log_loc_pop` for all patients. The covariate/level
machinery is **compiled-but-inert** (not commented out) — turning on baseline-covariate
response later is a flag flip plus a prior, no structural change (exactly how `tr`/`frac` work).

**Two-level gating, explicitly:**
- `enable_gr_decay = 0` (outer gate) → all `gr_decay_*` parameters sized 0, warp sites take the
  `t` branch → bit-identical to current model.
- `enable_gr_decay = 1, enable_pop_cov_gr_decay = 0` (inner gate) → `gr_decay_coef_qr_pop` sized
  0, pop intercept active, covariate response inert.

### Wiring κ into the growth accumulation

The live in-sample forward simulation lives in
`stan/modules/state_space/transformed_parameters.stan`, inside `profile("states")`. There are
**FOUR** mutually-exclusive branches (selected by `enable_patient_process_noise_tr` /
`enable_pop_process_noise_tr` / `need_states_full_grid`). Each accumulates the *growth*
component differently, so each needs the warp applied at the right place. **For the publication
model, process noise is OFF** (`publication_targets.R:371-372`), so branches 3–4 are the live
path there; branches 1–2 matter for sclc/pioneer configurations that enable noise.

Throughout, define the shared helper `growth_warp(t, κ) = φ(t)` (numerics below) and its
*difference* form `Δφ(t_a, t_b) = φ(t_b) − φ(t_a)`. The off-path guard is
`enable_gr_decay ? … : (linear form)` at each site.

**Branch 1 — patient-level process noise** (L62–90; cumsum @ L74):
current `states_full_grid[2][j, 2:] = init_g[j] + cumulative_sum(patient_growth_rate[j, :])`.
The per-step rate carries AR deviations (not constant), so warp at the **increment** level
using **φ-differences** so the cumsum telescopes exactly to the warped trajectory:
```
// decay weight per step, built from φ-differences (NOT a Riemann e^{-κt}·Δt quadrature):
inc[k] = patient_growth_rate[j,k] * (φ_j(t_{k+1}) − φ_j(t_k))      // deterministic part
state_g[j, 2:] = init_g[j] + cumulative_sum(inc_with_noise)
```
Because the deterministic increments are exact φ-differences, with constant rate they
telescope to `growth_rate · φ_j(t)` — **bit-identical** to the no-noise warped branches (this
fixes the cross-branch invariance the naive `e^{−κt_k}·Δt_k` Riemann sum would break, a ~2.5%
gap at κ=0.05 growing with horizon).

**Branch 2 — population process noise** (L91–135; outer-product @ L119–120):
current code exploits a **factorization** — `state_g[i,t] = init_g[i] + baseline_rate[i] ·
pop_cumsum_exp[t]` — valid only because the time factor is patient-independent. Gompertz's
decay weight `e^{−κ_i·t}` depends on **both** i and t, so **the outer-product factorization no
longer holds**. This branch requires a **structural rewrite to a per-patient loop** (like
Branch 1) computing `init_g[i] + cumulative_sum(baseline_rate[i] · pop_noise_factor · Δφ_i)`.
This is a larger change than "multiply each increment" — flag it explicitly in the plan. (Not
on the publication path; needed before sclc/pioneer use pop-noise + Gompertz together.)

**Branch 3 — no-noise full grid** (`need_states_full_grid`, L136–159; @ L145):
current `states_full_grid[2] = init_g * ones + patient_growth_rate[,1] * (time_since_first_visit − 1)`.
The elapsed time here is `(time_since_first_visit − 1)` (weeks since first visit), so the
substitution warps **that** argument:
`→ init_g * ones + patient_growth_rate[,1] .* growth_warp(time_since_first_visit − 1, κ)`
(note **φ(t−1)**, not φ(t) — must match the existing one-week offset or this branch desyncs
from the visit branch).

**Branch 4 — no-noise direct visit** (else, L160–187; @ L184):
current `states[,2] = init_g[j] + patient_growth_rate[j,1] * dt` where `dt = visit_week −
baseline_week`.
`→ init_g[j] + patient_growth_rate[j,1] * growth_warp(dt, κ_j)` (warp the per-visit elapsed
time). This is the live publication in-sample path.

All four reduce to the present code when `enable_gr_decay = 0`.

### LFO forecast path (where the explosion is most visible)

`stan/tumor/sf-ssls-lfo.stan:205` generates forecasts by calling `sf_log_space_trajectory_ncp`
with `growth_lag = negative_infinity(), transition = 1.0`. Those arguments make
`get_growth_lag_factor` return **1** at every step — i.e. the existing `time_varying_factor`
multiplier on `growth_rate` is currently a no-op there. **That multiplier is exactly the slot
for Gompertz `e^{−κ·t}`.** So contrary to the first draft, we *do* touch the helper, but
cleanly: thread κ into `sf_log_space_trajectory_ncp` and have it fold the Gompertz weight into
`time_varying_factor` (replacing the disabled lag factor), again via φ-differences so the
discrete forecast trajectory telescopes to `growth_rate · φ(t)`. κ_i reaches the LFO loop the
same way `static_log_level_per_patient[i]` already does (forecast-local index → unified
`forecast_patient_idx`; see the static wiring at `sf-ssls-lfo.stan:226,233`).
`_lfo_endpoints_generated_quantities.stan` uses cutoff-compact indexing — apply the same
treatment there. **Add a numerical check** that forecast SLD plateaus at
`init_log_growth + growth_rate/κ` as the horizon grows.

The change is confined to the live `transformed_parameters.stan` (4 branches) plus the LFO
helper + forecast paths. The legacy `calc_states` in `sf.stanfunctions` (only `legacy/` callers)
is left alone.

#### Numerical evaluation: the `growth_warp` helper

Single shared function in `sf.stanfunctions`:
```stan
real growth_warp(real t, real kappa) {
  if (kappa * t < 1e-8) {          // small-|κt| branch: guards 0/0 AND keeps dφ/dκ finite for autodiff
    return t * (1 - 0.5 * kappa * t);   // 2-term Taylor: φ ≈ t − κt²/2
  }
  return -expm1(-kappa * t) / kappa;     // expm1 fixes NUMERATOR precision only
}
```
Two distinct numerical concerns, often conflated:
- **`expm1`** addresses *numerator* catastrophic cancellation (`1 − e^{−x}` for small x) — but
  does **not** fix the `0/0` division.
- **The small-`|κt|` series branch** is what guards the `0/0` division *and* keeps the gradient
  `∂φ/∂κ` finite for Stan's autodiff (the raw quotient's derivative is `0/0` at κ=0). Use the
  series **inside** the threshold, not just `expm1`. Threshold constant: `|κ·t| < 1e-8`.

The flag-off path (`enable_gr_decay = 0`) never calls this helper — it takes the literal `t`
branch — so exact κ = 0 costs nothing; the helper's small-κ guard only matters for tiny but
nonzero κ inside the flag-on path. The φ-difference form used in the cumsum branches is
`growth_warp(t_b, κ) − growth_warp(t_a, κ)`.

## Data flow

```
baseline covariates (design matrix, QR) ──┐
                                          ▼
gr_decay_log_loc_pop  +  (QR pop cov)  +  (level intercepts)  +  (level slopes)
                                          │   [all inert this branch except the pop intercept]
                                          ▼
                          gr_decay_log_loc_patient  ──exp──▶  κ_i  (per patient, > 0)
                                          │
   ┌──────────────────┬──────────────────┼──────────────────┬──────────────────┐
   ▼                  ▼                   ▼                   ▼                  ▼
B1 patient-noise   B2 pop-noise      B3 no-noise grid   B4 no-noise visit   LFO forecast
cumsum of          per-patient       growth_rate·       growth_rate·        time_varying_factor
Δφ-increments      loop (factoriz-   φ_i(t−1)           φ_i(dt)             = e^{−κt} via Δφ
+ noise            ation broken;                                            (replaces disabled lag)
                   Δφ-cumsum)
   └──────────────────┴──────────────────┼──────────────────┴──────────────────┘
                                          ▼   (all Δφ-forms telescope to growth_rate·φ(t))
                          state_growth(t)  →  log_sum_exp with state_decrease  →  log-SLD mean
```

## Identifiability — weakly identified, not aliased (the crux)

This is the make-or-break, and the honest statement differs from the first draft.

**Where κ enters the likelihood.** Expanding the warp,
`state_g(t) = init_g + growth_rate·t − (growth_rate·κ/2)·t² + O(κ²t³)`. So κ enters **only as
the curvature (t² and higher) of the growth arm, with coefficient `growth_rate·κ`**. The
correct alias to worry about is therefore **κ vs `growth_rate`** (they co-determine the shape
of the same arm), *not* κ vs the decrease rate as the first draft argued.

**Why this is *weak* identifiability, not the *structural* non-identifiability that killed
static.** The static compartment had a likelihood gradient of **exactly zero** —
`state_decrease = init − rate·t` reproduced a static compartment *identically* as rate→0, so no
amount of data could separate them. κ is different: `∂state_g/∂κ ≈ −growth_rate·t²/2` is
**nonzero** wherever a growth arm is observed. The information is *weak* (scales with `t²` and
with `growth_rate`, which must itself be identified), but it is **present** — so the posterior
**will move off the prior where signal exists** (long-follow-up growers, pooled across the
cohort) and **fall back to the prior where the data is silent**. That is the design behaving
correctly, not failing.

**Honest precondition.** κ is informed only where `growth_rate` is identified **and** follow-up
reaches the regrowth (post-nadir) phase long enough to show curvature. Our own
[[lfo-forecast-regrowth-bias]] finding says that is a minority of patients on SCLC/sclc
data: growth-rate identification is already weak (cor ≈ 0.35; median post-nadir follow-up
≈ 7 wk; only ~25/78 with ≥12 wk). So on the *current* target data κ will be **substantially
prior-influenced** — this is the *expected* operating regime, accepted deliberately (the user's
goal is forecast realism / bounded SLD, not estimating a biological deceleration constant).
Critically, accepting prior influence does **not** preclude the data from speaking: a genuinely
weakly-informative prior lets a real signal pull κ off-center wherever the cohort shows it.

### The role of the prior — weakly-informative, signal can override

`r/priors.R` gets `gr_decay_log_loc_pop_mean` / `_sd` set so that:
- the prior is **centered on a mild positive κ** — implied growth-rate half-life on the order of
  the longest forecast horizon of interest — so that when the data is silent, forecasts plateau
  instead of exploding; and
- the prior is **wide enough that the posterior can contract toward 0** (arms look linear) **or
  toward strong decay** (arms bend). It must not be so tight that it *precludes* the data from
  moving κ — that would defeat the "let the data win where it can" intent.

**Concrete defaults (to be finalized against the sim's threshold result, but pinned so R wiring
and the gate are implementable now):** model `log κ` with
`gr_decay_log_loc_pop_mean = log(0.02)` (κ ≈ 0.02 /week → growth-rate half-life ≈ 35 weeks,
plateau offset `growth_rate/κ` reached over a multi-year horizon) and
`gr_decay_log_loc_pop_sd = 0.75` (95% prior κ ∈ ≈ [0.0045, 0.087] /week — spans
"barely-decelerating" to "plateaus within a year", so the data can move within a wide band).
These are cross-referenced to the sim's κ_true (below). Initializer
(`r/sclc/initializers_fixed.R` — the **production** path used by `publication_targets.R`,
*not* legacy `r/initializers.R`) draws the length-1 `gr_decay_log_loc_pop` via
`as.array(rnorm(1, mean, sd))` when `enable_gr_decay = 1`, length-0 when off (mirroring how the
static param init is conditioned).

### Recoverability simulation — characterize the threshold, not rubber-stamp

The sim's job is **not** a pass/fail rubber-stamp engineered with a visible plateau. It must be
**adversarial against our empirical follow-up**, answering: *given the censoring/nadir-timing
distribution we actually have, how strong must a real deceleration be before the posterior
detectably moves off the prior?*

1. Simulate N patients drawing follow-up length and nadir timing from the **empirical
   distribution** (most censored while still declining or ≤ ~7 wk post-nadir; a minority with
   ≥12 wk), with **per-patient `growth_rate` heterogeneity** (so κ rides on a realistically
   weakly-identified rate), and realistic measurement noise.
2. Sweep a **range** of κ_true (from ≈ prior center up to clearly-decelerating).
3. Assemble stan-data with `enable_gr_decay = 1`, priors from `r/priors.R`; fit
   `stan/tumor/sf-ssm-log-space.stan` short (2 chains, 500/500) at each κ_true.
4. Report, per κ_true: posterior mean/CI for `gr_decay_kappa_pop`, **prior-vs-posterior
   contraction** (does the posterior move off, or sit on, the prior?), and divergences.

**Outcome is diagnostic, not gating:** the κ_true at which the posterior reliably separates from
the prior **is** the answer to "is there signal?" If even strong deceleration can't move the
posterior under empirical censoring, κ is purely a prior regularizer on this data (still useful
for bounded forecasts); if moderate deceleration moves it, the data can speak. Either way the
mechanism ships — the sim tells us *which regime we are in*, it does not veto the work.

### Required deliverable: prior-sensitivity / contraction analysis

Because the prior will materially influence κ on the target data, the branch must ship a
**prior-sensitivity analysis**: vary `gr_decay_log_loc_pop_mean`/`_sd`, and report (a) the
forecast-SLD plateau and (b) in-window fit (e.g. log-lik / LOO) at each setting. This makes the
prior's influence explicit and auditable rather than hidden, and is the honest companion to
shipping a weakly-identified parameter.

## Error handling / edge cases

- **Flag off (`enable_gr_decay = 0`):** literal `t` branch at every site, no helper call,
  bit-identical to current model. The backward-compat gate is the full existing testthat suite
  staying green with the flag defaulted OFF.
- **Tiny κ (flag on):** `growth_warp` small-`|κt|` series branch guards the 0/0 division and
  keeps `∂φ/∂κ` finite for autodiff; reduces to `t`.
- **Cross-branch invariance:** all four branches + LFO use φ-differences, so a constant-rate
  patient gets *identical* `growth_rate·φ(t)` regardless of which branch computes it. A test
  asserts this (it is the invariance the naive Riemann sum would break).
- **Pop process noise + Gompertz together:** the factorized outer-product (Branch 2) is replaced
  by a per-patient Δφ-cumsum loop — a structural rewrite, not a tweak.
- **Patient with single growth visit:** `growth_warp(dt)` with one interval is well-defined.
- **Sized-zero params when off:** every `gr_decay_*` *parameter* is sized `enable_gr_decay ? … : 0`
  (an off model samples/stores nothing extra); *hyperparameters* are data — always declared,
  flag-independent shape, unused when off (matching existing-module convention).

## Testing

- **Stan unit (`growth_warp`):** a tiny test model confirms φ(t)→t as κ→0, φ(t)→1/κ as t→∞,
  the small-`|κt|` series branch, and the φ-difference telescoping identity
  `Σ Δφ = φ(t)` (mirrors `test-stan-calc-log-burden-mean.R`).
- **R unit (hierarchy reduces correctly):** with the pop-intercept-only config,
  `gr_decay_log_loc_patient` is constant across patients = `gr_decay_log_loc_pop`.
- **Cross-branch invariance:** constant-rate patient yields identical `growth_rate·φ(t)` in all
  four `transformed_parameters` branches and the LFO path.
- **LFO plateau check:** forecast SLD → `init_log_growth + growth_rate/κ` as horizon → ∞.
- **Backward-compat (the headline gate):** `enable_gr_decay = 0` → full testthat suite green;
  CmdStanR compiles `sf-ssm-log-space.stan` and `sf-ssls-lfo.stan` with no errors.
- **Recoverability sim (diagnostic, not gate):** characterizes the κ-detection threshold under
  empirical censoring — reports which regime (data-informed vs prior-driven) we are in.
- **Prior-sensitivity (deliverable):** forecast plateau + in-window fit across prior settings.

## Out of scope (this branch)

- Baseline-covariate response on κ (machinery built, flag OFF).
- Higher hierarchy levels on κ (trial/arm/patient REs — machinery built, modes NONE).
- Reviving the legacy `growth_lag` start-delay.
- Stochastic (AR1/OU) time-varying κ. Gompertz is deterministic decay; a stochastic κ is a
  possible later step but is explicitly deferred.
- Logistic/Verhulst carrying-capacity or Simeoni exponential→linear alternatives (considered,
  not chosen — Gompertz is the lowest-surgery fit to this model's existing log-linear growth
  accumulation, reusing the `time_varying_factor` slot. Identifiability of *any* deceleration
  parameter is limited by the same short post-nadir follow-up; Gompertz is not uniquely
  advantaged there, but it is the cheapest to wire and the cleanest weakly-informative
  forecasting regularizer).

## Files touched (summary)

**New:** `stan/modules/gr_decay/{flags,hyperparams,parameters,priors,transformed_data,transformed_parameters,generated_quantities}.stan`

**New:** `growth_warp()` (+ φ-difference helper) in `stan/modules/state_space/sf.stanfunctions`.

**Modified:**
- `stan/tumor/sf-ssm-log-space.stan` — `#include` the new module's fragments (alongside `tr`/`frac`/`init`)
- `stan/modules/state_space/transformed_parameters.stan` — warp **all four** growth-accumulation
  branches (B1 patient-noise Δφ-cumsum @L74; B2 pop-noise **structural per-patient rewrite**
  @L119; B3 no-noise grid φ(t−1) @L145; B4 no-noise visit φ(dt) @L184)
- `stan/modules/state_space/sf.stanfunctions` — thread κ into `sf_log_space_trajectory_ncp`
  (fold Gompertz weight into the existing `time_varying_factor`, replacing the disabled lag)
- `stan/tumor/sf-ssls-lfo.stan` (call @L205), `stan/tumor/_lfo_endpoints_generated_quantities.stan`
  (cutoff-compact indexing) — apply decay on the forecast path; κ_i reaches the loop via the
  same `forecast_patient_idx` discipline as `static_log_level_per_patient`
- `r/priors.R` — `gr_decay_log_loc_pop_mean = log(0.02)`, `_sd = 0.75` + (data) covariate/level hyperparams
- `r/sclc/initializers_fixed.R` — conditional length-1 `gr_decay_log_loc_pop` init (the
  **production** initializer for publication; *not* legacy `r/initializers.R`)
- **Target-file defaulting discipline:** `enable_static_init` is set only in
  `publication_targets.R:389`; sclc/pioneer rely on the data-prep / stan-data default
  for unset flags. `enable_gr_decay` must follow the **identical** discipline — the plan
  verifies how `prepare_tumor_stan_data()` defaults unset flags and sets `enable_gr_decay = 1L`
  (+ pop-only hierarchy config) **only** in `publication_targets.R` for this branch, leaving
  the safe `0` default to flow everywhere else. Enumerate and confirm: `publication_targets.R`,
  `sclc_targets.R`, `pioneer_targets.R`.
- `r/process_noise/recoverability_sim.R` — add the adversarial κ-threshold sim (empirical
  censoring, growth_rate heterogeneity, κ_true sweep, prior-vs-posterior contraction)
- model-spec documentation page — document the Gompertz growth-rate decay

**New tests:** `growth_warp` Stan test (limits + telescoping) + R hierarchy-reduction test +
cross-branch invariance test + LFO plateau check (+ the sim and prior-sensitivity are diagnostics/deliverables).
