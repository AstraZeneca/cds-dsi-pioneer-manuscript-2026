# Design: (level, velocity) multistate burden-coupling basis

**Branch:** `karim/laplace` · **Date:** 2026-06-09 · **Status:** design (pre-implementation)

## 1. Motivation

The multistate (MS) hazard couples burden-driven transitions (0→1 progression,
0→3 dropout/bridge, 0→2 death) to the tumour trajectory through a set of
**time-varying covariate features**. The publication model currently uses
**three** features (`n_time_varying_covar = 3`):

1. **level** — standardized log-burden `g(w)`
2. **log-decrease-rate** — `patient_log_decrease_rate`, a bi-exponential SSM
   *component parameter*
3. **log-growth-rate** — `patient_log_growth_rate`, the other component parameter

Features 2–3 are problematic:

- **No surrogate analog.** The Laplace surrogate that marginalizes background-trial
  tumour latents represents burden as a quadratic `g(w) = b₀ + b₁w + b₂w²`. That
  representation has no `decrease_rate`/`growth_rate` to read off — the rates are
  internals of the *bi-exponential* state, not of the *output trajectory*. This
  makes the forecast/background featurization fundamentally asymmetric and blocks
  a shared `tv_coef` across forecast and borrowed patients.
- **Static in costume.** With `enable_patient_process_noise_tr = FALSE` (the
  publication setting) the rates are *constant per patient*
  (`rep_row_vector(...)`) — baseline-rate covariates wearing a time-varying suit.

**Replacement:** a **(level, velocity)** basis, where velocity is the derivative
of the log-burden trajectory:

| Feature | Meaning | Forecast path | Surrogate (background) path |
|---|---|---|---|
| level | `g(w)` log-burden | `log_sum_exp(state_d, state_g)` (existing) | quadratic `b₀+b₁w+b₂w²` |
| velocity | `g'(w)` rate of change | central difference of latent trajectory | analytic `b₁+2b₂w` |

**Why velocity is the right second feature:**

- It is a derivative of the *output* trajectory, so it **exists in both
  representations** (bi-exponential forecast and quadratic surrogate) and is
  **identical in form** — "d/dw of the log-burden path." The same `tv_coef` is
  therefore meaningful across forecast and background patients.
- It is **linear in the marginalized latents `(b₁, b₂)`** on the surrogate path
  (`b₁ + 2b₂w`), so feeding it into the log-hazard keeps the survival term
  log-concave — the property the joint-surrogate gate validated. The
  bi-exponential rates do not have this property under the surrogate.
- Clean clinical reading: "where burden is (level) + which way it is moving
  (velocity)" — tumour kinetics, akin to growth-rate constants in the oncology
  literature.

**Avoid** features nonlinear in the latents (time-to-nadir, nadir depth, relative
velocity `g'/g`, doubling time) — they break the surrogate's log-concavity.

This change is **§9 step 0** of the joint-surrogate dependency chain
(`docs/superpowers/specs/2026-06-09-joint-surrogate-index-contract.md`) — the
gated prerequisite for the surrogate functor — **but it is also independently
useful to the publication model** regardless of the surrogate, and is specced for
the model's own sake.

### 1.1 Why velocity is admissible inside the Laplace surrogate functor

This is the load-bearing reason the basis was chosen, verified against the
current functor (`stan/modules/laplace_surrogate/surrogate.stanfunctions`,
`likelihood.stan`) and the index contract (§1, §4, §6 O1).

The surrogate marginalizes, per background patient, `theta = [b₁, b₂]` (the
quadratic burden's slope/curvature; intercept `b₀` is pinned —
`surrogate_ll` uses `b0 = beta_pop[1]`, `hessian_block_size = 2`). Stan's
`laplace_marginal_tol` runs an inner Newton solve that requires the functor's
log-density to be **log-concave in `theta`**.

When the contract adds the MS hazard term to the functor it reads, at week `w`:

```
loghaz_T(w) = base_static_T[i,w] + u_T              ← θ-independent + frailty
            + tv_coef_T[1]·std_level(g(w))          ← level coupling
            + tv_coef_T[2]·std_vel(g'(w))           ← velocity coupling
```
with `g(w) = b₀ + b₁·w + b₂·w²` and `g'(w) = b₁ + 2·b₂·w`.

Both `std_level(g(w))` and `std_vel(g'(w))` are **affine in `(b₁, b₂)`**
(constants `b₀`, week `w`, and the standardization scalars are all fixed within
the solve). Therefore `loghaz_T(w)` is affine in `theta`, the survival
accumulation `−exp(loghaz)` is concave (Hessian `−exp(·)·aaᵀ ⪯ 0`), the
event-week `+loghaz` is linear, and the prior `−½θᵀK⁻¹θ` is strictly concave —
so the functor stays log-concave and the Newton solve has a unique mode. **This
is exactly the property the bi-exponential rate features could not provide:**
`rate_d`/`rate_g` are `exp(·)` functions of the bridge's *input* parameters
`(tr_loc, frac_logit)`, which have no expression in the marginalized output
coordinates `(b₁, b₂)` — the term cannot even be written inside the functor,
let alone kept concave. Velocity, as `d/dw` of the *output* trajectory, lives
natively in `(b₁, b₂)`.

**Data vs latent split (important, mirrors the level feature):** the velocity
*value* `g'(w) = b₁ + 2·b₂·w` is **latent** — computed from `theta` inside the
functor body, per Newton step; it is never data. Only the **standardization
constants** `median_velocity_obs` / `iqr_velocity_obs` are `data` (parameter-free
observed-data summaries, §3.3), passed positionally as the functor's data-only
arguments, exactly as `median_log_burden_obs`/`iqr_log_burden_obs` already are
for the level feature. The standardized feature is
`(g'(w) − median_velocity_obs) / iqr_velocity_obs`: latent numerator, data
denominator.

**Division of labour with the contract.** This spec delivers (a) the forecast-path
velocity feature and (b) the `median_velocity_obs`/`iqr_velocity_obs` constants —
which the functor consumes as data. The functor edit that *reads* these
(extending `surrogate_ll` with the hazard term) belongs to the joint-surrogate
contract's implementation, not this spec. The contract's §4 functor signature and
§6 O1 already assume this basis as input; the compatibility above confirms the two
specs agree.

## 2. Scope and the mode flag

**Scope: tumour (publication) model only.** PSA / Pioneer is untouched and
stays on the existing 3-feature basis. The shared feature-builder files branch on
a flag that pioneer leaves `FALSE`. (PSA "velocity" has different units and
dynamics; migrating it is out of scope and deferred.)

**Flag-gated, not hard replacement.** Both bases coexist behind a single boolean,
selectable per fit. The velocity basis is required for the Laplace surrogate, so
it must be switchable on and off.

New flag in `stan/modules/multistate/flags.stan`:

```stan
int<lower=0, upper=1> enable_ms_velocity_basis;  // 0: [level, dec_rate, gro_rate]; 1: [level, velocity]
```

Semantics — **a pure mode switch**, mutually exclusive bases:

| `enable_ms_velocity_basis` | features emitted | `n_time_varying_covar` |
|---|---|---|
| `0` (default, current behaviour) | `[level, log_decrease_rate, log_growth_rate]` | `3` |
| `1` (new, required for Laplace) | `[level, velocity]` | `2` |

The flag flips **both** the feature set and the implied feature count. When the
velocity basis is on, the surrogate-incompatible bi-exponential rate features
never enter the covariate matrix. There is never a `[level, velocity, rate_d,
rate_g]` combination — the rates are exactly the thing being retired in this mode.

**Config contract.** R config sets `enable_ms_velocity_basis` and the matching
`n_tv_covar` **together**: `(TRUE, 2L)` or `(FALSE, 3L)`. A mismatch (flag `TRUE`
with `n_tv_covar = 3L`, or vice versa) is a **config error**. The implementation
adds a defensive check (Stan reject or R-side stop) so the pair cannot silently
disagree — otherwise the coefficient arrays size to the wrong dimension while the
builder emits a different number of features.

## 3. Velocity definition and standardization

### 3.1 Velocity = central difference of the latent log-burden trajectory

The log-burden trajectory `g(w) = log_sum_exp(state_d, state_g)` lives on a weekly
integer grid. Velocity is its central difference:

- Interior: `v(w) = (g(w+1) − g(w−1)) / 2`
- First column: forward difference `v = g(w+1) − g(w)`
- Last column: backward difference `v = g(w) − g(w−1)`

Central difference is chosen because it is **algebraically identical to the
quadratic surrogate's analytic derivative** `b₁ + 2b₂w` (the linear interpolation
of a parabola's secants recovers the tangent at the midpoint). Forecast and
surrogate velocities are therefore consistent **in form** — the property the
joint-surrogate gate requires. They are not bit-identical: the forecast velocity
differences a capped, bi-exponential, possibly noisy trajectory; the surrogate
velocity is the uncapped analytic polynomial derivative. `tv_coef` is shared
model-wide and the contract's §6 O1 / §9 step 2 feature-SCALE check covers this.

The level cap (`fmin(g, 10.0)`, applied in the existing builders to bound the
log-hazard during warmup) is applied to `g` **before** differencing, so a capped
excursion cannot produce a spurious velocity spike.

`log_sum_exp` is retained on the forecast path. This is deliberate and does **not**
affect marginalization: forecast-trial tumour latents are **sampled by HMC**, not
Laplace-marginalized. Log-concavity is only required of the Laplace approximation,
which applies to the *separate* surrogate functor (the quadratic-burden background
path), not to these sampled-latent builders.

### 3.2 Two forecast code paths, both branch on the flag

1. **Dense path** (`stan/_ms_burden_tv_covar.stan`, process-noise ON): velocity =
   central difference across the columns of `states_full_grid` (the materialized
   weekly latent trajectory). Cap-before-difference as above.

2. **Inline path** (`stan/_ms_burden_inline_tv_covar.stan`, process-noise OFF —
   the publication setting): no `states_full_grid` exists. Evaluate the analytic
   burden at three points `g(t−1), g(t), g(t+1)` via three `log_sum_exp` calls per
   step, then central difference. Cheap (3 evaluations/step). Applies to all three
   inline blocks (0→2 dense, 0→3 dense, 0→1 sparse visit-gated).

Both paths, in velocity mode, emit feature 2 = standardized velocity and emit **no**
feature 3. In legacy mode they emit features 2 (decrease rate) and 3 (growth rate)
exactly as today.

### 3.3 Standardization: observed-for-scale, latent-for-signal

Mirrors the existing level feature, where `median_log_sld_obs` / `iqr_log_sld_obs`
are observed-data scalars but the level covariate *value* uses the modeled
`states_full_grid`. The constants anchor units; they do not assert the observed
trajectory is correct.

**New transformed-data constants** `median_velocity_obs` / `iqr_velocity_obs`,
computed in `stan/modules/tumor/transformed_data.stan` (alongside the existing
`median_log_sld_obs` block) from **observed** consecutive-measured-visit log-SLD
deltas:

- For each patient, walk consecutive **measured** visits (mask `sum_tumor_size >
  0`, reusing `patient_visit_pos` + `t_patient_visits`).
- Per adjacent measured pair, compute a **per-week** velocity
  `(log_sld[v+1] − log_sld[v]) / (week[v+1] − week[v])`. Per-week normalization is
  essential: observed visits are irregularly spaced (weeks apart) while the latent
  grid is weekly, so dividing by the week gap puts the observed scale in the same
  units as the latent weekly-grid central difference.
- Collect all per-week deltas across patients; take `median` and `IQR` via the
  existing `quantile(..., {0.25, 0.5, 0.75})` idiom.

These are computed unconditionally in tumour transformed data (like the level
constants) so they are available whenever the tumour model includes the builder,
regardless of which basis is active.

**Aliasing.** `stan/tumor/_tumor_observed_covar_transformed_data.stan` aliases the
burden-agnostic names `median_velocity_obs` / `iqr_velocity_obs` to the tumour-side
values, parallel to the existing `median_log_burden_obs = median_log_sld_obs`
aliasing. (The PSA adapter does not, since PSA stays on the legacy basis.)

**New helper** in `stan/_burden.stanfunctions`:

```stan
real standardize_velocity(real velocity_raw, real median_velocity_obs, real iqr_velocity_obs) {
  return (velocity_raw - median_velocity_obs) / iqr_velocity_obs;
}
```

No cap — velocity is not subject to the burden-magnitude blowup; the level cap
upstream already bounds the trajectory the velocity is differenced from.

The standardized velocity covariate at week `w` is
`standardize_velocity(v(w), median_velocity_obs, iqr_velocity_obs)`.

## 4. Coefficient parameters and hazard application

**No structural change to the coefficient machinery — only its dimension changes
with `n_time_varying_covar`.** The coefficient arrays are already count-driven:

- `time_varying_coef_01/02/03` are sized `[n_time_varying_covar]`
  (`stan/modules/multistate/parameters.stan:43-44, 75-76, 149-150`). In velocity
  mode they become length-2 automatically.
- Hyperpriors (`stan/modules/multistate/hyperparams.stan:114-116`) and priors
  (`stan/modules/multistate/priors.stan:95`) loop over / size by
  `n_time_varying_covar` — they adapt with no edit.
- Hazard application (`stan/modules/multistate/transformed_parameters.stan:175,
  201, 366, 779`) is `log_cond_surv += Σ_k coef[k]·covar[k]` over `k in
  1:n_time_varying_covar` — index-agnostic, works for 2 features unchanged.
- R initializers (`r/initializers_ms.R:247-267`) and prior builders
  (`r/priors.R:43-56`) size off `n_time_varying_covar` — adapt automatically once
  `n_tv_covar = 2L` is set.

**The only places that "know" which feature is which are the two builder files.**
Everything downstream is purely count-driven. The flag's blast radius:

1. `stan/_ms_burden_tv_covar.stan` — branch on `enable_ms_velocity_basis`.
2. `stan/_ms_burden_inline_tv_covar.stan` — mirror the branch (all three blocks).
3. `stan/modules/tumor/transformed_data.stan` — compute the two new observed
   velocity constants.
4. `stan/tumor/_tumor_observed_covar_transformed_data.stan` — alias them.
5. `stan/_burden.stanfunctions` — `standardize_velocity` helper.
6. `stan/modules/multistate/flags.stan` — declare `enable_ms_velocity_basis`.
7. R config wiring (§5).

**Clinical reading.** In velocity mode the 0→1 / 0→3 coupling becomes "log-hazard
responds to where burden is (level) and which way it is moving (velocity)" — a
tumour-kinetics coupling, replacing the decrease/growth-rate-constant coupling.
This requires a one-paragraph update to the model-spec doc
(`quarto/sclc/website/documentation/multistate-specification.qmd`), tracked as
a doc task.

## 5. R config wiring

- Add `enable_ms_velocity_basis` to the multistate flag assembly in the tumour
  config (`targets/sclc_targets.R` flag block ~line 1003, and
  `targets/publication_targets.R`), threaded into the Stan data list like the
  other `enable_ms_*` flags.
- The publication config (`targets/publication_targets.R:230-243` tribble) gains a
  matched switch: `enable_ms_velocity_basis = TRUE` with `n_tv_covar = 2L`
  (velocity mode) **or** `FALSE` with `n_tv_covar = 3L` (legacy). Documented as a
  matched pair; the config-error check (§2) enforces it.
- **Surrogate / warm-start interaction.** The publication tribble carries
  warm-start metric files dimensioned to the parameter space. Switching
  `n_tv_covar` 3→2 drops one `time_varying_coef` element **per enabled
  transition**, changing the sampled-parameter dimension and invalidating any saved
  inv-metric. The velocity-mode run must `warmstart = FALSE` (cold start) for its
  first fit, exactly as the surrogate cold-start note at
  `publication_targets.R:234-242` already establishes for a dimension change.
- Pioneer tribble: set `enable_ms_velocity_basis = FALSE` for all rows (or rely
  on a `FALSE` default in the flag-assembly helper) — no behavioural change to PSA.

## 6. Validation plan (option C: gate now, comparative refit as follow-up)

**Landing gate — blocks the code merge:**

1. **Compiles in both modes.** `stanc --include-paths=stan --include-paths=stan/tumor
   stan/tumor/sf-ssm-log-space.stan` passes with `enable_ms_velocity_basis` both 0
   and 1. (Also syntax-check `sf-ssls-lfo.stan` / `-endpoints.stan` since they
   include the builders.)
2. **Converges + identifies (velocity mode).** A short publication MCMC run in
   velocity mode: Rhat < 1.01, adequate ESS, zero divergences, and
   `time_varying_coef_01` (and `_03` if enabled) is **identified** — posterior
   shifted off the prior, not stuck at the hyperprior mean.
3. **Scale sanity.** The standardized velocity feature is roughly unit-IQR across
   patients and time (confirms `median_velocity_obs` / `iqr_velocity_obs` anchor
   the right scale; a feature with IQR ≫ 1 or ≪ 1 signals a units bug, e.g. a
   missing per-week normalization).

**Documented follow-up — does NOT block landing, has its own acceptance criteria:**

4. **Comparative refit.** Side-by-side publication fits, velocity (2-feature) vs
   legacy (3-feature): compare PFS/OS KM fit quality and coefficient
   interpretability. **Acceptance:** the velocity basis must fit PFS/OS *at least as
   well* as the 3-feature basis before it is declared the publication default.
   Written up as a result in
   `docs/superpowers/specs/2026-06-09-laplace-surrogate-findings-and-blockers.md`
   (or a sibling findings doc).

## 7. Implementation task list (handed to writing-plans)

1. Declare `enable_ms_velocity_basis` in `stan/modules/multistate/flags.stan`.
2. Add `standardize_velocity` to `stan/_burden.stanfunctions`.
3. Compute `median_velocity_obs` / `iqr_velocity_obs` (observed per-week log-SLD
   deltas) in `stan/modules/tumor/transformed_data.stan`.
4. Alias the velocity constants in
   `stan/tumor/_tumor_observed_covar_transformed_data.stan`.
5. Branch `stan/_ms_burden_tv_covar.stan` (dense path) on the flag: velocity via
   central difference of `states_full_grid` (cap-before-difference) when on; legacy
   rates when off.
6. Mirror the branch in `stan/_ms_burden_inline_tv_covar.stan` (all three inline
   blocks): velocity via 3-point analytic `log_sum_exp` central difference when on.
7. Add the config-error check enforcing `(enable_ms_velocity_basis, n_tv_covar)` ∈
   {(1,2), (0,3)} (Stan reject and/or R-side stop).
8. Wire `enable_ms_velocity_basis` into the tumour flag assembly
   (`sclc_targets.R`, `publication_targets.R`) and set the matched
   `(flag, n_tv_covar, warmstart=FALSE)` in the publication tribble; default
   `FALSE` for pioneer.
9. Validation gate: compile both modes; short velocity-mode MCMC for
   convergence/identification/scale (§6 items 1–3).
10. Doc update: one paragraph in the MS model-spec qmd re-justifying `coef` as a
    kinetics (level + velocity) coupling.
11. (Follow-up, non-blocking) comparative refit §6 item 4.

## 8. Open risks

- **Boundary velocity at first/last column.** Forward/backward difference at the
  grid edges is a half-step estimate; if the 0→1 event mass concentrates at week 1
  the first-column velocity is a forward difference, slightly inconsistent with the
  interior central difference. Acceptable (the surrogate's `b₁+2b₂w` is exact at
  every point including the edge, but the forecast path's edge approximation only
  affects the first week). Note it; do not over-engineer.
- **Observed-velocity sample size.** Patients with a single measured visit
  contribute no delta; the `median/iqr` pool is over inter-visit pairs, not
  patients. Should be ample in the publication cohort but worth a guard against an
  empty pool (fall back / reject with a clear message).
- **`coef` sign/scale comparability across modes** is intentionally NOT claimed —
  the bases differ, so legacy and velocity coefficients are not directly
  comparable. The comparative refit (§6.4) compares *fit quality*, not coefficients.
