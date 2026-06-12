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
κ > 0 produces concave-down deceleration. That curvature is temporally distinct from the
decrease rate (which governs the early declining arm), so κ is not aliased by any existing
parameter. This is the live path; the static machinery stays committed but OFF.

## Goal

Replace the unbounded linear-in-time growth accumulation with Gompertz decay: the growth
rate attenuates as `growth_rate · e^{−κ_i·t}`, where κ_i is a **per-patient decay rate**
driven by a linear predictor (population intercept now; baseline-covariate slopes and
higher hierarchy levels available later). The whole mechanism is gated by a single flag
`enable_gr_decay`, default OFF, which recovers today's model **exactly**.

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
(no `exp`, no warp, zero added cost) — bit-identical to the current model. This mirrors the
`negative_infinity()` backward-compat trick used for the static compartment: the off-path is
a true no-op, not a κ value that approximates linearity.

### New module: `stan/modules/gr_decay/`

κ gets the **same N-level linear-model scaffolding** as `tr` and `frac` (it is modeled on the
`frac` module, which is the closest analog — both shape the growth process). Per the
"Adding Module Parameters" guidelines in `.claude/rules/stan-guidelines.md`, the module
contains the standard files:

| File | Contents (mirrors `frac`, `exp` link instead of `inv_logit`) |
|---|---|
| `flags.stan` | `enable_gr_decay`; `enable_pop_cov_gr_decay`; `enable_level_intercept_gr_decay[n_levels]`; `enable_level_cov_gr_decay[n_levels]` |
| `hyperparams.stan` | `gr_decay_log_loc_pop_mean`/`_sd`; per-level intercept/slope SD hyperpriors; QR-coef hyperparams — all guarded so they cost nothing when the flag is off |
| `parameters.stan` | `gr_decay_log_loc_pop`; QR-space `gr_decay_coef_qr_pop` (len 0 if disabled); the unified level intercept/slope SD + raw/cp draws, **all sized `… ? … : 0` on `enable_gr_decay`** so an off model samples nothing extra |
| `transformed_parameters.stan` | assembles `gr_decay_log_loc_patient` (pop intercept + QR pop covariate effects + level intercepts + level slopes), then `gr_decay_patient = exp(gr_decay_log_loc_patient)` — a length-`n_forecast_patients` vector of κ_i (length 0 / unused when flag off) |
| `priors.stan` | population-intercept prior + (guarded) covariate/level priors, identical loop structure to `frac/priors.stan` |
| `transformed_data.stan` | the flat-index precomputation parallel to `frac` (only the structures κ's hierarchy needs) |
| `generated_quantities.stan` | expose `gr_decay_kappa_pop = exp(gr_decay_log_loc_pop)` and the implied plateau offset `growth_rate_pop / kappa_pop` for interpretability |

**Link function:** κ is modeled on the **log scale** (`gr_decay_log_loc_*`) so κ_i = exp(linpred) > 0 is guaranteed for every patient. This matches `tr` (log-scale rate) rather than `frac` (logit-scale fraction); the hierarchy *plumbing* is copied from `frac`, the *link* from `tr`.

**Initial hierarchy configuration (this branch):** intercept-only, population level only.
`enable_pop_cov_gr_decay = 0`, `enable_level_intercept_gr_decay = rep(LEVEL_MODE_NONE, n_levels)`,
`enable_level_cov_gr_decay = rep(0, n_levels)`. So `gr_decay_log_loc_patient` reduces to the
scalar `gr_decay_log_loc_pop` broadcast across patients. The covariate/level machinery is built
and compiled but inert — turning on baseline-covariate response later is a flag flip plus
prior, no structural change (exactly how `tr`/`frac` already work).

### Wiring κ into the growth accumulation

The live forward simulation lives in
`stan/modules/state_space/transformed_parameters.stan`. There are **three** growth-accumulation
branches; each warps time differently because of how the per-step rate behaves:

1. **No-noise, direct visit computation** (≈ line 184):
   `state_g = init_g + growth_rate · dt`
   → `state_g = init_g + growth_rate · φ_i(dt)` (warp the elapsed time per visit).

2. **No-noise / no-grid vectorized full grid** (≈ lines 74, 119, 145):
   `state_g = init_g + growth_rate · t`
   → `state_g = init_g + growth_rate · φ_i(t)` (warp the grid time axis).

3. **Process-noise cumsum branches** (≈ lines 74, 119; pop and/or patient noise on):
   here the per-step rate is **not** constant (it carries the AR-style deviations), so we
   **cannot** simply warp cumulative time. Instead, multiply each increment by the decay
   weight `e^{−κ_i·t_k}` evaluated at that step's time **before** the cumulative sum:
   `growth_inc_det[k] = growth_rate · e^{−κ_i·t_k} · Δt_k`
   then `state_g = init_g + cumsum(growth_inc_with_noise)`. As κ_i → 0 the weights → 1 and
   this reduces to the current cumsum exactly. The decay weight composes cleanly with the
   existing noise increments already in those vectors.

All three reduce to the present code when `enable_gr_decay = 0`. The spec deliberately does
**not** push κ down into the `sf.stanfunctions` trajectory helpers (`calc_states`,
`sf_log_space_trajectory_ncp`) — those are not on the live path (only `legacy/` and the LFO
file call them). The change is confined to the live `transformed_parameters.stan` plus the LFO
forecast path, mirroring how the static work kept the surgery contained.

#### Numerical evaluation of φ and the decay weight

Evaluate `φ_i(t) = (1 − e^{−κt})/κ` with a guard for small κ to avoid 0/0:
- use `−expm1(−κ·t)/κ` (numerically stable for the numerator), and
- when `κ·t` is below a small threshold, fall back to the series `t·(1 − κt/2 + …)` or simply
  `t` (the flag-off path already covers exact κ = 0; the guard only matters for tiny but
  nonzero κ inside the flag-on path).

### LFO forecast path

`stan/tumor/sf-ssls-lfo.stan` and `stan/tumor/_lfo_endpoints_generated_quantities.stan`
generate long-horizon forecasts via the trajectory functions. The Gompertz weight must apply
there too (this is where "exploding SLD at long horizons" is most visible). The plan threads
κ_i (or `gr_decay_log_loc_patient` mapped through the unified `forecast_patient_idx`) into the
forecast accumulation exactly as the static `static_log_level_per_patient` vector was threaded
in the prior work — same indexing discipline (forecast-local `j` → unified
`forecast_patient_idx[j]`), same gating.

## Data flow

```
baseline covariates (design matrix, QR) ──┐
                                          ▼
gr_decay_log_loc_pop  +  (QR pop cov)  +  (level intercepts)  +  (level slopes)
                                          │   [all inert this branch except the pop intercept]
                                          ▼
                          gr_decay_log_loc_patient  ──exp──▶  κ_i  (per patient, > 0)
                                          │
            ┌─────────────────────────────┼──────────────────────────────┐
            ▼                              ▼                              ▼
   no-noise visit branch        vectorized full-grid branch     process-noise cumsum branch
   growth_rate·φ_i(dt)          growth_rate·φ_i(t)              cumsum(growth_rate·e^{−κ_i t_k}·Δt_k + noise)
            └─────────────────────────────┼──────────────────────────────┘
                                          ▼
                          state_growth(t)  →  log_sum_exp with state_decrease  →  log-SLD mean
```

## Identifiability — the gate

κ is the whole bet, and it is the thing that could repeat the static failure. The gate is a
**recoverability simulation** (same discipline as `r/process_noise/recoverability_sim.R`):

1. Simulate N patients (e.g. 300) from the generative form with a **known** population κ
   (e.g. κ_true giving a visible plateau within ~3× the observed follow-up), realistic visit
   schedule, and measurement noise.
2. Assemble stan-data with `enable_gr_decay = 1`, priors from `r/priors.R`.
3. Fit `stan/tumor/sf-ssm-log-space.stan` short (2 chains, 500/500).
4. Extract `gr_decay_kappa_pop`, compare posterior to truth, report divergences.

**Pass:** posterior 90% CI for κ_pop covers truth, with curvature visibly recovered.
**Fail (κ non-identifiable):** posterior collapses to 0 or its CI misses truth — then the
honest read is that in-window data is agnostic to curvature, and the mechanism is only useful
as a **prior-driven forecasting constraint** (a weakly-informative prior nudging κ slightly
positive so extrapolation plateaus regardless). The spec anticipates this: even a "fail" here
is informative and points to the prior-as-regularizer fallback rather than abandonment, which
is a materially better outcome than static (where the parameter was structurally aliased, not
merely weakly informed).

### Prior choice

`r/priors.R` gets `gr_decay_log_loc_pop_mean` / `_sd` defaulting to a weakly-informative prior
centered on a *mild* decay — small positive κ corresponding to a plateau several years out, so
that (a) in-window fit is barely perturbed and (b) long-horizon forecasts are bounded. The
exact center is set so the implied half-life of the growth rate is on the order of the longest
forecast horizon of interest. Initializer (`r/initializers.R`) draws the length-1 pop param
when the flag is on, length-0 when off (mirroring the static param handling).

## Error handling / edge cases

- **κ = 0 (flag off):** literal linear path, bit-identical to current model. The backward-compat
  gate is the full existing testthat suite staying green with the flag defaulted OFF.
- **Tiny κ (flag on):** `expm1`-based φ + small-κ guard prevents 0/0; reduces to `t`.
- **Process noise on:** decay applied at the increment level, composes with existing AR noise.
- **Patient with single growth visit:** φ_i(dt) with one interval is well-defined; no special case.
- **Sized-zero params when off:** every new parameter/data block is guarded so an off model
  samples and stores nothing extra (no wasted sampling, no posterior bloat).

## Testing

- **Stan unit (warped time):** a tiny test model exercising `φ_i(t)` confirms φ(t)→t as κ→0
  and φ(t)→1/κ as t→∞, plus the small-κ guard (mirrors `test-stan-calc-log-burden-mean.R`).
- **R unit (hierarchy reduces correctly):** with the pop-intercept-only config,
  `gr_decay_log_loc_patient` is constant across patients = `gr_decay_log_loc_pop`.
- **Backward-compat (the headline gate):** `enable_gr_decay = 0` → full testthat suite green;
  CmdStanR compiles `sf-ssm-log-space.stan` and `sf-ssls-lfo.stan` with no errors.
- **Recoverability sim (viability gate):** as above — decides whether κ is data-identifiable
  or prior-driven before any real fit.

## Out of scope (this branch)

- Baseline-covariate response on κ (machinery built, flag OFF).
- Higher hierarchy levels on κ (trial/arm/patient REs — machinery built, modes NONE).
- Reviving the legacy `growth_lag` start-delay.
- Stochastic (AR1/OU) time-varying κ. Gompertz is deterministic decay; a stochastic κ is a
  possible later step but is explicitly deferred.
- Logistic/Verhulst carrying-capacity or Simeoni exponential→linear alternatives (considered,
  not chosen — Gompertz is the lowest-surgery, most-identifiable fit to this model's existing
  log-linear growth accumulation).

## Files touched (summary)

**New:** `stan/modules/gr_decay/{flags,hyperparams,parameters,priors,transformed_data,transformed_parameters,generated_quantities}.stan`

**Modified:**
- `stan/tumor/sf-ssm-log-space.stan` — `#include` the new module's fragments (alongside `tr`/`frac`/`init`)
- `stan/modules/state_space/transformed_parameters.stan` — warp the three growth-accumulation branches
- `stan/tumor/sf-ssls-lfo.stan`, `stan/tumor/_lfo_endpoints_generated_quantities.stan` — apply decay on the forecast path
- `r/priors.R` — `gr_decay_log_loc_pop_mean`/`_sd` + (guarded) covariate/level hyperparams
- `r/initializers.R` — conditional length-1 pop param init
- `targets/publication_targets.R` (and other project target files as needed) — `enable_gr_decay = 0L` default + the pop-only hierarchy config
- `r/process_noise/recoverability_sim.R` — extend / add a κ-recoverability variant
- model-spec documentation page — document the Gompertz growth-rate decay

**New tests:** Stan warped-time test + R hierarchy-reduction test (+ the recoverability sim is the viability gate).
