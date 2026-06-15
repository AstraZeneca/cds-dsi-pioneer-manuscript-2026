# HANDOFF — Write the spec for the (level, velocity) MS coupling-basis change

**For:** a new session. **Branch:** `karim/laplace`. **Date:** 2026-06-09.
**Your task:** brainstorm + write a design spec (NOT implementation) for replacing
the multistate burden-coupling covariates with a **(level, velocity)
trajectory-derivative basis**, model-wide. This is "§9 step 0" — the gated
prerequisite for the Laplace joint-surrogate work, but it is **independently
useful to the publication model** regardless of the surrogate.

Start by invoking `superpowers:brainstorming`. Do NOT write code. Produce a spec
in `docs/superpowers/specs/YYYY-MM-DD-level-velocity-coupling-design.md`, get it
reviewed, and stop before implementation.

---

## Why this change exists (the motivation)

The MS hazards for the burden-coupled transitions currently read **three**
time-varying covariate features (`n_time_varying_covar = 3` in the publication):

1. standardized log-burden (the "level")
2. log-**decrease**-rate — a bi-exponential state-space component parameter
3. log-**growth**-rate — the other bi-exponential component parameter

Features 2–3 are properties of the *two-component bi-exponential SSM dynamics*.
They are problematic for two reasons:

- **They have no analog under the Laplace surrogate** that marginalizes
  background-trial patients' tumor latents (the surrogate represents burden as a
  quadratic `g(w)=b0+b1·w+b2·w²`, which has no `dec_rate`/`gro_rate` to read off).
- With `enable_patient_process_noise_tr = FALSE` (publication), those rates are
  **constant per patient** anyway (`rep_row_vector(...)`) — static baseline-rate
  covariates wearing a time-varying costume.

**The replacement:** use the **derivatives of the (log-)burden trajectory**:

| Feature | Meaning | Forecast computation | Surrogate (background) computation |
|---|---|---|---|
| level | `g(w)` log-burden | existing `states_full_grid` log_sum_exp | quadratic `b0+b1·w+b2·w²` |
| velocity | `g'(w)` rate of change | numerical derivative across `states_full_grid` columns | analytic `b1+2·b2·w` |

**Why this is the right basis:**
- velocity is **linear in the marginalized latents `(b1,b2)`** ⇒ feeding it into
  the log-hazard keeps the survival term log-concave (the property the joint
  surrogate gate validated). The bi-exponential rates do NOT have this property
  under the surrogate.
- velocity is computable **identically in form** for forecast and background
  patients — dissolving the forecast/background featurization asymmetry.
- velocity has a clean clinical reading (tumor kinetics: "where you are + which
  way you're moving"), akin to tumor growth-rate constants in the oncology lit.

**Avoid** features nonlinear in the latents — time-to-nadir, nadir depth,
relative velocity `g'/g`, doubling time — they break the surrogate's
log-concavity.

---

## Scope & blast radius (CRITICAL — this is cross-project)

This is **NOT a surrogate-only change.** It alters the core MS model and affects
**forecast patients too**, and the shared feature-builder file is `#include`d
across projects:

- **`stan/_ms_burden_tv_covar.stan`** — the dense-grid feature builder (the file
  that computes features 1–3 at `:64-80`). Included by:
  - `stan/psa/pioneer.stan` (PSA / Pioneer project)
  - `stan/tumor/sf-ssm-log-space.stan` (the publication full model)
  - `stan/tumor/sf-ssls-lfo.stan`, `sf-ssls-lfo-endpoints.stan`
  - `stan/tumor/_tumor_observed_covar_transformed_data.stan`
  - `stan/modules/multistate/transformed_parameters.stan`, `cond_surv_transform.stan`
- **`stan/_ms_burden_inline_tv_covar.stan`** — the inline (no-full-grid) path,
  same feature semantics; mirror any change here.
- **Pioneer uses BOTH `n_tv_covar = 2L` and `3L`** across its tribble variants
  (`targets/pioneer_targets.R:427` tribble, `:542` consumes `n_tv_covar`).
  So changing what "feature 2/3" *mean* touches pioneer's fitted models — the
  spec MUST decide pioneer's disposition (migrate it too? keep a legacy flag?
  PSA burden has its own dynamics — does "velocity" even mean the same thing for
  PSA?).
- **Publication** sets `n_tv_covar = 3L` (`targets/publication_targets.R:233`),
  consumed at `:261`.

**The spec must explicitly scope:** publication-only first, or all projects at
once? A legacy flag to keep the 3-feature basis available, or a hard replacement?
(Repo convention from CLAUDE.md: *no backward-compat aliases unless requested* —
so lean toward hard replacement, but the cross-project fits make this a real
decision to surface to the user.)

---

## Key implementation anchors (for the spec, not to edit yet)

- **Feature builder (dense):** `stan/_ms_burden_tv_covar.stan:52-82` — feature 1
  at `:64` (standardized log-burden), features 2/3 at `:66-82` (the rates to
  replace). `states_full_grid` carries the log-burden trajectory; **velocity =
  numerical derivative across its columns** (e.g. forward/central difference over
  the weekly grid).
- **Feature builder (inline):** `stan/_ms_burden_inline_tv_covar.stan:41-123`.
- **Standardize fn:** `standardize_log_burden(...)` at
  `stan/_burden.stanfunctions:34-40` (caps at `fmin(·,10)`, then
  `(log_burden_abs − median)/iqr`). **Velocity needs its OWN standardization** —
  log-burden's `median_log_burden_obs`/`iqr_log_burden_obs` do not apply to a
  velocity (different units/distribution). Add `median_velocity_obs` /
  `iqr_velocity_obs` (or analogous) computed from observed-trajectory numerical
  derivatives, on the forecast/data path. This is a NEW transformed-data constant
  pair.
- **Coefficient params:** `time_varying_coef_01` / `_03` are sized by
  `n_time_varying_covar` (`stan/modules/multistate/parameters.stan:43-44`;
  hyperpriors at `hyperparams.stan:114-116`; priors at `priors.stan:95`). Dropping
  3→2 features changes these array dimensions.
- **Hazard application:** `stan/modules/multistate/transformed_parameters.stan:176`
  (continuous) and `:202,206` (visit-gated) — `log_cond_surv += coef[k]·covar[k]`.
- **R config:** `n_tv_covar` in `targets/publication_targets.R:233` (and the
  pioneer tribble at `:427`/`:542`).
- **Standardization constants are computed in STAN transformed-data, not R.**
  `median_log_burden_obs`/`iqr_log_burden_obs` are aliases set per marker:
  `stan/tumor/_tumor_observed_covar_transformed_data.stan:48`
  (`= median_log_sld_obs`) and `stan/psa/_psa_observed_covar_transformed_data.stan:22`
  (`= median_log_psa_obs`); the upstream `median_log_*_obs` are built in
  `stan/_observed_covar_transformed_data.stan` (and the per-marker files) from the
  observed data passed in. **The new `median_velocity_obs`/`iqr_velocity_obs` must
  be computed there** (from observed-trajectory numerical derivatives) — likely
  the velocity samples need to be passed in as data from R, so the spec must trace
  what observed-covariate data the R targets currently assemble and add the
  velocity inputs. Grep `median_log_sld_obs` / `median_log_psa_obs` to find the
  build sites and their data dependencies.

---

## What the spec must decide (open questions to brainstorm)

1. **Scope:** publication-only vs all-projects (esp. pioneer/PSA). Legacy flag
   or hard replacement?
2. **Velocity definition precisely:** of *log*-burden or *absolute* burden? Which
   numerical-derivative scheme (forward/central/over what week spacing)? Units?
3. **Velocity standardization:** how/where to compute `median/iqr` of observed
   velocities (from raw assessment deltas? from fitted trajectories?). Must match
   between forecast feature and surrogate analytic velocity (the contract's gate
   adds a feature-SCALE comparison — see below).
4. **Clinical re-justification:** `coef`'s interpretation changes from
   "decrease/growth-rate coupling" to "kinetics (level + velocity) coupling."
   Needs a model-spec doc update (`quarto/.../multistate-specification.qmd`).
5. **Validation plan:** a **refit comparison** vs the current 3-feature model
   (does the kinetics basis fit PFS/OS as well or better? do the coefs identify?).
6. **PSA caveat:** does "velocity" transfer to the PSA burden marker, or does
   pioneer need a different treatment?

---

## How this connects to the surrogate work (context, do not action)

This change is **step 0** of the dependency chain in the joint-surrogate index
contract (`docs/superpowers/specs/2026-06-09-joint-surrogate-index-contract.md`,
§9). The surrogate's functor **assumes the (level, velocity) basis as input** and
cannot be validated until this lands. The contract's §6 O1 and §9 step 2 also
require a **forecast-vs-surrogate feature-SCALE check** (the forecast velocity is
a numerical derivative of a capped/bi-exponential/noisy trajectory; the surrogate
velocity is the uncapped analytic `b1+2b2w` — they are consistent in *form*, not
bit-identical, and `tv_coef` is shared model-wide). Read that contract's O1 and §9
for the downstream constraints your spec should be compatible with — but this spec
stands on its own and should be written for the model's own sake.

## Repo conventions the new session must follow
- Read `CLAUDE.md`, `.claude/rules/stan-guidelines.md`, `.claude/rules/r-guidelines.md`.
- Stan syntax check: `~/.cmdstan/cmdstan-2.39.0/bin/stanc --include-paths=stan
  --include-paths=stan/tumor stan/tumor/sf-ssm-log-space.stan` (NOTE: cmdstan
  **2.39** is pinned on this branch, not 2.38 — see `.Rprofile`).
- tidyverse + native pipe `|>`; no `install.packages`; no hardcoded subject IDs.
- Memory note: `[[laplace-velocity-coupling-basis]]` in the user's auto-memory
  captures this decision.
