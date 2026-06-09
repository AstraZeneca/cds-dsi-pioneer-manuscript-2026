# Publication spop-PFS pessimism: root cause = 0→1 baseline-hazard GP shape

**Date:** 2026-06-05
**Branch context:** investigated on `karim/pioneer-pub`; experiment fits in
`karim/pub-burden-only-bridge` (store `publication/burden-only`). Main reference
store: `publication/main`.

## TL;DR

The publication spop PFS KM is systematically pessimistic vs observed, worst in
weeks ~12–36 (amgen wk24: predict 0.48 vs observed 0.64). Root cause is the
**MSM 0→1 baseline hazard being too high/flat in early weeks**: a smooth,
zero-mean, trial-level GP (`enable_ms_level_baseline_hazard = c(trial = 3L, ...)`)
cannot represent the steeply *ramping* empirical hazard, so it sits at a
compromise level that over-predicts early progression. Confirmed it is NOT the
bridge coefficients, NOT informative censoring, NOT a frailty re-draw, NOT the
off-trial-death estimand.

## Evidence chain (each step ruled a candidate in or out)

1. **Burden-only ablation** (drop kinetic rate covariates from the 0→1 bridge):
   spop PFS curve essentially unchanged (≤0.013 survival diff at every week).
   → The B(t)/rate **bridge coefficients are not the driver**. `coef_01[1]`
   flipping sign (+0.215 full → −0.110 burden-only) is a conditional-on-no-target-PD
   quantity (RECIST routing selection), not a bug, and doesn't move the forecast.

2. **Channel decomposition** (amgen, spop): `combined ≈ ms_only` at every week;
   `target_only` (mechanistic RECIST) tracks observed well (wk24: 0.60 vs 0.64).
   → The MSM hazard channel wins the `min()` almost always and is the pessimist;
   the mechanistic channel is well-calibrated.

3. **Off-trial-death estimand check**: 57 `died_off_trial` patients (55 amgen)
   are routed to dropout (state 3) in the likelihood but counted as composite-PFS
   events in `km_trial_pfs`. Re-censoring them to the illness-death estimand closes
   only ~12% of the gap (mean |gap| 0.059 → 0.052) and barely touches wk24.
   → Real estimand subtlety but NOT the pessimism cure. **NOTE:** the spop forecast
   uses the PFS-from-OS graft (dropout death = PFS event at death time), so the
   observed composite KM is actually *consistent in kind* — do NOT censor these in
   the observed KM. (An earlier proposed fix to `prepare_publication_analysis_data`
   was reverted for this reason.)

4. **Sample vs spop split**: conditional `sample_ms` (wk24: 0.68) ≈ observed;
   unconditional `spop_ms` (wk24: 0.48) is the pessimist. Both read the SAME
   per-patient fitted hazard `ms_log_cond_surv[j]` (no frailty re-draw — confirmed
   in `multistate/transformed_parameters.stan:67,131-133`). Difference is only the
   RNG conditioning (sample pins survival to observed `pfs[p]`; spop runs from t=0).

5. **Weekly hazard bins** (amgen) — the decisive table:

   | window  | empirical train | sample_ms | spop_ms |
   |---------|-----------------|-----------|---------|
   | wk 1–12 | 0.0049          | 0.0082    | 0.0115  |
   | wk 13–24| 0.0140          | 0.0241    | 0.0490 (3.5×) |
   | wk 25–36| 0.0390          | 0.0647    | 0.0675  |
   | wk 37–52| 0.0852          | 0.0705    | 0.0705  |

   Empirical hazard **ramps 17×** (0.005→0.085, monotone). Model hazard is
   flat-but-high (0.012→0.049→0.068→0.071). Over-progression concentrated wk13–24.

6. **GP hyperparameters** (fitted, `tumor_ssls_res_posterior`, 0→1 pop GP):
   - `rho`   = 19.5 wk [9, 31]  (prior `inv_gamma(8,135)`, mean 19.3 — barely updated)
   - `alpha` = 1.62
   - `intercept` = −3.92 → baseline ≈ 0.020/wk (≈4× the wk1–12 empirical 0.0049)
   - grid step = 4 wk (50 knots), GP input in **raw weeks** (not normalized).

   → GP is smooth (ρ≈19) and **prior-dominated** (only 2 trials inform it). Being
   **zero-mean around a flat intercept**, it cannot sustain a monotone low→high ramp
   — the prior expects mean-reversion. Result: anchored at a flat ~0.02/wk that is
   too high early, never reaching the late-week 0.085.

## Why this is publication-specific (pioneer comparison)

The MSM + forecast code is **fully shared** (`multistate.stanfunctions`,
`pfs.stanfunctions`, `modules/multistate/*`, `state_space/generated_quantities.stan`,
`_ms_burden_tv_covar`). The divergence is **R-side config**:

| axis | publication | pioneer |
|------|-------------|------------|
| 0→1 baseline | **trial-level GP** (`trial = 3L`) | **intercept-only per arm** (no GP) |
| hazard time-shape from | the GP | the **PSA covariate B(t)** |
| 0→3 dropout | ON (corr. w/ 0→1) | OFF (`enable_ms_03 = FALSE`) |
| propensity borrowing | none | yes (Flatiron RWD) |

**There is no pioneer fix the publication regressed on** — it's a deliberate
architectural choice. pioneer avoids the baseline-GP-shape problem entirely by
keeping the baseline flat and letting the strong PSA kinetic covariate carry the
hazard ramp. Publication's SLD covariate is weaker (burden-only ablation showed it
carries little hazard signal), so it leaned on a GP baseline — which under-identifies
with 2 small trials.

## Fix options (for `karim/pioneer-pub`, must be flag-gated — shared Stan code)

- **(A) GP mean trend** (recommended): change the 0→1 baseline GP mean from a
  constant `intercept` to `intercept + slope·g(t)` (keep the intercept; add ONE
  slope param). `g(t) = log(t)` (Weibull-like ramp, matches the steep-early-then-
  flatten empirical shape) preferred over `g(t)=t` (Gompertz). Lets the monotone
  climb live in `slope` (no mean-reversion penalty); GP models residual wiggle.
  Touches `modules/multistate/{parameters,transformed_parameters,priors}.stan` +
  hyperparams — **gate behind a new flag** so pioneer behavior is unchanged.
- **(B) Shorten lengthscale + fatten left tail**: re-center `rho` prior shorter
  AND lower its `alpha` (e.g. the existing "Relaxed" block α=5) so the GP *can*
  wiggle on ~6–8 wk. Cheap (R prior change only) but partial — a zero-mean GP still
  resists a sustained ramp. Lengthscale is prior-dominated (2 trials), so this is
  "you get the ρ you assume." Likely insufficient alone.
- **(C) pioneer-style intercept-only baseline**: drop the GP, lean on B(t).
  Risky here — SLD's B(t) is a weak hazard driver (per burden-only ablation), so
  the ramp may not materialize. Not recommended as sole fix.

**Recommendation:** (A), optionally with (B)'s prior relaxation riding along.
All require a refit; validate with the same spop-vs-observed PFS KM comparison
(expect the wk12–36 gap to close once early hazard can sit low).

## Implementation sketch for Fix (A) — chosen design

Decisions (2026-06-05):
- **Baseline temporal structure = one unit, separate from covariates.** The trend
  and GP residual are **sibling terms in a single delimited "baseline temporal"
  block**: `log_pop_lambda_01(t) = intercept + slope·g(t) + GP_resid(t)`, computed
  where `log_pop_lambda_01` is built (transformed_parameters ~L55-67), BEFORE and
  separate from the time-varying covariate loop (~L144-145). The trend is NOT a
  covariate and must not go through `ms_time_varying_covar`.
- **0→1 only** for now (the diagnosed transition). Other transitions keep
  intercept+GP unchanged.
- **New per-transition flag**, OFF by default → pioneer and all existing fits
  remain bit-identical; publication opts in.
- `g(t) = log(t)` (Weibull-like; matches steep-early-then-flatten empirical ramp).
  Standardize: `g(t) = log(t) - mean(log(1:max_all_t))` so `slope` is centered and
  the intercept keeps its current meaning.
- **HIERARCHICAL, mirroring the GP exactly.** The baseline GP is hierarchical:
  population GP (`gp_01_pop_*` → `log_pop_lambda_01`) PLUS per-level GP residual
  (`gp_01_level_*` per enabled group → `log_level_lambda_01_residual[g]`), active
  at whatever levels `enable_ms_level_gp[MS_SLOT_01, lv]` turns on (publication:
  trial level). The trend mirrors this:
  - **population slope** `log_lambda_trend_01_pop_slope` → added to `log_pop_lambda_01`.
  - **per-level slope deviation** `log_lambda_trend_01_level` (one per enabled group,
    NCP-pooled by `log_lambda_trend_01_level_sd[lv]`, reusing the SAME group
    infrastructure as `log_lambda_gp_01_level_intercept`:
    `n_enabled_groups_ms_baseline_01`, `raw_level_pos_ms_baseline_slot`,
    scattered via `patient_ms_baseline_flat_idx_slot[MS_SLOT_01, p, lv]`).
  - The per-group slope deviation is folded INTO the existing level loop
    (transformed_parameters L129-147): for each enabled group `g` at level `lv`,
    add `slope_group[g] * g_t` across weeks to `log_level_lambda_01_residual[g]`,
    so it rides the existing patient-scatter (L150-152) for free and is
    automatically hierarchical over whatever levels the GP is enabled at.
  - Trend per-level slope is active ONLY at levels where the baseline GP is active
    (i.e. `enable_ms_level_gp[MS_SLOT_01, lv] == 1`), so the patient level (GP off,
    intercept-only) gets NO trend — matching the GP.

### Stan changes (`modules/multistate/`) — HIERARCHICAL (mirrors GP)

`flags.stan` — add:
```stan
int<lower=0, upper=1> enable_ms_baseline_trend_01;  // 0→1 baseline log-time trend
```

`transformed_data.stan` — add the log-time centering constant:
```stan
real ms_log_t_centering = enable_ms_01 ? mean(log(linspaced_vector(max_all_t, 1, max_all_t))) : 0;
```

`hyperparams.stan` — add (only consumed when flag on):
```stan
real log_lambda_trend_01_pop_mean;          // prior mean, population slope (e.g. 0)
real<lower=0> log_lambda_trend_01_pop_sd;   // prior sd,  population slope (e.g. 0.5)
array[n_levels] real<lower=0> log_lambda_trend_01_level_sd_sd;  // half-normal scale for per-level slope SD (e.g. 0.25)
```

`parameters.stan` — add (size 0 when flag off → bit-identical when off):
```stan
// Population slope
array[enable_ms_baseline_trend_01 ? 1 : 0] real log_lambda_trend_01_pop_slope;
// Per-group slope deviation (NCP), reusing the baseline-01 enabled-group count;
// per-level SD for pooling. Both size 0 when flag off.
vector[enable_ms_baseline_trend_01 ? n_raw_groups_ms_baseline_01 : 0] raw_log_lambda_trend_01_level;
array[enable_ms_baseline_trend_01 && any_re_level_slot[MS_SLOT_01] ? n_levels : 0] real<lower=0> log_lambda_trend_01_level_sd;
```

`transformed_parameters.stan` — TWO insertions inside `if (enable_ms_01) { ... }`:

(1) **population slope** — sibling to the pop GP in the baseline temporal block
(this part is ALREADY written, L75-82) but generalize the param name to
`log_lambda_trend_01_pop_slope[1]`.

(2) **per-level slope deviation** — fold into the existing level loop (L89-154),
mirroring how `log_lambda_gp_01_level_intercept[g]` is built and used. For each
enabled group `g` at level `lv` WHERE the trend is active (gate on
`enable_ms_baseline_trend_01 && enable_ms_level_gp[MS_SLOT_01, lv]`):
```stan
// scale the per-group raw slope by the per-level SD (NCP), same as the GP intercept
real slope_g = raw_log_lambda_trend_01_level[<raw idx for g>]
               * log_lambda_trend_01_level_sd[lv];
// add the trend deviation across weeks to this group's residual
for (t in 1:max_all_t)
  log_level_lambda_01_residual[g, t] += slope_g * (log(t) - ms_log_t_centering);
```
The existing patient-scatter at L150-152 then propagates it for free.

`priors.stan` — add (guarded):
```stan
if (enable_ms_baseline_trend_01) {
  log_lambda_trend_01_pop_slope[1] ~ normal(log_lambda_trend_01_pop_mean, log_lambda_trend_01_pop_sd);
  raw_log_lambda_trend_01_level ~ std_normal();  // NCP
  if (any_re_level_slot[MS_SLOT_01])
    for (lv in 1:n_levels)
      log_lambda_trend_01_level_sd[lv] ~ normal(0, log_lambda_trend_01_level_sd_sd[lv]);  // half-normal (lower=0)
}
```

Note: reuse the existing `raw_level_pos_ms_baseline_slot[MS_SLOT_01]` /
`enabled_level_pos_ms_baseline_slot[MS_SLOT_01]` position arrays to index the
per-group raw slopes exactly as the GP intercept does — do NOT introduce a new
group-counting scheme.

### R changes

- `targets/publication_targets.R` — set `enable_ms_baseline_trend_01 = 1L` in the
  publication stan-data assembly (near the `enable_ms_level_baseline_hazard` config,
  ~L249). pioneer/sclc leave it unset → defaults 0.
- `r/priors.R` (or `get_tumor_priors`) — add defaults
  `log_lambda_trend_01_pop_mean = 0`, `log_lambda_trend_01_pop_sd = 0.5`,
  `log_lambda_trend_01_level_sd_sd = rep(0.25, n_levels)`.
- `r/initializers*.R` — when flag on, init `log_lambda_trend_01_pop_slope` (small
  positive, ~0.1), `raw_log_lambda_trend_01_level` (zeros, length
  `n_raw_groups_ms_baseline_01`), `log_lambda_trend_01_level_sd` (small, ~0.1);
  size-0 / omit when off.
- Default the flag to `0L` wherever stan data is assembled so non-publication
  models are untouched (check `_full_model_data.stan` / assembly defaults).

### Validation

Refit `tumor_ssls_res_posterior_*` into a fresh store; rebuild
`tumor_ssls_km_rvar_posterior_*`; rerun the spop-vs-observed PFS KM table and the
weekly-hazard-bin check. Success = wk1–24 hazard drops toward empirical
(0.005/0.014) and wk12–36 spop PFS gap closes. Also confirm `slope` posterior is
positive and well-identified (>0 with interval excluding 0).

### Guards / gotchas

- Size-0 parameter arrays when flag off → bit-identical to current for pioneer.
- Trend must be added to `log_pop_lambda_01` ONLY (population baseline), not to
  level residuals or covariates — keeps the "baseline temporal unit" clean.
- `calc_gp_pred` mean stays the constant `intercept`; the trend is a sibling
  addition (per chosen design B), NOT folded into the GP mean argument.
- Watch identifiability between `slope` and the GP `eta` (both shape the trajectory).
  The GP residual should shrink toward zero once the trend absorbs the ramp; if the
  GP fights the trend, consider tightening the GP `alpha` prior alongside.

## Secondary issue found (fix regardless)

The `rho` prior `inv_gamma(8,135)` (mean 19.3 wk) is fine in *value* but the GP
input is in **raw weeks**, not normalized — so any prior re-use across trials with
different `max_all_t` is fragile. Consider normalizing the GP input or documenting
the week-scale assumption. Low priority; not the pessimism driver.

## Key files

- `targets/publication_targets.R:249` — `enable_ms_level_baseline_hazard = c(trial = 3L, patient = 0L)`
- `r/multistate.R:480` — `decompose_ms_level_baseline_hazard()` (mode 3 → GP + RE-NCP intercept)
- `r/priors.R:64-74` — 0→1 GP pop hyperpriors (intercept/alpha/rho)
- `r/multistate.R:446` — `ms_gp_grid_step = 4L`
- `stan/modules/multistate/transformed_parameters.stan:55-135` — baseline GP build + level scatter
- `stan/modules/multistate/priors.stan:14` — `rho ~ inv_gamma(...)`
- `stan/pfs.stanfunctions:1956,2004-2009` — spop `min()` combine + PFS-from-OS dropout graft
- `r/publication/prepare_analysis_data.R` — off-trial policy (do NOT censor; PFS-from-OS handles it)
