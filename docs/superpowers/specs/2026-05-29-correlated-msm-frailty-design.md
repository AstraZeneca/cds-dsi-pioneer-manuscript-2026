# Decomposed MSM Level Effects with Correlated Frailty Intercepts

**Date:** 2026-05-29
**Status:** Approved — ready for implementation planning
**Branch:** `karim/pioneer-pub`

## Problem

The superpopulation (spop) PFS Kaplan–Meier curve is systematically **pessimistic**
relative to the observed KM — it predicts progression substantially earlier than
patients actually experience. Diagnostics on the state-3-covariate fit (run
`202605282352`) localize the bias to one cohort:

| ms_pattern | n | obs median PFS | spop median PFS | gap (wk) | P(0→1) | P(dropout 0→3) |
|---|---|---|---|---|---|---|
| **died_off_trial** | 57 | **42** | 22 | **−20.5** | **77%** | **17%** |
| progressed_alive | 57 | 32 | 23 | −7 | 82% | 13% |
| progressed_died | 298 | 25 | 20 | −4 | 83% | 12% |
| admin_censored | 53 | 21 | 22 | −1 | 84% | 11% |
| died_on_trial | 32 | 17 | 19 | +2.5 | 79% | 12% |

The `died_off_trial` patients followed the 0→3→2 path (dropout, then off-trial death
much later) and have the **longest** real PFS. The model collapses them to a fast
0→1 progression at week 22 instead of routing them to dropout.

### Root cause

Per-patient dropout probability is **flat (~11–17%) across every cohort** — the model
assigns essentially the same 0→3 hazard to everyone. Inspection of the fitted
coefficients confirms why:

- **0→3 covariate coefficients are statistically inert.** Every
  `time_invariant_coef_qr_03[1..6]` straddles zero; `time_varying_coef_03[1..3]` have
  posterior SDs (0.13–0.38) barely shrunk from the `Normal(0, 0.5)` prior. With only
  **57 dropout events** there is not enough signal to identify 9 covariate coefficients
  against the 0→3 baseline hazard.
- **0→1 covariate coefficients ARE identified** (e.g. `time_invariant_coef_qr_01[5]`
  = +0.26 [+0.14, +0.38]; SDs ~0.07, half the 03 SDs), because 0→1 has hundreds of
  events.

The 0→1 and 0→3 transitions use the **same design matrix** (`age + male + ecog + hgb
+ ldh_log + albumin` plus latent log-SLD / decrease / growth). Sharing the regressors
does not share the *information*: in a cause-specific competing-risks model each
transition's coefficients are identified only from that cause's events.

Critically, **there is currently no patient-level random effect on any hazard**
(`enable_ms_level_baseline_hazard = c(trial = 3L, patient = 0L)`). A patient who is
censored for 0→1 for a long time (because they dropped out and progressed slowly) has
**no parameter to absorb that "slow-progressor" signal**.

## Goal

Let the abundant 0→1 survival history inform each patient's dropout propensity, so
`died_off_trial` patients route 0→3 (long PFS) in the spop re-roll instead of 0→1
(short PFS) — without manufacturing covariate signal that the 57 events cannot support.

## Solution overview

A **negatively-correlated patient-level frailty** on the 0→1 and 0→3 hazards. A patient
who is intrinsically slow on 0→1 is, through the correlation, dropout-prone on 0→3.
The patient-level 0→1 intercept is well-fed by all 497 patients' progression/censoring
histories; the negative correlation transmits that evidence to 0→3 without relying on
the sparse dropout events to identify covariate coefficients.

Rather than bolt on a special-purpose frailty, we **generalize the existing level-effect
machinery**. Today a single per-level flag `enable_ms_level_baseline_hazard[lv]` (0–4)
conflates three independent decisions and is shared across all transitions. We decompose
it and add cross-transition correlation as the natural extension of the existing
independent-intercept design (correlation = the `Ω ≠ I` generalization).

## Design

### Section 1 — Configuration schema

Replace the single shared flag with **three orthogonal, per-transition × per-level
arrays**. All transitions are treated identically; nothing is hardcoded per transition.

```stan
// 1. Hierarchical baseline GP residual per transition × level
array[N_TRANS, n_levels] int<lower=0,upper=1> enable_ms_level_gp;

// 2. Level intercept mode per transition × level
//    0 = none, 1 = FE, 2 = RE (NCP), 3 = RE (CP)
array[N_TRANS, n_levels] int<lower=0,upper=3> ms_level_intercept_mode;

// 3. Correlation grouping (Approach B): same positive code at a level
//    => shared MVN/LKJ block; 0 = independent singleton (today's behavior)
array[N_TRANS, n_levels] int<lower=0> ms_level_intercept_corr_group;
```

`N_TRANS` indexes the **intercept slots**, one per additive hazard-intercept channel.
Because the 1→2 transition can carry two GP intercepts (`12_s` sojourn-clock and
`12_t` clock-forward, used together in the extended time-scale mode), the honest slot
set has **six entries**: `{01, 02, 03, 12_s, 12_t, 32}`. The correlation operates on
individual scalar intercepts, so each slot is a potential MVN row. (For this analysis
only `{01, 03}` are correlated, so the `12_*` distinction is moot here, but the schema
must enumerate slots precisely.)

This decomposes the three decisions the old flag conflated:

| Decision | Old (mode 0–4, shared) | New |
|---|---|---|
| Hierarchical baseline GP | bundled in mode 3 | `enable_ms_level_gp[k, lv]` |
| Level intercept + FE/RE/CP | bundled in mode 1/2/4 | `ms_level_intercept_mode[k, lv]` |
| Correlation across transitions | did not exist | `ms_level_intercept_corr_group[k, lv]` |

**Validation (`fatal_error` in `transformed_data`):**
- Any transition with `corr_group > 0` at a level must have `intercept_mode ∈ {2,3}` (RE).
- All transitions sharing a `corr_group` code at a level must share the same group
  partition at that level (equal enabled group counts). This is the real
  well-definedness check — fully general, **not** a hardcoded transition exclusion.
- A `corr_group` with a single member collapses to an independent singleton.

**Level-generic.** The machinery is indexed by generic `(transition, level)`. "Patient"
and "trial" are just levels whose groups are patients / trials. Any level can host
correlated intercepts. Our analysis *configures* the patient level; the code never
assumes which level that is.

**Configuration for this analysis** (not baked into the code): at the patient level,
`intercept_mode[01]=2`, `intercept_mode[03]=2`, `intercept_mode[02]=none`;
`corr_group[01,patient]=1`, `corr_group[03,patient]=1`. Trial level unchanged
(`02` and the sojourn transitions left at `corr_group=0`).

### Section 2 — Stan parameters & transformed parameters

`transformed_data.stan` group-sizing generalizes from one shared mode vector to
per-transition mode vectors: `enable_ms_level_gp` sizes the GP `eta`/intercept blocks,
`ms_level_intercept_mode` sizes the intercept raw/cp buckets
(`compute_ms_transition_group_counts` / `split_cp_ncp_pos` take per-transition modes).

For each level `lv` and each correlation group `g` with `d ≥ 2` members:

```stan
// parameters.stan
cholesky_factor_corr[d]                      L_ms_intercept_corr[lv, g];   // LKJ Cholesky
matrix[d, n_forecast_groups_per_level[lv]]   z_ms_intercept[lv, g];        // NCP std-normal

// transformed_parameters.stan — one row per member transition, one col per GROUP at level lv
matrix[d, n_grp_lv] u =
  diag_pre_multiply(sigma_members, L_ms_intercept_corr[lv,g]) * z_ms_intercept[lv,g];
// scatter row m -> member transition's level-lv intercept vector, gathered to
// forecast units via the existing patient_ms_baseline_flat_idx[:, lv]
```

- `sigma_members` = the existing per-transition level-`lv` intercept SDs
  (`log_lambda_gp_0k_level_intercept_sd[lv]`), unchanged — **no new scale hyperparameters.**
- Singletons / `corr_group=0` / `d=1` keep the **exact current scalar `σ·raw` NCP path**.
  When no correlation is configured the generated parameter vector and density are
  **bit-identical to today.**
- **Addition site is unchanged.** Each transition's per-patient intercept already adds
  as a flat shift `log_cond_surv_0k += rep_matrix(intercept_term, max_all_t)` (e.g.
  `transformed_parameters.stan:730` for 0→3, `:185` for 0→1). The frailty rides the
  intercept channel that already exists; only its *generation* changes
  (correlated MVN vs independent scalar).

NCP is mandatory: sample standard-normal `z` and LKJ-Cholesky `L`, form
`u = diag_pre_multiply(σ, L) z`. This avoids the centered funnel where `u_i` and `σ`
couple — critical with only 57 dropout events feeding `σ_03`.

### Section 3 — Priors & identification

The only new prior is the LKJ on each correlation block; σ scales reuse existing priors.

```stan
// priors.stan — one per configured (lv, g) block
L_ms_intercept_corr[lv,g] ~ lkj_corr_cholesky(2);   // eta = 2
to_vector(z_ms_intercept[lv,g]) ~ std_normal();
```

**`eta = 2`** gently concentrates toward the identity (ρ near 0) without forbidding
strong correlation, and is **symmetric about zero** — it does not bake in the expected
negative sign; the data reveal it. Conservative when one member is data-poor.

**Constraint: correlation ⇒ Gaussian marginals.** The correlated members use the MVN
(Gaussian copula via Cholesky) path; we do not mix the existing Student-t marginal
option (`enable_student_t_hierarchy`) into a correlated block (would require
multivariate-t, out of scope). Singletons keep whatever marginal they have today.

**Identification:**
- `σ_01` — well-fed (497 patients' progression/censoring); tight.
- `σ_03` — weakly fed (57 events); a **learned half-normal**
  (`log_lambda_gp_03_level_intercept_sd[lv] ~ normal(0, ..._sd_sd[lv])`), so it
  **self-shrinks toward 0 if no signal** — worst case is "no worse than today," not
  divergence.
- `ρ_{01,03}` — identified jointly from how patients' long 0→1-censored histories
  co-occur with 0→3 events; LKJ(2) regularizes it off the ±1 boundary.

### Section 4 — R plumbing & spop consistency

Four R touch-points (per the "Adding Module Parameters" checklist):

1. **`r/priors.R`** — `get_multistate_priors()` gains the LKJ shape (`eta = 2`) per
   configured `(level, group)`. Existing σ hyperprior defaults unchanged.
2. **`r/initializers.R`** — init `L_ms_intercept_corr` at the identity Cholesky
   (`diag(d)`) and `z_ms_intercept` at `rnorm(., sd = 0.3)` (zero-spread inits cause
   the known lp=-1e50 stuck-chain failure). Chains start uncorrelated, discover
   correlation in warmup — stable starting geometry.
3. **`targets/publication_targets.R`** — the three new config arrays for all 5
   transitions × `n_levels`; patient level configured as in Section 1.
4. **Config translation helper** (not an alias) — a one-time pipeline-assembly function
   mapping the legacy single `enable_ms_level_baseline_hazard[lv]` to the three
   decomposed arrays, so sclc/pioneer reproduce **bit-identically**. The legacy
   flag is removed from the Stan code; only the three arrays are passed.

**Spop consistency — no GQ change.** The GQ simulator
(`calculate_all_patients_endpoints_rng`, ~868 lines) consumes the assembled
`log_cond_surv_0k` matrices and reuses them for the spop re-roll. The patient-level
frailty is baked into `log_cond_surv_0k` in `transformed_parameters`, so it
automatically rides into the spop re-roll with its posterior value — exactly like
`patient_log_decrease_rate` today. **Zero lines change in the simulator or
`r/util.R`'s CIF recomputation.** This is a payoff of the model's separation between
hazard assembly and event simulation.

**Tests** (`tests/testthat/`):
- **Dimensioning** — the `array[N_TRANS, n_levels]` configs size correlated/singleton
  blocks correctly (extend `test-time-varying-coef-01-dimensioning.R`).
- **Degeneracy** — `corr_group` all-zero ⇒ generated quantities bit-identical to the
  legacy single-flag model. *This is the linchpin: it lets the change land on `main`
  without re-fitting sclc/pioneer.*
- **Validation** — `fatal_error` fires when a corr_group member isn't RE, or members
  don't share a group partition.

### Section 5 — Build, validation sequence & cost

Gated order (each gate passes before the next):

1. **Stan decomposition first, no correlation** — split the flag, wire per-transition
   group-sizing, prove the degeneracy test (bit-identical to legacy). The big, risky
   refactor isolated as a diff where *nothing should change numerically*.
2. **Add the correlated block** — params, MVN assembly, LKJ(2), validation.
3. **Stan syntax check** — `stanc` on `sf-ssm-log-space.stan` and `pioneer.stan`.
4. **Unit tests** — dimensioning + degeneracy + validation.
5. **Prior-predictive gate** — sample the new block from priors, push to spop CIF,
   confirm plausible routing spread (cheap; catches misspecification before MCMC spend).
6. **Full posterior fit + downstream** — Domino via `pioneer-toolkit:start-job` with
   `-D`; then `$diagnostic_summary()` and the `died_off_trial` spop-vs-observed KM
   comparison (gap was −20.5 wk; success = meaningful shrinkage).

**Cost & risk:** the refactor (step 1) is the bulk — touches `transformed_data`,
`parameters`, `transformed_parameters`, `priors`, `hyperparams` for all 5 transitions
plus R priors/initializers/tests; degeneracy test is the safety net. The correlated
block (step 2) is small and additive. New params: one `d=2` Cholesky (1 correlation)
plus a `2 × n_patients` `z` — negligible vs the state-space. Warmup may lengthen
slightly as `ρ`/`σ_03` adapt. If `σ_03 → 0`, the frailty self-disables (no worse than
today) — itself a publishable finding.

## Out of scope

- Multivariate-t marginals for correlated blocks (Gaussian only).
- Correlating covariate *slopes* across transitions (this design correlates intercepts).
- Re-fitting sclc / pioneer (degeneracy test guarantees no change for them).

## Success criterion

The `died_off_trial` cohort's spop median PFS moves materially toward the observed
median (from a −20.5 wk gap), with clean convergence
(`$diagnostic_summary()`: divergences ≪ 1%, E-BFMI > 0.2, rhat ≈ 1).
