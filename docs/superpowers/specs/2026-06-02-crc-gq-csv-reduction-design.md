# CRC Fit CSV Size Reduction — Design

**Date:** 2026-06-02
**Branch:** `karim/crc`
**Status:** Design — awaiting review

> **Status update 2026-06-03:** The Stan refactor (Option B) described below was
> superseded by a simpler R-side CSV projection approach (no Stan changes, no refit).
> Implementation plan: `docs/superpowers/plans/2026-06-03-crc-csv-slim-projection.md`.
> Validated: chain 1 slim (14.78 GB), 0 NA cells across 3,015,000 KM draw cells.

## Problem

CRC publication KM/forecast plots render blank. Root cause (traced, not guessed):
all CRC `*_rvar` targets have the correct shape (2000 draws, 51930 vars) but **all
values are NA**. The NA originates at the CSV read step. CRC fit CSVs are
**~94.7 GB/chain** (hash `00f4a5`); `cmdstanr::read_cmdstan_csv()` (cmdstanr 0.9.0,
CmdStan 2.38.0) cannot parse files this large and silently returns NA draws.
SCLC chains are ~19 GB and parse fine.

### Why CRC is 94 GB

Each draw is **8,268,147 columns** wide (verified by extracting the 266 MB header
line and tallying prefixes). The width is dominated by `patient × time`
transformed-parameter matrices with `n_patients = 2128`, `max_all_t = 200`,
`max_t_width = 218`:

| Block | cols/draw | block | feeds likelihood? | feeds GQ? |
|---|---:|---|:---:|:---:|
| `log_cond_surv_*` (6 transitions) | 1,906,688 | transformed params | ✅ MS likelihood | ✅ endpoint RNG |
| `ms_time_varying_covar_01` | 1,276,800 | transformed params | ✅ (builds hazard) | — |
| `states_full_grid` | 927,808 | transformed params | ✅ (builds `states`) | ✅ |
| `patient_log_*_rate_residual` | 927,808 | transformed params | ✅ (builds states) | ✅ |
| `log_level_lambda_*_residual` | 855,680 | transformed params | ✅ (builds `log_cond_surv`) | — |
| `forecast_patient_states` | 691,424 | generated quantities | — | ✅ (is an endpoint output) |
| `forecast_patient_log*` | 1,382,848 | generated quantities | — | ✅ |

## Constraints (verified)

1. **Stan writes every top-level `transformed parameters` variable to CSV.** No
   annotation suppresses output. The only ways to keep a quantity out of the CSV:
   (a) declare it inside a **local `{ }` block** (scope ends at the brace), or
   (b) compute it in the **`model` block** as a local.
2. **A variable cannot be shared between the `model` block and the `generated
   quantities` block except via `transformed parameters`.** So any quantity used by
   both the likelihood AND the GQ endpoints is structurally forced into the CSV
   unless it is computed **twice** (once as a model-block local, once in GQ).
3. **CmdStan 2.38 / cmdstanr 0.9 have no output-variable-selection option** at
   sample time (`sig_figs`, `save_warmup`, `save_metric` exist; no variable filter).
4. **Forecast/background patient routing is dead code** in every live pipeline
   (`forecast_split_level = 0` always ⇒ `n_forecast_patients == n_patients`).
   It is NOT a usable lever to subset patients, and wiring it into the tumor
   likelihood would drop the historical trials from the fit entirely. Out of scope.
5. **All 2128 patients must stay in the fit** — the historical CRC trials inform the
   shared hierarchy parameters. The parameter space is therefore **unchanged**, so
   the existing job-1850 adapted-metric JSONs remain the right size and **reusable**
   for warm-starting.

## Approach (Option B)

Reduce CSV width by keeping the giant `patient × time` matrices **out of the CSV**,
without changing which patients are fit or the parameter space. Two patterns,
applied per the block's consumer set:

### Pattern 1 — Localize (block feeds likelihood only, NOT GQ)

Move the computation into a **local `{ }` block** inside `transformed parameters`,
emitting only the small downstream quantity the rest of the model needs.

- `ms_time_varying_covar_01` (1.28M cols) — used only to build the MS hazard. Wrap
  in a local scope that produces the hazard contribution; don't expose the matrix.
- `log_level_lambda_*_residual` (856K cols) — used only to build `log_cond_surv`.
  Localize within the MS transformed-parameters scope.
- `states_full_grid` (928K cols) — the **model block uses `states`** (the small
  `[n_forecast_visits, 2]` extracted-at-visits form), NOT the full grid. The grid is
  an intermediate. Build the grid in a local scope, extract `states`, let the grid
  go out of scope. The GQ block recomputes the grid (Pattern 2 cost applies there).
- `patient_log_*_rate_residual` (928K cols) — intermediate to `states`; localize
  alongside the grid build.

### Pattern 2 — Dual-compute (block feeds BOTH likelihood and GQ)

Compute as a **model-block local** for the likelihood; **recompute** in the GQ block
for the endpoints. The recomputation runs once per saved draw (2000×), not per
leapfrog — negligible cost vs sampling. The correctness requirement: the model-local
and the GQ recomputation must be **bit-identical** (shared helper function or
copy-verified code), or the likelihood and reported endpoints silently diverge.

- `log_cond_surv_*` (1.91M cols) — handed in full to both the MS likelihood
  (`ms_final_state ~ multistate(..., log_cond_surv_*, ...)`) and the GQ endpoint RNG
  (`calculate_all_patients_endpoints_rng(..., log_cond_surv_*, ...)`). Move the
  hazard-matrix construction to a model-block local; recompute in GQ.

### Expected result

Removing Pattern-1 blocks (~4.0M cols) + Pattern-2 blocks (~1.9M cols) ≈ 5.9M of
8.27M cols → CSV drops from ~94 GB to **~28–33 GB/chain**, comfortably above SCLC's
working 19 GB. `forecast_patient_*` GQ outputs (the actual forecasts) stay — they
are the deliverable and are inherently `patient × time`.

## Affected files

| File | lines | change |
|---|---:|---|
| `stan/modules/state_space/transformed_parameters.stan` | 188 | localize `states_full_grid`, `patient_log_*_rate_residual` |
| `stan/modules/multistate/transformed_parameters.stan` | 933 | localize `ms_time_varying_covar_01`, `log_level_lambda_*_residual`; move `log_cond_surv_*` construction to model-block local |
| `stan/modules/multistate/cond_surv_transform.stan` | 18 | follows `log_cond_surv` relocation |
| `stan/_ms_burden_tv_covar.stan` | 85 | `ms_time_varying_covar_01` localization |
| `stan/tumor/sf-ssm-log-space.stan` | — | model block gains `log_cond_surv_*` local computation; GQ block recomputes |
| `stan/modules/state_space/generated_quantities.stan` | 89 | recompute `states_full_grid` for GQ |
| `stan/tumor/_tumor_endpoints_generated_quantities.stan` | 486 | recompute `log_cond_surv_*` for endpoint RNG |

**Shared by PSA model:** `state_space` and `multistate` transformed-parameters are
included by `stan/psa/pioneer.stan` too. Changes must preserve PSA behavior
(PSA also has `n_forecast_patients == n_patients`). Stan syntax check both models.

## Verification

1. **Stan compiles** — both `sf-ssm-log-space.stan` and `pioneer.stan` pass
   `stanc` with includes.
2. **Likelihood unchanged** — a short CRC run (few iters, subsample) produces the
   same `lp__` trajectory as the current model on identical seed/inits/data
   (proves the localize/dual-compute refactor is value-preserving).
3. **CSV width** — new header column count ≈ 2.3M (down from 8.27M); chain size
   ~30 GB.
4. **Read succeeds** — `select_draws()` on the new fit returns non-NA draws.
5. **Endpoints match** — GQ-recomputed `log_cond_surv`/`states_full_grid` yield KM
   curves identical (within MC error) to a reference computed the old way.
6. **CRC plots populate** — `quarto render quarto/publication/crc` shows non-blank KM.

## Out of scope

- Forecast/background patient routing (dead code; would drop historical trials).
- Propensity weighting (off in all runs).
- Changing the hierarchy (`trial_arm`, patient) or `patient_trial`/reporting.
- SCLC (unaffected; its 19 GB chains read fine).

## Warm-starting (bundled into this change)

Option B touches only `transformed parameters`, never the `parameters {}` block, so
the unconstrained parameter vector — and thus the mass-matrix length — is
**identical** before and after. The job-1850 adapted-metric JSONs are therefore
size-compatible and are wired into the CRC fit **as part of this change**:

- Copy the four job-1850 per-chain `*_metric.json` files from the fit dir to `data/`
  using the convention `data/inv_metric_crc_chain<N>.json`.
- Add a `metric_files` value (per-chain vector) to the CRC row and pass it via
  `inv_metric` in the `tumor_ssls_res` fit target (per the pioneer warm-start
  rule in `.claude/rules/pioneer.md`).

One refit then delivers both the CSV reduction and the warm-start. The dual-compute
refactor (Pattern 2) is confirmed acceptable; the bit-identical model-local vs
GQ-recompute requirement is the primary correctness gate for review.
