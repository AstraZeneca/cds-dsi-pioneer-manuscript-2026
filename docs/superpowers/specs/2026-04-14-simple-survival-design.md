# Simple Survival Model — Design Spec

**Date**: 2026-04-14
**Branch**: `karim/simple-survival`
**Project**: Pioneer

---

## Overview

A landmark overall survival model using the existing `ms-standalone.stan` infrastructure,
restricted to the 0→2 (direct death) transition only. Two PSA-derived covariates are added
alongside the standard baseline time-invariant covariate set, and the study is landmarked
at Week 13 (Day 85, Cycle 4 Day 1).

---

## Motivation

The full pioneer joint model ties PSA dynamics to survival outcomes. This simpler model
asks a more direct question: conditional on surviving to Week 13 (when the first
post-treatment PSA is observed), how much of the remaining survival hazard is explained by
(a) baseline PSA and (b) the first post-treatment PSA response, over and above the standard
clinical covariates?

---

## Landmark Analysis

**Landmark time**: Week 13 (Day 85 ≈ Cycle 4 Day 1).

**Rationale**: The trial protocol schedules the first post-treatment PSA assessment in the
window Day 64–Day 85. There are no protocol-scheduled PSA draws between Days 49 and 64,
making Day 85 the natural, well-defined first post-treatment PSA timepoint that all
on-study patients share.

**Procedure**:
1. Exclude patients who died or were censored before Week 13 — they do not have the
   post-treatment PSA covariate and do not contribute to this analysis.
2. Rebase all event/censoring times: `t_rebased = t_original − 13` (weeks).
3. Baseline PSA is drawn from the window [Day −90, Day +7] relative to treatment start.
4. First post-treatment PSA is the first observed PSA from Day 64 onwards.

---

## Model

### Stan

`stan/ms-standalone.stan` — **no changes**. All configuration is via data-block flags.

**Flag overrides** (relative to current ms-standalone defaults):

| Flag | Value | Effect |
|------|-------|--------|
| `enable_ms_01` | 0 | Disable 0→1 (progression) |
| `enable_ms_02` | 1 | Enable 0→2 (death) — only active transition |
| `enable_ms_12` | 0 | Disable 1→2 |
| `enable_ms_03` | 0 | Disable 0→3 (dropout) |
| `enable_ms_32` | 0 | Disable 3→2 |
| `enable_ms_pop_time_invariant_cov` | 1 | Enable time-invariant covariate effects on hazard |
| `enable_ms_pop_time_varying_cov` | 0 | No time-varying covariates |

### Covariates

`n_time_invariant_covar` = standard pioneer covariate count + 1.

**Covariate matrix** (columns, in order):
1. All existing pioneer time-invariant covariates processed through the standard recipe:
   - Demographics: `age`, `race`, `arm`
   - Baseline labs: `baseline_psa`, `baseline_albumin`, `baseline_ALP`, `baseline_AST`,
     `baseline_hemoglobin`, `baseline_ldh`, `baseline_neutrophils`, `baseline_nlr`
   - Metastatic sites: `liver_mets`, `bone_mets`, `lung_mets`, `lymph_mets`, `visceral_mets`
   - Treatment history: `prev_lines`, `chemo_flag`, `prior_arpi`
2. **New**: `log_first_posttreat_psa` — log PSA from first visit ≥ Day 64, standardised

Processing: same centering/scaling/QR decomposition as the existing pioneer recipe.
Continuous predictors are scaled; binary/ordinal indicators are centred only.

---

## R Implementation

### New function: `prepare_simple_survival_stan_data()`

Location: `r/pioneer/prepare_analysis_data.R` (or a new
`r/pioneer/prepare_simple_survival_data.R`).

**Steps**:
1. Extract baseline PSA per patient: latest PSA in [Day −90, Day +7].
2. Extract first post-treatment PSA per patient: first PSA from Day 64 onwards.
3. Apply landmark:
   - Drop patients with `ms_time_02 <= 13` and `ms_censored_02 == 1` (censored before landmark)
   - Drop patients with `ms_final_state == 2` and `ms_time_02 <= 13` (died before landmark)
   - Rebase: `ms_time_02 = ms_time_02 - 13`
4. Build covariate matrix: run existing pioneer recipes on the landmarked patient set,
   then append `log(first_posttreat_psa)` (centred/scaled) as a final column before QR.
5. Call existing `prepare_ms_standalone_stan_data()` with:
   - Landmarked patient data
   - Updated covariate matrix
   - Flag overrides listed above

### No changes to `ms-standalone.stan`

### No changes to existing `prepare_ms_standalone_stan_data()`

---

## Targets Pipeline

New entries in `targets/pioneer_targets.R` (or a dedicated
`targets/simple_survival_targets.R`):

```r
tar_target(simple_survival_patient_data,
  prepare_simple_survival_patient_data(pioneer_patient_data, visit_data)),

tar_target(simple_survival_stan_data_posterior,
  prepare_simple_survival_stan_data(simple_survival_patient_data, fit_multistate_data = 1L)),

tar_target(simple_survival_stan_data_prior,
  prepare_simple_survival_stan_data(simple_survival_patient_data, fit_multistate_data = 0L)),

tar_target(simple_survival_fit_posterior,
  sample_and_save(ms_standalone_exe, simple_survival_stan_data_posterior, ...)),

tar_target(simple_survival_fit_prior,
  sample_and_save(ms_standalone_exe, simple_survival_stan_data_prior, ...)),
```

Reuses: `ms_standalone_exe` (already compiled), `sample_and_save()`, existing patient/visit
data targets.

---

## What Is Not Changing

- `stan/ms-standalone.stan` — untouched
- `prepare_ms_standalone_stan_data()` — untouched
- Existing pioneer model variants — untouched
- Stan multistate module — untouched

---

## Out of Scope

- Post-progression (1→2) or dropout (0→3/3→2) transitions
- Time-varying PSA trajectory as a covariate
- Propensity weighting (can be added later if needed)
- Landmark sensitivity analysis (different landmark times)
