# MS-Cohort Subsetting Design Spec

**Branch:** `karim/partial-ms-fit` (from `karim/re-cp`)
**Date:** 2026-04-30
**Goal:** Enable the pioneer model to fit the PSA (biomarker) likelihood on trial + RWD while fitting the multistate (MS) likelihood on a configurable subset of patients (initially, trial only).

## Motivation

The Student-t-off ablation (`karim/re-cp` branch) showed that PSA-only fits sample cleanly on both Flatiron-only and combined trial+RWD data once the Student-t hierarchy is disabled (jobs #705, #715). The next step toward the full joint model is to add the MS likelihood back — but we first want to test whether MS misbehaves when fed RWD patients, whose event times are coarser than trial patients'. The cleanest experiment is: **PSA on trial+RWD, MS on trial only**.

Rather than hardcode a trial-only gate, make the mechanism generic: configure which *hierarchy level* and which *group IDs at that level* contribute to the MS likelihood. This is analogous to how propensity weighting specifies a target group (but propensity weights all patients; MS subsetting excludes them).

## Design Decisions (user-confirmed)

| Decision | Choice |
|---|---|
| Subset mechanism | Generic: `ms_split_level` + `ms_target_groups[]` |
| Subset identification | Use existing hierarchy level 1 (trial_id); pass group IDs |
| Default behavior (no gating) | `ms_split_level = 0` → MS runs on all forecast patients (backwards-compatible) |
| Relation to propensity | Fully decoupled — MS subsetting is independent of propensity weighting |
| Scope | MS likelihood only (propensity unchanged) |
| MS generated quantities | Compute **only** for subset patients (compacted output arrays) |
| Visit gating | Keep `visit_gated_01 = TRUE` for trial-only variant |
| PSA standalone model | Does not expose new fields |
| Tribble API | New column `ms_split` (`"all"` / `"trial"`) drives behavior |

## Stan Data Contract

### New fields in `stan/modules/multistate/data.stan`

```stan
int<lower=0, upper=n_levels> ms_split_level;
int<lower=0> n_ms_target_groups;
array[n_ms_target_groups] int ms_target_groups;
```

### Semantics

- `ms_split_level == 0` → no gating; MS uses all `forecast_patient_idx`. Back-compat.
- `ms_split_level > 0` → patient `i` contributes to MS iff
  `patient_level_groups[i, ms_split_level] ∈ ms_target_groups`.

### Transformed data (in `multistate/transformed_data.stan`)

```stan
array[n_patients] int<lower=0, upper=1> ms_is_target;
int<lower=0, upper=n_forecast_patients> n_ms_patients;
array[/* dynamic */] int ms_patient_idx;
```

- When `ms_split_level == 0`: `ms_is_target[i] = 1` for all `i ∈ forecast_patient_idx`; `ms_patient_idx = forecast_patient_idx`; `n_ms_patients = n_forecast_patients`.
- When `ms_split_level > 0`: iterate `forecast_patient_idx`, check group membership, set `ms_is_target` and build `ms_patient_idx` compactly.
- Validate: `fit_multistate_data == 1 → n_ms_patients > 0`.

## Likelihood (`multistate/likelihood.stan`)

Replace `forecast_patient_idx` with `ms_patient_idx` in the `~ multistate(...)` call. No other logic changes.

## Generated Quantities (`multistate/generated_quantities.stan`)

- Output arrays sized `n_ms_patients` instead of `n_forecast_patients`.
- Loop over `ms_patient_idx`.
- Emit `ms_patient_idx` as a sidecar record so downstream code can map subset rows back to patient IDs.

## Propensity Module

**Untouched.** `likelihood_weight[n_patients]` still computed over all patients. The MS subset and propensity target group are logically independent.

## R Changes

### `r/pioneer/prepare_analysis_data.R`

`prepare_pioneer_stan_data()` gains two args:

```r
prepare_pioneer_stan_data <- function(
  ...,
  ms_split_level = 0L,
  ms_target_groups = integer(0)
) {
  ...
  list_assign(
    ms_split_level     = ms_split_level,
    n_ms_target_groups = length(ms_target_groups),
    ms_target_groups   = as.integer(ms_target_groups)
  )
}
```

`prepare_psa_standalone_stan_data()` sets both to 0/empty — the Stan contract accepts these even for PSA standalone (though the fields are unused there).

### Trial-only wiring (`targets/pioneer_targets.R`)

When the tribble has `ms_split == "trial"`:
- `ms_split_level = 1L` (trial level in the hierarchy)
- `ms_target_groups = <all integer trial IDs except Flatiron>`

The set of "trial" IDs is derived from `analysis_data`: `unique(trial_ids[studyid != "FLATIRON"])`. If this needs to be computed at target-evaluation time, a helper target can produce the vector.

### Downstream KM / endpoint targets

Because MS GQ output is now compacted, downstream code needs a mapping from output row → patient ID. Emit a sidecar target (e.g. `pioneer_ms_patient_idx_<suffix>`) that records `ms_patient_idx` alongside the draws. KM / RMST code joins on this.

## Targets Tribble Extension

New column `ms_split` added to the `tar_map()` tribble. Default `"all"` for existing rows. New variant:

```r
"combined_ms_trial_no_st",
  TRUE,    # hist
  ...,
  ms_split = "trial",
  enable_student_t_hierarchy = FALSE,
  ...
```

The `ms_split` column feeds into the data-prep step via a helper:

```r
ms_args_for_split <- function(split, analysis_data) {
  if (split == "all") {
    list(ms_split_level = 0L, ms_target_groups = integer(0))
  } else if (split == "trial") {
    trial_ids <- analysis_data |>
      distinct(trial_id = as.integer(as_factor(studyid)), studyid) |>
      filter(studyid != "FLATIRON") |>
      pull(trial_id)
    list(ms_split_level = 1L, ms_target_groups = trial_ids)
  } else {
    stop("unknown ms_split: ", split)
  }
}
```

## Out of Scope

- Arm-level or subgroup-level MS gating (generic mechanism supports this; not used yet).
- Propensity-weighted MS contributions for RWD (would be a new weighting option, not a subset).
- Changes to PSA standalone — it still uses its own Stan data prep.

## Open Questions (to resolve during implementation)

1. How does `ms_patient_idx` interact with existing patient-indexed GQ outputs that aren't compacted (e.g., PSA draws)? Likely: PSA draws stay at `n_forecast_patients`; MS draws become `n_ms_patients`. Downstream code uses two different index arrays.
2. When propensity is enabled AND MS subset is trial-only: propensity weights are computed (and applied to PSA) for all patients, but MS sees only trial patients. Is this the intended semantics? (User answer: yes, fully decoupled.)
3. Validate that `ms_target_groups` actually contains IDs present in `patient_level_groups[, ms_split_level]`. Fail loudly in transformed data if the subset is empty.

## Risks

- **Bit-exact regression**: existing variants with `ms_split_level = 0` must produce identical output. The sidecar refactor of GQ outputs risks unexpected shape changes. Mitigate with a bit-exact regression test on one existing combined variant before merging.
- **Downstream KM breakage**: any target that assumes MS draws have shape `n_forecast_patients` will silently mis-join when MS output is compacted. Audit all downstream targets that consume MS draws before running the new variant.
- **Trial ID derivation**: the `trial_ids != FLATIRON` mapping must happen against the same factor ordering Stan sees. Use `patient_level_groups` or the exact factor used to build it — not a separately re-factored vector.
