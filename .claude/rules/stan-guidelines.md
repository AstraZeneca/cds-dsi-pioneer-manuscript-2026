---
paths:
  - "**/*.stan"
---

## Stan Coding Guidelines

- Use built-in zero constructors: `zeros_vector()`, `zeros_int_array()`
- Don't pass array/vector sizes as arguments; use `size()` internally
- Ignore linter warnings about code sections (modular `#include` architecture places code fragments across sections)
- Use non-centered parameterization (NCP) for hierarchical parameters
- Use `fatal_error()` instead of `reject()` for data validation errors in transformed data

## Multistate Architecture Rules

- **All multistate parameters MUST use N-level hierarchy** — never hardcode per-trial (e.g., `log_lambda[patient_trial[i]]`). Use population intercept + `patient_ms_baseline_flat_idx[i, lv]` level shifts instead.
- **Conditional parameter sizing**: Level GP arrays must be `array[enable_ms_XX ? n_levels : 0]` (not unconditionally `array[n_levels]`). The 0→1 transition had this bug — verify every new transition is consistent.
- **`no_oe` model**: Sets `fit_multistate_data=FALSE` — multistate state fields are irrelevant for it. No separate SLD-only state variable needed.
- **`update_dropout_state()`**: Call in targets pipeline after assembling stan data to swap `ms_final_state` → `ms_final_state_dropout` when `enable_ms_03=TRUE`.
- **Pattern E patients**: `death & !progression_before_death & (death_week - patient_max_t > admin_censor_buffer)` — died off-trial; classified as state 3 when `enable_ms_03=TRUE`.
- **Classification-first routing** (`r/sclc/multistate.R`): All multistate Stan fields are derived via a two-step pipeline — never compute `ms_final_state` or transition times inline anywhere else.
  1. `classify_ms_patients(analysis_data)` → adds `ms_pattern` column (factor, 6 levels); called in `prepare_analysis_data()`
  2. `derive_ms_fields(analysis_data, ms_mode)` → named list of all Stan ms fields; called in the targets pipeline as `c(derive_ms_fields(all_analysis_data, ms_mode))`
- **Six patient patterns** (what happened, objective): `admin_censored`, `true_dropout`, `progressed_alive`, `progressed_died`, `died_on_trial`, `died_off_trial`
- **ms_mode** (how model treats it, subjective): `"none"`, `"pfs"`, `"illness_death"`, `"full"` — the pattern→state mapping table is in `multistate.R`
- **`prepare_tumor_stan_data()` does NOT include ms fields** — they are added by `derive_ms_fields()` in the targets pipeline
- **`ms_prog_deterministic`**: Computed in `prepare_analysis_data()` mutate (needs `visit_data`); stored as a column in `analysis_data` and consumed by `derive_ms_fields()`
- **`multistate_lpmf` weight argument**: Always pass the `likelihood_weight` vector — omitting it causes test failures. Signature: `multistate_lpmf(state | ..., likelihood_weight)`

## Propensity Module

`stan/modules/propensity/` implements propensity-weighted borrowing from RWD (pioneer). It follows a **6-file pattern** (no `data.stan` — data lives in `_base_data.stan`):
- `flags.stan` — `enable_propensity_weighting`, `propensity_split_level`, `propensity_target_group`
- `hyperparams.stan`, `parameters.stan` (`beta_propensity`, `beta_propensity_intercept`), `priors.stan`
- `transformed_data.stan` — locates contiguous target patient range; computes `propensity_log_marginal_odds`
- `transformed_parameters.stan` — computes `likelihood_weight[i] = min(1, P(X|trial)/P(X|RWD))`

**Weight formula**: `log_density_ratio = logit(P(trial|X)) - logit(P(trial))`. Capped at 1.0 so no RWD patient outweighs a trial patient. Uses the original `covar_design_matrix` (not QR-decomposed) to keep `beta_propensity` coefficients interpretable.

**Target patients must be contiguous** in the patient array — validated with `fatal_error()` if violated.

## Adding Module Parameters

1. Add feature flag in `modules/<module>/flags.stan`
2. Add hyperparameters in `modules/<module>/hyperparams.stan`
3. Add parameters in `modules/<module>/parameters.stan`
4. Implement transforms in `modules/<module>/transformed_parameters.stan`
5. Add priors in `modules/<module>/priors.stan`
6. Update `r/priors.R` with defaults
7. Update initializer in `r/initializers.R`
8. Update `targets/sclc_targets.R` with flag value
