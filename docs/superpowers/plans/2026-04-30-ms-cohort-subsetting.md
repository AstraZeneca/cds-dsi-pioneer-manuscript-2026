# MS-Cohort Subsetting Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a generic `ms_split_level` + `ms_target_groups[]` mechanism so the multistate likelihood can be restricted to a subset of hierarchy groups (initially: trial patients only). PSA likelihood keeps using all patients. Default behavior unchanged.

**Architecture:** Add two new fields to Stan data contract (`ms_split_level`, `ms_target_groups[]`), compute `ms_patient_idx[n_ms_patients]` in transformed data, swap `forecast_patient_idx → ms_patient_idx` in `multistate/likelihood.stan` and `multistate/generated_quantities.stan`. Add sidecar target `ms_patient_idx` for downstream KM/endpoint joins. Add tribble column `ms_split` with helper function that derives trial IDs.

**Tech Stack:** Stan 2.38, R + targets, cmdstanr, testthat

**Reference spec:** `docs/superpowers/specs/2026-04-30-ms-cohort-subsetting.md`

---

## File Structure

### Modified files

| File | Responsibility |
|---|---|
| `stan/modules/multistate/data.stan` | Add `ms_split_level`, `n_ms_target_groups`, `ms_target_groups[]` |
| `stan/modules/multistate/transformed_data.stan` | Compute `ms_is_target[]`, `ms_patient_idx[]`, `n_ms_patients` |
| `stan/modules/multistate/likelihood.stan` | Swap `forecast_patient_idx → ms_patient_idx` |
| `stan/modules/multistate/generated_quantities.stan` | Size GQ outputs by `n_ms_patients`; loop over `ms_patient_idx` |
| `r/pioneer/prepare_analysis_data.R` | Add `ms_split_level`, `ms_target_groups` args to `prepare_pioneer_stan_data()` |
| `targets/pioneer_targets.R` | Add `ms_split` column; helper `ms_args_for_split()`; new variant row; `ms_patient_idx` sidecar target |
| `tests/testthat/test-ms-patient-idx.R` | NEW — unit tests for transformed-data subset logic via Stan's `expose_functions` |

### Out of scope

- `stan/psa/pioneer.stan` top-level model block — no changes needed (likelihood.stan + GQ carry the change)
- `stan/modules/propensity/*` — fully decoupled, untouched
- PSA standalone prep — untouched

---

## Task 1: Add Stan data contract fields

**Files:**
- Modify: `stan/modules/multistate/data.stan`

- [ ] **Step 1: Read current data.stan to understand placement**

Run: `Read stan/modules/multistate/data.stan`
Expected: see current fields like `ms_final_state[n_patients]`, `ms_time_01[n_patients]`, etc.

- [ ] **Step 2: Add new fields at the end of the data block**

Insert after the last existing MS data field:

```stan
// MS-cohort subsetting (generic level/group selector)
// ms_split_level == 0 → MS uses all forecast_patient_idx (default, backwards-compatible)
// ms_split_level > 0  → restrict MS likelihood to patients whose group ID at this
//                       hierarchy level matches any entry in ms_target_groups[].
int<lower=0, upper=n_levels> ms_split_level;
int<lower=0> n_ms_target_groups;
array[n_ms_target_groups] int ms_target_groups;
```

- [ ] **Step 3: Stan compile check**

Run: `~/.cmdstan/cmdstan-2.38.0/bin/stanc --include-paths=stan --include-paths=stan/psa stan/psa/pioneer.stan`
Expected: compiles with no errors.

- [ ] **Step 4: Commit**

```bash
git add stan/modules/multistate/data.stan
git commit -m "feat(multistate): add ms_split_level + ms_target_groups data fields"
```

---

## Task 2: Compute ms_patient_idx in transformed_data

**Files:**
- Modify: `stan/modules/multistate/transformed_data.stan`

- [ ] **Step 1: Read current transformed_data.stan**

Run: `Read stan/modules/multistate/transformed_data.stan`
Expected: see existing transformed-data logic (e.g. `ms_time_*_diff`, etc.).

- [ ] **Step 2: Add ms_is_target + ms_patient_idx computation**

Place this block after the existing transformed-data logic but before the block closes:

```stan
// Compute ms_patient_idx: compacted index of patients contributing to MS likelihood
array[n_patients] int<lower=0, upper=1> ms_is_target = rep_array(0, n_patients);
int n_ms_patients = 0;
if (ms_split_level == 0) {
  // No gating: MS uses every forecast patient
  for (j in 1:n_forecast_patients) {
    ms_is_target[forecast_patient_idx[j]] = 1;
  }
  n_ms_patients = n_forecast_patients;
} else {
  for (j in 1:n_forecast_patients) {
    int p = forecast_patient_idx[j];
    int grp = patient_level_groups[p, ms_split_level];
    for (k in 1:n_ms_target_groups) {
      if (grp == ms_target_groups[k]) {
        ms_is_target[p] = 1;
        n_ms_patients += 1;
        break;
      }
    }
  }
}

// Validate: if fitting MS, at least one patient must contribute
if (fit_multistate_data == 1 && n_ms_patients == 0) {
  fatal_error("MS likelihood enabled but n_ms_patients == 0 (check ms_split_level/ms_target_groups)");
}

array[n_ms_patients] int ms_patient_idx;
{
  int write_idx = 1;
  for (j in 1:n_forecast_patients) {
    int p = forecast_patient_idx[j];
    if (ms_is_target[p] == 1) {
      ms_patient_idx[write_idx] = p;
      write_idx += 1;
    }
  }
}
```

- [ ] **Step 3: Stan compile check**

Run: `~/.cmdstan/cmdstan-2.38.0/bin/stanc --include-paths=stan --include-paths=stan/psa stan/psa/pioneer.stan`
Expected: compiles.

- [ ] **Step 4: Commit**

```bash
git add stan/modules/multistate/transformed_data.stan
git commit -m "feat(multistate): compute ms_patient_idx + ms_is_target in transformed_data"
```

---

## Task 3: Swap index in MS likelihood

**Files:**
- Modify: `stan/modules/multistate/likelihood.stan`

- [ ] **Step 1: Read current likelihood.stan**

Run: `Read stan/modules/multistate/likelihood.stan`
Expected: ~25 lines, single `if (fit_multistate_data)` block calling `~ multistate(...)` with `forecast_patient_idx`.

- [ ] **Step 2: Replace forecast_patient_idx with ms_patient_idx in the ~ multistate call**

Change:
```stan
ms_final_state[forecast_patient_idx] ~ multistate(
  likelihood_weight[forecast_patient_idx],
  ...
);
```
to:
```stan
ms_final_state[ms_patient_idx] ~ multistate(
  likelihood_weight[ms_patient_idx],
  ...
);
```

(Apply to every argument that currently uses `forecast_patient_idx`.)

- [ ] **Step 3: Stan compile check**

Run: `~/.cmdstan/cmdstan-2.38.0/bin/stanc --include-paths=stan --include-paths=stan/psa stan/psa/pioneer.stan`
Expected: compiles.

- [ ] **Step 4: Commit**

```bash
git add stan/modules/multistate/likelihood.stan
git commit -m "feat(multistate): gate likelihood by ms_patient_idx"
```

---

## Task 4: Update MS generated quantities

**Files:**
- Modify: `stan/modules/multistate/generated_quantities.stan`

- [ ] **Step 1: Read current generated_quantities.stan**

Run: `Read stan/modules/multistate/generated_quantities.stan`
Expected: arrays sized by `n_forecast_patients`, loops over `forecast_patient_idx`.

- [ ] **Step 2: Resize output arrays to n_ms_patients and loop over ms_patient_idx**

Swap every occurrence of `n_forecast_patients` → `n_ms_patients` and `forecast_patient_idx[j]` → `ms_patient_idx[j]` **only for MS-related outputs**. Do not change PSA GQ outputs.

- [ ] **Step 3: Stan compile check**

Run: `~/.cmdstan/cmdstan-2.38.0/bin/stanc --include-paths=stan --include-paths=stan/psa stan/psa/pioneer.stan`
Expected: compiles.

- [ ] **Step 4: Commit**

```bash
git add stan/modules/multistate/generated_quantities.stan
git commit -m "feat(multistate): size MS GQ outputs by n_ms_patients"
```

---

## Task 5: R prep function signature

**Files:**
- Modify: `r/pioneer/prepare_analysis_data.R`

- [ ] **Step 1: Read current prepare_pioneer_stan_data signature**

Run: `Grep "^prepare_pioneer_stan_data" r/pioneer/prepare_analysis_data.R`
Expected: find the function start. Read its signature and body structure.

- [ ] **Step 2: Add ms_split_level + ms_target_groups args**

Add to function signature:
```r
ms_split_level    = 0L,
ms_target_groups  = integer(0),
```

Add to `list_assign(...)` at the end:
```r
ms_split_level     = as.integer(ms_split_level),
n_ms_target_groups = length(ms_target_groups),
ms_target_groups   = as.integer(ms_target_groups),
```

- [ ] **Step 3: Commit**

```bash
git add r/pioneer/prepare_analysis_data.R
git commit -m "feat(pioneer): add ms_split_level + ms_target_groups to stan data prep"
```

---

## Task 6: Bit-exact regression test

**Goal:** Run existing `combined` variant with `ms_split_level = 0` (default), confirm identical log-posterior per iteration vs. `main` branch.

- [ ] **Step 1: Get reference log_prob from main**

```bash
git checkout main -- <minimal-reproduction-script>  # or use a cached reference
```

Alternative: run a short 10-iter fit on `main`, record `lp__` from the first 10 iterations of each chain.

- [ ] **Step 2: Run same fit on this branch**

Rebuild just the `psa_standalone_fit_posterior_flatiron_only_no_st` (or a small combined variant) with `ms_split_level = 0`. Record `lp__`.

- [ ] **Step 3: Assert bit-exact equivalence**

```r
stopifnot(all.equal(ref_lp, new_lp, tolerance = 1e-12))
```

- [ ] **Step 4: Commit test or test runner**

```bash
git add tests/testthat/test-ms-subset-regression.R
git commit -m "test(multistate): bit-exact regression for ms_split_level=0"
```

---

## Task 7: Sidecar ms_patient_idx target

**Files:**
- Modify: `targets/pioneer_targets.R`

- [ ] **Step 1: Find the existing MS draws target pattern**

Run: `Grep "ms_standalone_draws_endpoints|pioneer_draws_endpoints" targets/pioneer_targets.R`
Expected: find existing target that reads the MS GQ output.

- [ ] **Step 2: Add sidecar ms_patient_idx target**

Inside the nested `tar_map`, add a target that extracts `ms_patient_idx` from the fit:

```r
tar_target(
  pioneer_ms_patient_idx,
  posterior::as_draws_array(fit$draws(variables = "ms_patient_idx"))[1, 1, ] |> as.integer()
),
```

(The `ms_patient_idx` is a deterministic transformed-data value; a single draw is sufficient.)

Downstream endpoint targets (`pioneer_draws_endpoints_*`, `km_*`) should join on this vector when interpreting MS GQ output.

- [ ] **Step 3: Commit**

```bash
git add targets/pioneer_targets.R
git commit -m "feat(targets): add ms_patient_idx sidecar target"
```

---

## Task 8: Add ms_split column + helper + variant

**Files:**
- Modify: `targets/pioneer_targets.R`

- [ ] **Step 1: Add ms_split column to tribble**

Add column with default `"all"` for every existing row.

- [ ] **Step 2: Define ms_args_for_split() helper**

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

- [ ] **Step 3: Wire ms_args_for_split() into the prep step**

In the data-prep target, inject:
```r
ms_args <- ms_args_for_split(ms_split, analysis_data)
prepare_pioneer_stan_data(
  ...,
  ms_split_level   = ms_args$ms_split_level,
  ms_target_groups = ms_args$ms_target_groups
)
```

- [ ] **Step 4: Add new variant row**

Add to tribble:
```r
"combined_ms_trial_no_st", TRUE, ..., ms_split = "trial",
  enable_student_t_hierarchy = FALSE, ...
```

- [ ] **Step 5: Commit**

```bash
git add targets/pioneer_targets.R
git commit -m "feat(targets): add ms_split column + combined_ms_trial_no_st variant"
```

---

## Task 9: End-to-end smoke test on Domino

- [ ] **Step 1: Push branch**

```bash
git push -u origin karim/partial-ms-fit
```

- [ ] **Step 2: Launch Domino job via /start-job skill**

Target: `pioneer_fit_posterior_combined_ms_trial_no_st` (or equivalent).
Store: `ms-trial-only-no-st` (new store for this experiment).
Hardware: `ods-cpu-r6i-8x` (250 GB) — sufficient for PSA+MS fit.

- [ ] **Step 3: Monitor via /stan-diagnose**

Expected:
- Sampling E-BFMI ≥ 0.3 (comparable to #715)
- Zero divergences (or ≤2 per chain during warmup)
- Maxtree hits < 10 per chain

- [ ] **Step 4: Commit any fixes + document the run**

If issues found during smoke test, fix + re-run. When green:

```bash
git add memory/
git commit -m "docs: record ms_trial_only smoke test results"
```

---

## Task 10: Downstream KM integration

**Goal:** Make sure existing KM / endpoint plots work with compacted MS output.

- [ ] **Step 1: Identify downstream targets consuming MS draws**

Run: `Grep "ms_standalone_draws_endpoints|km_trial" targets/pioneer_targets.R r/`
Expected: list of targets that currently assume full-shape MS output.

- [ ] **Step 2: Update each to join with ms_patient_idx**

For each downstream function that ingests MS draws: add a `patient_idx` column based on `ms_patient_idx` so plots/KM logic can filter/join correctly.

- [ ] **Step 3: Run end-to-end from fit → endpoints → KM plot**

Ensure no silent mis-joins; visually inspect KM for sanity.

- [ ] **Step 4: Commit**

```bash
git add r/ targets/pioneer_targets.R
git commit -m "feat(targets): propagate ms_patient_idx through endpoint/km targets"
```

---

## Remember

- Exact file paths always
- Complete code in every step — if a step changes code, show the code
- Exact commands with expected output
- DRY, YAGNI, TDD, frequent commits
- **Never run data-prep pipelines**. All changes are to Stan + targets wiring.
- **Always push branch before launching Domino job** (skill enforces this).
- **Never hardcode commit SHA** when launching jobs — use `$(git rev-parse HEAD)`.
