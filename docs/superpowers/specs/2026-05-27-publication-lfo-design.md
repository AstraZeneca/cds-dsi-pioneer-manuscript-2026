# Publication LFO Pipeline Design

**Date:** 2026-05-27
**Branch:** `karim/pub-lfo`
**Scope:** Add leave-future-out cross-validation to the publication targets pipeline, using sclc as the template.

---

## Context

The publication pipeline (`targets/publication_targets.R`) runs the joint tumor-survival model on the `lilly_cxcr4` target trial + historical controls. It currently has no LFO infrastructure. The sclc pipeline has a fully working LFO setup that we port here with minimal adaptation.

Key differences between publication and sclc data:
- Publication visit data uses `ady` (study day); sclc targets adds `visit_calendar_day = trt_calendar_day + ady - 1` to nested visit data.
- Publication patient data has no `trtsdt`, `patient_first_visit`, or `patient_last_visit`; it has pre-computed `calendar_day` (1-indexed from first patient), `patient_min_t`, `patient_max_t`.
- Publication target trial is `"lilly_cxcr4"` (factor level 1 by construction — first `bind_rows` arg in `prepare_publication_analysis_data()`).

---

## Changes

### 1. Generalize `get_lfo_cutoffs()` in `r/sclc/accuracy.R`

Add a `target_trial` parameter (default `"sclc"`) and a `visit_cal_day_var` parameter (default `visit_calendar_day`). Replace hardcoded `"sclc"` filter, `trtsdt`, `patient_first_visit`, `patient_last_visit` references with the parameterized equivalents.

**New signature:**
```r
get_lfo_cutoffs <- function(all_analysis_data, lfo_step,
                             target_trial = "sclc")
```

**Logic change:** Derive `first_cutoff_day` and `last_cutoff_day` from `calendar_day` (patient-level) on the filtered target trial rows. Recover `cutoff_date` (a real Date) from the origin date:
- If `trtsdt` is present: origin = `min(all_analysis_data$trtsdt) - 1` (preserves sclc behavior).
- Otherwise (publication): origin = `as.Date("1970-01-01") + min(all_analysis_data$calendar_day) - 1` (synthetic but internally consistent — `cutoff_date` is only used for labeling/joining, not Stan data).

The visit-counting step still references `visit_calendar_day` (now present in both pipelines).

**Backward compatibility:** Default args mean sclc calls `get_lfo_cutoffs(all_analysis_data, lfo_step)` unchanged.

### 1b. Add `visit_calendar_day` to publication nested visit data in `r/publication/prepare_analysis_data.R`

Publication visit data has `ady` (per-patient study day, scale −29…440) but NOT `visit_calendar_day` (absolute calendar day). `apply_calendar_cutoff` and `get_lfo_cutoffs` both filter/count on `visit_calendar_day` by default.

Add this `mutate` step inside `prepare_publication_analysis_data()` after the `nest_join`, alongside the other `visit_data` derivations:

```r
visit_data = map2(visit_data, calendar_day, \(d, trt_cal_day) {
  mutate(d, visit_calendar_day = trt_cal_day + ady - 1)
})
```

This mirrors the sclc targets step exactly (`visit_calendar_day = trt_calendar_day + ady - 1`). With `visit_calendar_day` present, all `apply_calendar_cutoff` / `get_lfo_cutoffs` calls use their defaults — no `calendar_day_var = ady` override needed anywhere.

### 2. `publication_targets.R` — top-level additions

Add env-var config and crew controllers, mirroring sclc:

```r
lfo_workers <- as.integer(Sys.getenv("LFO_WORKERS", 24))
controller_lfo <- crew_controller_local(name = "lfo", workers = lfo_workers)

lfo_groups <- Sys.getenv("LFO_GROUPS", 24)
lfo_save_warmup <- Sys.getenv("LFO_SAVE_WARMUP", "false") == "true"
```

Add `controller_lfo` to the `crew_controller_group(...)` call in `tar_option_set()`.

### 3. `publication_targets.R` — LFO Stan model targets

Three new top-level targets (outside any `tar_map`):

```r
tar_target(lfo_tumor_ssls_model_file,
           here("stan", "tumor", "sf-ssls-lfo.stan"), format = "file")
tar_target(lfo_tumor_ssls_include_files,
           find_stan_includes(lfo_tumor_ssls_model_file), format = "file")
tar_target(lfo_tumor_ssls_exe_hash,
           build_model_exe_hash(lfo_tumor_ssls_model_file,
                                lfo_tumor_ssls_include_files,
                                publication_artifacts_path,
                                include_paths = c(here("stan"), here("stan", "tumor"))))
```

### 4. `publication_targets.R` — LFO fit targets

A new section appended after the existing `tar_map` block. These are top-level targets (not inside a `tar_map` — publication has only one model variant and no DCO map):

```r
tar_target(lfo_cutoffs,
           get_lfo_cutoffs(all_analysis_data, lfo_step,
                           target_trial = levels(all_analysis_data$trial)[1]))

tar_target(
  cutoff_all_analysis_data,
  apply_calendar_cutoff(
    all_analysis_data,
    lfo_cutoffs$cutoff_calendar_day,
    require_post_baseline = TRUE
  ) |> mutate(cutoff_calendar_day = lfo_cutoffs$cutoff_calendar_day),
  pattern = map(lfo_cutoffs),
  iteration = "list"
)

tar_group_count(grouped_lfo_cutoffs, lfo_cutoffs, lfo_groups)

tar_target(
  tumor_ssls_lfo,
  lfo(
    base_tumor_ssls_stan_data,
    lfo_tumor_ssls_exe_hash$exe_file,
    grouped_lfo_cutoffs,
    lfo_cutoffs,
    publication_output_path,
    "lfo_tumor_ssls",
    NULL,
    fit_output_timestamp,
    verbose = TRUE,
    fit_only = FALSE,
    exact = TRUE,
    iter_warmup = 500,
    iter_sampling = 500,
    save_warmup = lfo_save_warmup,
    parallel_chains = 4,
    adapt_delta = 0.8,
    threads_per_chain = base_tumor_ssls_stan_data$n_shards,
    future_window = 2,
    initializer_factory = function(stan_data, save_dir, run_id) {
      source(initializers_fixed_file)
      create_tumor_ssls_initializer_fixed(stan_data, save_dir, run_id)
    }
  ),
  resources = tar_resources(crew = tar_resources_crew(controller = "lfo")),
  pattern = map(grouped_lfo_cutoffs)
)

tar_target(
  tumor_ssls_lfo_clean,
  clean_lfo_results(tumor_ssls_lfo) |>
    select(refit_n, n, starts_with("E_"), fit) |>
    left_join(lfo_cutoffs, by = "n")
)
tar_target(tumor_ssls_lfo_clean_lite, select(tumor_ssls_lfo_clean, !fit))
tar_target(
  full_oos_confusion_matrix,
  get_oos_confusion_marix(
    tumor_ssls_lfo_clean,
    recover_data = all_analysis_data |>
      select(trial, visit_data) |>
      unnest(visit_data) |>
      transmute(trial, response, pred_response = response)
  ),
  pattern = map(tumor_ssls_lfo_clean)
)
tar_target(
  oos_confusion_matrix,
  full_oos_confusion_matrix |>
    group_by(n, m, response) |>
    mutate(
      total = rvar_sum(oos_recist_confusion_matrix),
      prop = oos_recist_confusion_matrix / total
    ) |>
    group_by(response, pred_response) |>
    summarize(
      count = rvar_sum(oos_recist_confusion_matrix),
      mean_pred = rvar_weighted_mean(prop, w),
      cell_size = n(),
      .groups = "drop"
    ) |>
    mutate(se = sd(mean_pred) / sqrt(cell_size))
)
```

**Note:** `lfo_step` does not currently exist in `publication_targets.R` — add `tar_target(lfo_step, 30L)` alongside the LFO Stan model targets in Section 3.

---

## What is NOT included

- No `model-validation.qmd` Quarto page (out of scope for this task).
- No `cutoff_km_trial_pfs` targets (not needed for LFO itself).
- No LFO endpoints targets (commented out in sclc too).

---

## Open Questions

None — design is complete.
