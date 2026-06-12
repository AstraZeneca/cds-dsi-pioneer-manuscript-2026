# Publication LFO Pipeline Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add leave-future-out (LFO) cross-validation targets to the publication pipeline (`targets/publication_targets.R`), mirroring the working sclc LFO setup.

**Architecture:** Three files change: (1) `r/publication/prepare_analysis_data.R` gains `visit_calendar_day` on nested visit data so LFO cutoff comparisons work; (2) `r/sclc/accuracy.R`'s `get_lfo_cutoffs()` is generalized with a `target_trial` parameter so it works for any trial name; (3) `targets/publication_targets.R` gains crew controllers, Stan model file targets, and LFO fit targets that mirror sclc exactly. The Stan LFO model (`sf-ssls-lfo.stan`) already exists and is unchanged — it handles the training/test boundary internally via `cutoff_calendar_day`.

**Tech Stack:** R, targets, crew, cmdstanr, testthat, tidyverse

---

### Task 1: Add `visit_calendar_day` to publication nested visit data

The LFO machinery (both `apply_calendar_cutoff` and `get_lfo_cutoffs`) filters visit rows using `visit_calendar_day` — an absolute calendar day on a shared scale across all patients. Publication visit data only has `ady` (per-patient study day). This task adds `visit_calendar_day = calendar_day + ady - 1` inside `prepare_publication_analysis_data()`, exactly mirroring what sclc does in `targets/sclc_targets.R:718`.

**Files:**
- Modify: `r/publication/prepare_analysis_data.R:46-76`
- Test: `tests/testthat/test-publication-lfo.R` (create)

- [ ] **Step 1: Write a failing test for `visit_calendar_day` presence**

Create `tests/testthat/test-publication-lfo.R`:

```r
library(testthat)

source(here::here("r/util.R"))
source(here::here("r/sclc/priors.R"))
source(here::here("r/sclc/prepare_analysis_data.R"))
source(here::here("r/publication/prepare_analysis_data.R"))

make_minimal_pub_data <- function() {
  target_patient <- tibble(
    studyid = "S1", usubjid = "P1", trial = "lilly_cxcr4",
    calendar_day = 100L, calendar_week = 15L,
    pfs = 10L, death = FALSE, death_week = NA_integer_,
    progression_before_death = NA, right_censored = 1L,
    interval_censored = 0L, age = 60, sex = 0L, ecog = 1L,
    hgb = 120, ldh_log = 5, albumin = 40,
    race = "White", stage = "IV", prev_lines = 1L, smoker = TRUE,
    patient_min_t = 1L, patient_max_t = 10L, patient_t_width = 10L
  )
  historical_patient <- target_patient |>
    mutate(usubjid = "P2", trial = "amgen_darbe", calendar_day = 1L, calendar_week = 1L)

  target_visit <- tibble(
    studyid = "S1", usubjid = "P1",
    visitnum = 1L, ady = 1L, week = 1L, mmsumdiam = 50, response = "SD"
  )
  historical_visit <- target_visit |>
    mutate(usubjid = "P2")

  list(
    target_patient = target_patient,
    historical_patient = historical_patient,
    target_visit = target_visit |> determine_visit_data_response(),
    historical_visit = historical_visit |> determine_visit_data_response()
  )
}

test_that("prepare_publication_analysis_data adds visit_calendar_day to nested visit data", {
  d <- make_minimal_pub_data()
  result <- prepare_publication_analysis_data(
    d$target_patient, d$historical_patient,
    d$target_visit, d$historical_visit
  )
  visit_cols <- result |> pull(visit_data) |> purrr::map(colnames) |> purrr::reduce(intersect)
  expect_true("visit_calendar_day" %in% visit_cols)
})

test_that("visit_calendar_day equals calendar_day + ady - 1", {
  d <- make_minimal_pub_data()
  result <- prepare_publication_analysis_data(
    d$target_patient, d$historical_patient,
    d$target_visit, d$historical_visit
  )
  # Target patient: calendar_day=100, ady=1 => visit_calendar_day = 100 + 1 - 1 = 100
  target_row <- result |> filter(fct_match(trial, "lilly_cxcr4"))
  expect_equal(
    target_row$visit_data[[1]]$visit_calendar_day,
    100L + target_row$visit_data[[1]]$ady - 1L
  )
})
```

- [ ] **Step 2: Run to confirm it fails**

```bash
Rscript -e 'testthat::test_file("tests/testthat/test-publication-lfo.R")'
```

Expected: FAIL — `"visit_calendar_day" %in% visit_cols` is FALSE.

- [ ] **Step 3: Add `visit_calendar_day` derivation in `prepare_publication_analysis_data()`**

In `r/publication/prepare_analysis_data.R`, inside the `mutate(...)` block (after the `pfs` adjustment on line 75, before the closing `|>`), add:

```r
      # Lift per-patient ady onto the shared absolute calendar-day scale,
      # needed by apply_calendar_cutoff / get_lfo_cutoffs
      visit_data = map2(visit_data, calendar_day, \(d, trt_cal_day) {
        mutate(d, visit_calendar_day = trt_cal_day + ady - 1L)
      }),
```

The full `mutate()` block in `prepare_publication_analysis_data()` should now end:

```r
    mutate(
      trial = as_factor(trial),
      map_dfr(visit_data, \(d) determine_pfs(d, 1)),
      male = sex,
      across(
        c(ecog, hgb, ldh_log, albumin),
        \(x) replace_na(x, median(x, na.rm = TRUE))
      ),
      original_right_censored = right_censored,
      potential_followup = patient_max_t,
      ms_prog_deterministic = as.integer(replace_na(
        map_lgl(visit_data, \(d) any(fct_match(d$det_response, "PD"), na.rm = TRUE)) &
          replace_na(progression_before_death, FALSE),
        FALSE
      )),
      pfs = if_else(!right_censored, pfs + interval_censored + 1L, pfs),
      visit_data = map2(visit_data, calendar_day, \(d, trt_cal_day) {
        mutate(d, visit_calendar_day = trt_cal_day + ady - 1L)
      }),
    ) |>
```

- [ ] **Step 4: Run tests to confirm they pass**

```bash
Rscript -e 'testthat::test_file("tests/testthat/test-publication-lfo.R")'
```

Expected: 2 tests PASS.

- [ ] **Step 5: Commit**

```bash
git add r/publication/prepare_analysis_data.R tests/testthat/test-publication-lfo.R
git commit -m "feat(publication): add visit_calendar_day to nested visit data for LFO support"
```

---

### Task 2: Generalize `get_lfo_cutoffs()` to accept a `target_trial` parameter

The current `get_lfo_cutoffs()` in `r/sclc/accuracy.R` hardcodes `"sclc"` and uses `trtsdt` / `patient_first_visit` / `patient_last_visit` — columns absent from publication data. This task replaces those with `calendar_day` / `patient_max_t`, adds a `target_trial` parameter, and keeps the sclc call site unchanged (default `"sclc"`).

The `cutoff_date` column (a real `Date`, used only for labelling/joining in display targets) is recovered via: if `trtsdt` is present in the data, use `min(trtsdt) - 1` as origin (preserving sclc's current behaviour); otherwise use a synthetic origin of `as.Date("1970-01-01") + min(all_analysis_data$calendar_day) - 1`.

**Files:**
- Modify: `r/sclc/accuracy.R:29-63`
- Test: `tests/testthat/test-publication-lfo.R` (add cases)

- [ ] **Step 1: Write failing tests for generalized `get_lfo_cutoffs()`**

Add to `tests/testthat/test-publication-lfo.R`:

```r
source(here::here("r/accuracy.R"))

test_that("get_lfo_cutoffs works for a non-sclc target_trial", {
  # Two patients: target trial lilly_cxcr4 (calendar_day 100, 130),
  # historical amgen_darbe (calendar_day 1)
  d <- tibble(
    trial = as_factor(c("lilly_cxcr4", "lilly_cxcr4", "amgen_darbe")),
    calendar_day = c(100L, 130L, 1L),
    patient_max_t = c(10L, 5L, 20L),
    visit_data = list(
      tibble(ady = c(1L, 7L), week = c(1L, 2L),
             visit_calendar_day = 100L + c(1L, 7L) - 1L),
      tibble(ady = c(1L, 7L), week = c(1L, 2L),
             visit_calendar_day = 130L + c(1L, 7L) - 1L),
      tibble(ady = c(1L, 7L), week = c(1L, 2L),
             visit_calendar_day = 1L + c(1L, 7L) - 1L)
    )
  )
  result <- get_lfo_cutoffs(d, lfo_step = 30L,
                             target_trial = "lilly_cxcr4")
  expect_s3_class(result, "data.frame")
  expect_true(all(c("n", "cutoff_date", "cutoff_calendar_day") %in% names(result)))
  expect_true(nrow(result) >= 1L)
  # All cutoff_calendar_days should be within the target trial's calendar_day range
  expect_true(all(result$cutoff_calendar_day >= min(d$calendar_day[d$trial == "lilly_cxcr4"])))
})

test_that("get_lfo_cutoffs default target_trial='sclc' still works", {
  # Sclc data has trtsdt — origin recovery branch must fire
  d <- tibble(
    trial = as_factor(c("sclc", "historical")),
    trtsdt = as.Date(c("2020-01-01", "2019-06-01")),
    calendar_day = c(215L, 1L),
    patient_max_t = c(10L, 20L),
    visit_data = list(
      tibble(ady = c(1L, 7L), week = c(1L, 2L),
             visit_calendar_day = 215L + c(1L, 7L) - 1L),
      tibble(ady = c(1L, 7L), week = c(1L, 2L),
             visit_calendar_day = 1L + c(1L, 7L) - 1L)
    )
  )
  result <- get_lfo_cutoffs(d, lfo_step = 30L)
  expect_s3_class(result, "data.frame")
  expect_true(nrow(result) >= 1L)
  # cutoff_date should be a real Date
  expect_s3_class(result$cutoff_date, "Date")
})
```

- [ ] **Step 2: Run to confirm it fails**

```bash
Rscript -e 'testthat::test_file("tests/testthat/test-publication-lfo.R")'
```

Expected: new tests FAIL — `unused argument (target_trial = "lilly_cxcr4")`.

- [ ] **Step 3: Rewrite `get_lfo_cutoffs()` in `r/sclc/accuracy.R`**

Replace lines 29–63 with:

```r
get_lfo_cutoffs <- function(all_analysis_data, lfo_step, target_trial = "sclc") {
  target_data <- filter(all_analysis_data, fct_match(trial, target_trial))

  first_cutoff_day <- min(target_data$calendar_day)
  last_cutoff_day  <- max(target_data$calendar_day) +
    max(target_data$patient_max_t) * 7L

  # Recover a real Date origin for the cutoff_date label column.
  # cutoff_date is only used for display/joining — it does not enter Stan.
  if ("trtsdt" %in% names(all_analysis_data)) {
    origin <- min(all_analysis_data$trtsdt) - 1L
  } else {
    origin <- as.Date("1970-01-01") + min(all_analysis_data$calendar_day) - 1L
  }

  all_visits_data <- all_analysis_data |>
    select(visit_data) |>
    unnest(visit_data)

  get_lfo_cutoff_days(
    origin + first_cutoff_day,
    origin + last_cutoff_day,
    first_cutoff_day,
    lfo_step
  ) |>
    mutate(
      n_visits_added = map_int(cutoff_calendar_day, \(cutoff_day) {
        all_visits_data |> filter(visit_calendar_day <= cutoff_day) |> nrow()
      }),
      n_future_visits = nrow(all_visits_data) - n_visits_added,
      n_visits_added  = n_visits_added - lag(n_visits_added)
    ) |>
    filter(n_future_visits > 0)
}
```

- [ ] **Step 4: Run all LFO-related tests to confirm they pass**

```bash
Rscript -e 'testthat::test_file("tests/testthat/test-publication-lfo.R")'
```

Expected: all 4 tests PASS.

- [ ] **Step 5: Run the full test suite to confirm no regressions**

```bash
Rscript -e 'testthat::test_dir("tests/testthat")'
```

Expected: all tests pass (same count as before).

- [ ] **Step 6: Commit**

```bash
git add r/sclc/accuracy.R tests/testthat/test-publication-lfo.R
git commit -m "feat: generalize get_lfo_cutoffs() to accept target_trial parameter"
```

---

### Task 3: Add LFO crew controllers and env-var config to `publication_targets.R`

The LFO targets run in parallel across cutoffs using a dedicated crew controller. This task adds the controller and its env-var-controlled worker count, mirroring the sclc pattern.

**Files:**
- Modify: `targets/publication_targets.R:27-56`

- [ ] **Step 1: Add `lfo_workers`, `controller_lfo`, `lfo_groups`, `lfo_save_warmup` to `publication_targets.R`**

After line 29 (`controller_many_samples <- ...`), add:

```r
lfo_workers <- as.integer(Sys.getenv("LFO_WORKERS", 24))
controller_lfo <- crew_controller_local(name = "lfo", workers = lfo_workers)

lfo_groups <- Sys.getenv("LFO_GROUPS", 24)
lfo_save_warmup <- Sys.getenv("LFO_SAVE_WARMUP", "false") == "true"
```

- [ ] **Step 2: Add `controller_lfo` to the `crew_controller_group(...)` call**

The current call (lines 45–51) is:

```r
  controller = crew_controller_group(
    controller_default,
    controller_fit,
    controller_many_samples
  ),
```

Change it to:

```r
  controller = crew_controller_group(
    controller_default,
    controller_fit,
    controller_many_samples,
    controller_lfo
  ),
```

- [ ] **Step 3: Verify the file parses without error**

```bash
Rscript --no-init-file -e 'source("targets/publication_targets.R")' 2>&1 | grep -i "error\|Error" | head -5
```

Expected: no errors.

- [ ] **Step 4: Commit**

```bash
git add targets/publication_targets.R
git commit -m "feat(publication): add LFO crew controller and env-var config"
```

---

### Task 4: Add LFO Stan model file targets and `lfo_step` to `publication_targets.R`

Before running LFO, the pipeline needs to track the Stan LFO model file (so targets invalidates when the Stan code changes) and compile it. This mirrors the three `lfo_tumor_ssls_*` targets in sclc. We also add `lfo_step = 30L` (the cutoff increment in days).

**Files:**
- Modify: `targets/publication_targets.R:80-103`

- [ ] **Step 1: Add `lfo_step` and the three LFO Stan model targets after the existing model targets**

After the closing `),` of the `tumor_ssls_exe_hash` target (around line 103), add:

```r
  tar_target(lfo_step, 30L),

  tar_target(
    lfo_tumor_ssls_model_file,
    here("stan", "tumor", "sf-ssls-lfo.stan"),
    format = "file"
  ),
  tar_target(
    lfo_tumor_ssls_include_files,
    find_stan_includes(lfo_tumor_ssls_model_file),
    format = "file"
  ),
  tar_target(
    lfo_tumor_ssls_exe_hash,
    build_model_exe_hash(
      lfo_tumor_ssls_model_file,
      lfo_tumor_ssls_include_files,
      publication_artifacts_path,
      include_paths = c(here("stan"), here("stan", "tumor"))
    )
  ),
```

- [ ] **Step 2: Verify the Stan LFO model file exists at the expected path**

```bash
ls stan/tumor/sf-ssls-lfo.stan
```

Expected: file is listed.

- [ ] **Step 3: Verify the file parses without error**

```bash
Rscript --no-init-file -e 'source("targets/publication_targets.R")' 2>&1 | grep -i "error\|Error" | head -5
```

Expected: no errors.

- [ ] **Step 4: Commit**

```bash
git add targets/publication_targets.R
git commit -m "feat(publication): add lfo_step and LFO Stan model file targets"
```

---

### Task 5: Add LFO fit and results targets to `publication_targets.R`

This is the core LFO pipeline: compute cutoffs, run fits in parallel across cutoff groups, clean results, and aggregate the OOS confusion matrix. These targets are appended after the existing `tar_map` / combined targets block, at the end of the `publication_targets` list (before the final closing `)`).

**Files:**
- Modify: `targets/publication_targets.R` — end of file, before final `)`

- [ ] **Step 1: Add the LFO targets block at the end of the `publication_targets` list**

Find the last `tar_target(...)` in the list (the `all_tumor_ssls_patient_decrease_prop_bpi` target, ending around line 682). After its closing `),`, add:

```r
  # LFO cross-validation -------------------------------------------------------

  tar_target(
    lfo_cutoffs,
    get_lfo_cutoffs(
      all_analysis_data,
      lfo_step,
      target_trial = levels(all_analysis_data$trial)[1]
    )
  ),

  tar_target(
    cutoff_all_analysis_data,
    apply_calendar_cutoff(
      all_analysis_data,
      lfo_cutoffs$cutoff_calendar_day,
      require_post_baseline = TRUE
    ) |>
      mutate(cutoff_calendar_day = lfo_cutoffs$cutoff_calendar_day),
    pattern = map(lfo_cutoffs),
    iteration = "list"
  ),

  tar_group_count(grouped_lfo_cutoffs, lfo_cutoffs, lfo_groups),

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
  ),

  tar_target(
    tumor_ssls_lfo_clean,
    clean_lfo_results(tumor_ssls_lfo) |>
      select(refit_n, n, starts_with("E_"), fit) |>
      left_join(lfo_cutoffs, by = "n")
  ),

  tar_target(
    tumor_ssls_lfo_clean_lite,
    select(tumor_ssls_lfo_clean, !fit)
  ),

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
  ),

  tar_target(
    oos_confusion_matrix,
    full_oos_confusion_matrix |>
      group_by(n, m, response) |>
      mutate(
        total = rvar_sum(oos_recist_confusion_matrix),
        prop  = oos_recist_confusion_matrix / total
      ) |>
      group_by(response, pred_response) |>
      summarize(
        count     = rvar_sum(oos_recist_confusion_matrix),
        mean_pred = rvar_weighted_mean(prop, w),
        cell_size = n(),
        .groups   = "drop"
      ) |>
      mutate(se = sd(mean_pred) / sqrt(cell_size))
  ),
```

- [ ] **Step 2: Verify the file parses without error**

```bash
Rscript --no-init-file -e 'source("targets/publication_targets.R")' 2>&1 | grep -i "error\|Error" | head -5
```

Expected: no errors.

- [ ] **Step 3: Validate the targets pipeline recognises the new targets**

```bash
TAR_PROJECT=publication Rscript --no-init-file -e '
  source(".Rprofile")
  source("targets/publication_targets.R")
  cat("New LFO targets:\n")
  print(grep("lfo|cutoff", sapply(publication_targets, \(t) t$settings$name), value = TRUE))
' 2>&1 | grep -v "^Warning\|conflicted\|prefer\|Attaching\|package\|mask"
```

Expected output includes: `lfo_cutoffs`, `cutoff_all_analysis_data`, `grouped_lfo_cutoffs`, `tumor_ssls_lfo`, `tumor_ssls_lfo_clean`, `tumor_ssls_lfo_clean_lite`, `full_oos_confusion_matrix`, `oos_confusion_matrix`.

- [ ] **Step 4: Run the full test suite one final time**

```bash
Rscript -e 'testthat::test_dir("tests/testthat")'
```

Expected: all tests pass.

- [ ] **Step 5: Commit**

```bash
git add targets/publication_targets.R
git commit -m "feat(publication): add LFO fit and results targets"
```
