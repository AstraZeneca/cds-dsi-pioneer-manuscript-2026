# OOS RECIST Progression-Only Stamp Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the LFO out-of-sample RECIST confusion matrix stamp PD only on genuine progression (0→1), not on death (0→2) or dropout (0→3), so the *predicted* RECIST axis matches the radiographic, scan-only *observed* axis. PFS/OS and all other endpoints keep treating death as an event — the change is confined to the confusion-matrix prediction.

**Architecture:** The shared RNG endpoint function `calculate_all_patients_endpoints_rng` already computes, per patient, a progression-only PFS pair (`sample_t01`, `sample_c01`) that bundles target-lesion PD and the MS 0→1 hazard while *excluding* the death (0→2) and dropout (0→3) competing risks. We expose this pair as two new return-tuple elements, thread it through all six call sites, and switch the LFO RECIST stamp in `sf-ssls-lfo.stan` from `sample_ms_pfs`/`sample_ms_right_censored` to the new progression-only arrays. Everything else that consumes the tuple ignores the two new fields.

**Tech Stack:** Stan (CmdStan 2.38, modular `#include`), R (testthat, cmdstanr), targets pipeline.

---

## Background — Why This Change

Diagnosis established in the investigation that produced this plan (all verified against SCLC LFO job #1894, store `lfo`):

1. **MS hazard is well-calibrated** — forecast MS-PFS KM median ≈28–30 wk vs observed 28 wk for eval trial `lilly_cxcr4`. The hazard is NOT too high. (Concern 1 cleared.)
2. **The confusion matrix over-predicts PD** (53% predicted vs 19% observed) because the *predicted* RECIST axis stamps PD from `sample_ms_pfs`, which is operational PFS = progression OR death OR dropout. The *observed* axis is scan-only RECIST (`pull(response)` from `visit_data`) and never encodes death-after-last-scan as PD. 21/53 eval-trial deaths have no recorded PD (last scan PR/CR/SD). The axes use incompatible definitions of PD. (Concern 2 — the defect.)
3. **Visit gating does NOT misalign weeks** — predicted and observed RECIST are indexed off the same observed-visit array; scoring is correctly confined to real scan weeks. No code change needed for week alignment. (Concern 3 cleared.)

**User's decision:** "The model-generated RECIST prediction should not include death as a PD but only for the purposes of the confusion matrix: PFS and everything else should still treat death as a progress event. Don't inject any weeks."

**Why `sample_t01`/`sample_c01` is the correct signal** (`stan/pfs.stanfunctions:2042-2050`):
```stan
int sample_t01; int sample_c01;
if (!right_censored[p]) {
  sample_t01 = pfs[p];
  sample_c01 = 0;
} else {
  sample_t01 = min(sample_target_pfs[j], sample_ms_pfs[j]);   // sample_ms_pfs here is RAW MS 0→1, pre-graft
  sample_c01 = sample_target_right_censored[j] && sample_ms_right_censored[j];
}
```
At this point in the function (before the death/dropout graft at lines 2078–2094), `sample_ms_pfs[j]` is still the raw 0→1 multistate hazard sample. So `sample_t01` = first of {target-lesion PD, MS 0→1 progression}, and `sample_c01` = censored only if BOTH progression channels are censored. This:
- includes SLD target progression ✓
- includes MS 0→1 progression, which captures non-target/new-lesion progression (the 0→1 transition is fit on overall progression) ✓
- excludes direct death 0→2 ✓ (that is `sample_cause==2`)
- excludes dropout 0→3 ✓ (`sample_cause==3`)
- keeps PD for post-progression death 1→2 ✓ (0→1 fired first, so `sample_t01` is an event)

This is exactly the user's rule. Note `sample_t01`/`sample_c01` are NOT currently saved or returned — they are locals consumed only by `classify_sample_exit`.

---

## File Structure

- **`stan/pfs.stanfunctions`** — defines `calculate_all_patients_endpoints_rng`. Add 2 elements to the return tuple (`sample_prog_pfs`, `sample_prog_right_censored`); declare + populate them from `sample_t01`/`sample_c01`. This is the only file where the new values are *computed*.
- **`stan/multistate.stanfunctions`** — `compute_ms_os_km_rng` canonical overload calls the function and unpacks its tuple. Add 2 ignored slots to its unpack.
- **`stan/_ms_standalone_lfo_os_km_generated_quantities.stan`** — call site; add 2 ignored unpack slots.
- **`stan/modules/state_space/burden_endpoints.stan`** — call site (pioneer); add 2 ignored unpack slots.
- **`stan/tumor/_tumor_endpoints_generated_quantities.stan`** — call site (full tumor model); add 2 ignored unpack slots.
- **`stan/tumor/_lfo_endpoints_generated_quantities.stan`** — call site (LFO); unpack the 2 new values into real arrays `sample_prog_pfs`, `sample_prog_right_censored` so they are in scope for the stamp.
- **`stan/tumor/sf-ssls-lfo.stan:248-271`** — switch the RECIST stamp to use `sample_prog_pfs`/`sample_prog_right_censored`.
- **`tests/testthat/test-stan-lfo-ms-all-transition.R`** — add a behavioral test that a death-only (0→2) forecast does NOT stamp PD, while a 0→1 progression does.
- **`quarto/sclc/website/documentation/clinical-endpoints-specification.qmd`** — document that the OOS RECIST confusion matrix uses progression-only (0→1) prediction, distinct from operational PFS.

**Ordering rationale:** The Stan compiler validates the full tuple arity at every call site, so the return-tuple change (Task 1) breaks compilation until ALL call sites are updated (Tasks 2–6). Tasks 1–6 must be completed and compiled together before the stamp switch (Task 7) is meaningful. Commit after each call-site fix, but expect `stanc` to fail until Task 6 is done — each task says so explicitly.

---

## Task 1: Add progression-only return values to the endpoint RNG function

**Files:**
- Modify: `stan/pfs.stanfunctions` (return-tuple type block ~1618–1646; declarations near 2030; population near 2050; final `return (...)` near 2097)

- [ ] **Step 1: Read the three regions to edit**

Run: open `stan/pfs.stanfunctions` and locate:
- The `tuple(...)` return-type block ending at the line `array[] int   // spop_dropout_week (pre-graft cause-3 exit week; 0 otherwise)` (~1645).
- The local-array declarations block for sample outputs (search for `array[n_patients] int sample_ms_pfs;` — declarations live in the function body near the top of the per-patient loop setup).
- The per-patient loop where `sample_t01`/`sample_c01` are computed (~2042–2050).
- The final `return (` statement (~2097).

- [ ] **Step 2: Add the two new return-tuple type entries**

In the `tuple(...)` return-type block, after the last entry `array[] int   // spop_dropout_week ...`, add a trailing comma to that line and two new lines:

```stan
  array[] int,  // spop_dropout_week (pre-graft cause-3 exit week; 0 otherwise)
  array[] int,  // sample_prog_pfs (0→1 progression-only PFS; excludes 0→2 death, 0→3 dropout)
  array[] int   // sample_prog_right_censored (censored unless a 0→1 progression occurred)
```

- [ ] **Step 3: Declare the two new output arrays**

Find the declaration of `sample_is_dropout` (search `array[n_patients] int sample_is_dropout`) and add immediately after it:

```stan
  array[n_patients] int sample_prog_pfs = zeros_int_array(n_patients);
  array[n_patients] int sample_prog_right_censored = rep_array(1, n_patients);
```

If the existing declarations use a size other than `n_patients` (e.g. `n_forecast_patients`), match the size used by `sample_is_dropout` exactly — read that declaration and copy its size expression.

- [ ] **Step 4: Populate the new arrays inside the per-patient loop**

In the per-patient loop, immediately after the line `sample_is_dropout[j] = (sample_cause == 3);` (~2060), add:

```stan
    // Progression-only PFS for the OOS RECIST confusion matrix: the 0→1 event
    // (target-lesion PD OR MS 0→1 hazard), excluding 0→2 death and 0→3 dropout.
    // Captured here, before the death/dropout PFS-from-OS graft below, so it
    // reflects radiographic progression only. PFS/OS endpoints are unaffected.
    sample_prog_pfs[j]            = sample_t01;
    sample_prog_right_censored[j] = sample_c01;
```

Note: `j` is the loop's output index (same index used by `sample_is_dropout[j]`). Confirm by reading the surrounding lines that `sample_is_dropout[j]` uses `j` — if it uses a different index variable, use that same one.

- [ ] **Step 5: Add the two values to the final return statement**

Find the final `return (` of the function. It ends with `spop_dropout_week);` (or the matching last element). Add the two new arrays as the final elements:

```stan
    spop_dropout_week,
    sample_prog_pfs,
    sample_prog_right_censored);
```

(Adjust the preceding line to add a comma after `spop_dropout_week`.)

- [ ] **Step 6: Syntax-check (expected to FAIL at call sites)**

Run: `~/.cmdstan/cmdstan-2.38.0/bin/stanc --include-paths=stan --include-paths=stan/tumor stan/tumor/sf-ssm-log-space.stan`
Expected: FAIL — error about tuple arity mismatch at a `calculate_all_patients_endpoints_rng` call site (the full model unpacks 27 elements, now 29). This confirms the return type changed. Call sites are fixed in Tasks 2–6.

- [ ] **Step 7: Commit**

```bash
git add stan/pfs.stanfunctions
git commit -m "feat(lfo): return progression-only PFS from endpoint RNG for OOS RECIST

Adds sample_prog_pfs / sample_prog_right_censored (0->1 progression only,
excluding 0->2 death and 0->3 dropout) to calculate_all_patients_endpoints_rng.
Call sites updated in follow-up commits.

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 2: Update the `compute_ms_os_km_rng` canonical overload unpack

**Files:**
- Modify: `stan/multistate.stanfunctions` (~1772, the `) = calculate_all_patients_endpoints_rng(` unpack)

- [ ] **Step 1: Locate the tuple unpack**

Run: open `stan/multistate.stanfunctions`, find the assignment whose RHS is `calculate_all_patients_endpoints_rng(` (~1772). The LHS is a `(... ) =` tuple-unpack listing variable names or `_` discards.

- [ ] **Step 2: Add two discard slots to the unpack**

The function returns 29 elements now. The overload only uses `sample_os` for KM and discards the rest. Add two trailing `_` discards (or named throwaway vars matching the file's style) as the final two unpack targets, immediately after the current last element (`spop_dropout_week` or its discard).

If the file uses explicit throwaway declarations (not `_`), declare two:
```stan
  array[n_patients] int ignore_sample_prog_pfs;
  array[n_patients] int ignore_sample_prog_right_censored;
```
and append them to the unpack. If it uses `_`, append two `_`.

Read the existing unpack to see which convention is in use, and match it exactly.

- [ ] **Step 3: Syntax-check (still expected to FAIL at OTHER call sites)**

Run: `~/.cmdstan/cmdstan-2.38.0/bin/stanc --include-paths=stan --include-paths=stan/tumor stan/tumor/sf-ssm-log-space.stan`
Expected: FAIL, but the error should now point at a *different* call site (the full tumor model in `_tumor_endpoints_generated_quantities.stan`), not `multistate.stanfunctions`. This confirms this unpack is fixed.

- [ ] **Step 4: Commit**

```bash
git add stan/multistate.stanfunctions
git commit -m "fix(stan): unpack 2 new endpoint-RNG return values in compute_ms_os_km_rng

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 3: Update the full tumor model call site

**Files:**
- Modify: `stan/tumor/_tumor_endpoints_generated_quantities.stan` (~235, the `) = calculate_all_patients_endpoints_rng(` unpack at ~216–235)

- [ ] **Step 1: Locate the unpack**

Run: open `stan/tumor/_tumor_endpoints_generated_quantities.stan`, find the multi-line tuple-unpack ending at `) = calculate_all_patients_endpoints_rng(` (~216). The LHS lists ~27 variable names ending with `spop_is_dropout, sample_is_dropout` and possibly `spop_dropout_week`.

- [ ] **Step 2: Add two discard slots**

This call site does NOT need progression-only PFS (the full model's OOS RECIST is handled differently). Append two `_` discards (matching Stan 2.38 tuple-discard syntax) as the final two unpack targets. If the file's last unpacked element is `spop_dropout_week`, the unpack becomes `... spop_dropout_week, _, _) = calculate_all_patients_endpoints_rng(`.

If `spop_dropout_week` is currently discarded with `_` already, just add two more `_`. Read the actual last elements before editing — match the exact current tail.

- [ ] **Step 3: Syntax-check (full model should now PASS)**

Run: `~/.cmdstan/cmdstan-2.38.0/bin/stanc --include-paths=stan --include-paths=stan/tumor stan/tumor/sf-ssm-log-space.stan`
Expected: PASS (the full model's only call site is now consistent). If it still fails, the error names the remaining offending file — continue.

- [ ] **Step 4: Commit**

```bash
git add stan/tumor/_tumor_endpoints_generated_quantities.stan
git commit -m "fix(stan): unpack 2 new endpoint-RNG return values in tumor endpoints GQ

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 4: Update the pioneer burden_endpoints call site

**Files:**
- Modify: `stan/modules/state_space/burden_endpoints.stan` (~50, the `) = calculate_all_patients_endpoints_rng(` unpack)

- [ ] **Step 1: Locate the unpack**

Run: open `stan/modules/state_space/burden_endpoints.stan`, find the unpack ending at `) = calculate_all_patients_endpoints_rng(` (~50).

- [ ] **Step 2: Add two discard slots**

Append two `_` discards as the final two unpack targets (pioneer does not use progression-only PFS). Match the current tail exactly as in Task 3.

- [ ] **Step 3: Syntax-check the pioneer model**

Run: `~/.cmdstan/cmdstan-2.38.0/bin/stanc --include-paths=stan --include-paths=stan/psa stan/psa/pioneer.stan`
Expected: PASS for this call site (other unrelated errors, if any, are pre-existing — but this specific tuple-arity error must be gone).

- [ ] **Step 4: Commit**

```bash
git add stan/modules/state_space/burden_endpoints.stan
git commit -m "fix(stan): unpack 2 new endpoint-RNG return values in burden_endpoints

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 5: Update the ms-standalone LFO OS-KM call site

**Files:**
- Modify: `stan/_ms_standalone_lfo_os_km_generated_quantities.stan` (~161, the `) = calculate_all_patients_endpoints_rng(` unpack)

- [ ] **Step 1: Locate the unpack**

Run: open `stan/_ms_standalone_lfo_os_km_generated_quantities.stan`, find the unpack ending at `) = calculate_all_patients_endpoints_rng(` (~161).

- [ ] **Step 2: Add two discard slots**

Append two `_` discards as the final two unpack targets. Match the current tail exactly.

- [ ] **Step 3: Syntax-check (deferred to Task 6 compile)**

This file is included by ms-standalone models. Run the broad LFO compile check now:
Run: `~/.cmdstan/cmdstan-2.38.0/bin/stanc --include-paths=stan --include-paths=stan/tumor stan/tumor/sf-ssls-lfo.stan`
Expected: FAIL — but the error should now point only at `_lfo_endpoints_generated_quantities.stan` (Task 6), not this file.

- [ ] **Step 4: Commit**

```bash
git add stan/_ms_standalone_lfo_os_km_generated_quantities.stan
git commit -m "fix(stan): unpack 2 new endpoint-RNG return values in ms-standalone OS-KM GQ

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 6: Unpack progression-only PFS into real arrays in the LFO endpoints GQ

**Files:**
- Modify: `stan/tumor/_lfo_endpoints_generated_quantities.stan` (declarations near top ~13–27; unpack at ~216–230)

This is the ONE call site that actually USES the new values (the stamp in `sf-ssls-lfo.stan` reads them).

- [ ] **Step 1: Declare the receiving arrays**

Near the other endpoint output declarations (search for `array[n_cutoff_observed_patients] int<lower = 0, upper = 1> spop_is_dropout, sample_is_dropout;` ~27), add:

```stan
array[n_cutoff_observed_patients] int<lower = 0> sample_prog_pfs;
array[n_cutoff_observed_patients] int<lower = 0, upper = 1> sample_prog_right_censored;
```

Match the size expression (`n_cutoff_observed_patients`) used by the sibling arrays like `sample_ms_pfs` in this file — read line 17 (`array[n_cutoff_observed_patients] int<lower = 0> sample_ms_pfs, spop_ms_pfs;`) to confirm.

- [ ] **Step 2: Unpack the two new values**

Find the LHS tuple-unpack of `calculate_all_patients_endpoints_rng` (~216–230), ending with `spop_is_dropout, sample_is_dropout, spop_dropout_week) =`. Replace the closing of the LHS so the two new arrays are unpacked into the declared variables:

```stan
   spop_is_dropout, sample_is_dropout,
   spop_dropout_week,
   sample_prog_pfs, sample_prog_right_censored) =
    calculate_all_patients_endpoints_rng(
```

Read the exact current tail first; if `spop_dropout_week` is the current last element, append `, sample_prog_pfs, sample_prog_right_censored` after it.

- [ ] **Step 3: Syntax-check (should now PASS)**

Run: `~/.cmdstan/cmdstan-2.38.0/bin/stanc --include-paths=stan --include-paths=stan/tumor stan/tumor/sf-ssls-lfo.stan`
Expected: PASS — all call sites consistent. The stamp still uses `sample_ms_pfs` (changed in Task 7), so behavior is unchanged but it compiles.

- [ ] **Step 4: Commit**

```bash
git add stan/tumor/_lfo_endpoints_generated_quantities.stan
git commit -m "feat(lfo): unpack progression-only PFS in LFO endpoints GQ

Makes sample_prog_pfs / sample_prog_right_censored available to the OOS
RECIST stamp. Stamp switch follows.

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 7: Switch the OOS RECIST stamp to progression-only

**Files:**
- Modify: `stan/tumor/sf-ssls-lfo.stan:248-271`

- [ ] **Step 1: Read the current stamp block**

Run: open `stan/tumor/sf-ssls-lfo.stan`, read lines 248–272. The relevant lines:
```stan
        int cutoff_patient_idx = patient_to_cutoff_idx[i];
        int forecast_ms_pfs = sample_ms_pfs[cutoff_patient_idx];
        int forecast_ms_censored = sample_ms_right_censored[cutoff_patient_idx];
```

- [ ] **Step 2: Change the source arrays and the comment**

Replace those three lines (255–257) with:

```stan
        int cutoff_patient_idx = patient_to_cutoff_idx[i];
        // Confusion-matrix prediction uses PROGRESSION-ONLY PFS (0→1: target-lesion
        // PD or MS 0→1 hazard). Death (0→2) and dropout (0→3) must NOT stamp PD here
        // because the observed RECIST axis is scan-only and never records death as PD.
        // PFS/OS endpoints continue to treat death as an event (unchanged).
        int forecast_ms_pfs = sample_prog_pfs[cutoff_patient_idx];
        int forecast_ms_censored = sample_prog_right_censored[cutoff_patient_idx];
```

Leave the variable names `forecast_ms_pfs` / `forecast_ms_censored` as-is (renaming would touch lines 259–269 unnecessarily); only the right-hand-side source arrays change. The downstream logic (`if (!forecast_ms_censored)` … stamp PD from `forecast_ms_pfs` week onward) is correct as written.

- [ ] **Step 3: Syntax-check**

Run: `~/.cmdstan/cmdstan-2.38.0/bin/stanc --include-paths=stan --include-paths=stan/tumor stan/tumor/sf-ssls-lfo.stan`
Expected: PASS.

- [ ] **Step 4: Commit**

```bash
git add stan/tumor/sf-ssls-lfo.stan
git commit -m "fix(lfo): stamp OOS RECIST PD on progression only, not death/dropout

The confusion-matrix prediction previously stamped PD from sample_ms_pfs
(operational PFS = progression OR death OR dropout), while the observed
RECIST axis is scan-only and never records death as PD. This asymmetry
inflated predicted PD ~3x (53% vs 19% observed). Switch the stamp to
sample_prog_pfs (0->1 progression only). PFS/OS endpoints unchanged.

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 8: Behavioral test — death does not stamp PD, progression does

**Files:**
- Modify: `tests/testthat/test-stan-lfo-ms-all-transition.R`

- [ ] **Step 1: Read the existing test scaffold**

Run: open `tests/testthat/test-stan-lfo-ms-all-transition.R` and `tests/testthat/helper-lfo.R`. Identify the helper that compiles `sf-ssls-lfo.stan` and assembles LFO stan data (the #92 parity test uses it). Note the function names used to build a minimal eval-trial patient with a known multistate outcome.

- [ ] **Step 2: Write the failing test**

Add a test that constructs two single-patient OOS scenarios at a fixed cutoff and inspects `oos_recist`:

```r
test_that("OOS RECIST stamps PD on 0->1 progression but NOT on 0->2 death", {
  skip_on_cran()
  skip_if_not(cmdstan_available(), "cmdstan not installed")

  mod <- compile_lfo_model()  # helper from helper-lfo.R

  # Scenario A: patient whose multistate outcome is a DIRECT DEATH (0->2)
  #   with NO target-lesion progression in the forecast window.
  #   Expectation: forecast oos_recist contains NO PD attributable to the
  #   multistate stamp (only SLD-trajectory RECIST applies).
  data_death <- build_lfo_oos_fixture(
    ms_outcome = "death_02",      # 0->2 direct death, no prior progression
    sld_progresses = FALSE        # SLD trajectory stays sub-PD threshold
  )
  fit_death <- mod$sample(data = data_death, chains = 1,
                          iter_warmup = 0, iter_sampling = 1,
                          fixed_param = TRUE, seed = 1)
  oos_death <- fit_death$draws("oos_recist", format = "draws_matrix")
  expect_false(any(oos_death == 4),
    info = "0->2 death must not stamp PD in oos_recist")

  # Scenario B: patient whose multistate outcome IS a 0->1 progression.
  #   Expectation: oos_recist contains PD from the progression week onward.
  data_prog <- build_lfo_oos_fixture(
    ms_outcome = "progression_01",
    sld_progresses = FALSE        # isolate the MS-stamp channel
  )
  fit_prog <- mod$sample(data = data_prog, chains = 1,
                         iter_warmup = 0, iter_sampling = 1,
                         fixed_param = TRUE, seed = 1)
  oos_prog <- fit_prog$draws("oos_recist", format = "draws_matrix")
  expect_true(any(oos_prog == 4),
    info = "0->1 progression must stamp PD in oos_recist")
})
```

If `build_lfo_oos_fixture` and `compile_lfo_model` do not already exist in `helper-lfo.R`, add them in this step modeled on the existing #92 parity-test fixtures (read those fixtures and adapt — the parity test already builds an eval-trial patient with a controllable multistate outcome and SLD path). The fixture must:
- set `lfo_eval_trial` to the patient's trial,
- place the cutoff before the forecast window,
- for `ms_outcome = "death_02"`: set `ms_final_state`/`ms_time_02`/`ms_censored_02` so the sampled cause is direct death with `ms_censored_01 = 1` (0→1 censored),
- for `ms_outcome = "progression_01"`: set `ms_time_01`/`ms_censored_01` so a 0→1 event occurs in-window,
- for `sld_progresses = FALSE`: choose `patient_log_decrease_rate`/`patient_log_growth_rate` inits (via `fixed_param` data or fixed inits) so the forecast SLD never crosses the +20% nadir / +5mm PD threshold.

- [ ] **Step 3: Run the test to verify it FAILS on the pre-Task-7 binary**

To prove the test discriminates, temporarily check out the stamp file from before Task 7 is NOT required; instead confirm the test passes on the fixed model. But first verify the test is meaningful by running it against the current (fixed) build:

Run: `Rscript -e 'testthat::test_file("tests/testthat/test-stan-lfo-ms-all-transition.R")'`
Expected: PASS (model already fixed in Task 7). 

To confirm the test would have caught the bug, run once with the stamp reverted to `sample_ms_pfs` (git stash the Task-7 change), observe Scenario A FAILS (`oos_death` contains PD), then restore:
```bash
git stash   # temporarily revert nothing committed? Task 7 is committed — instead:
git revert --no-commit HEAD~1  # only if Task 7 is the prior commit; otherwise edit the two RHS arrays back to sample_ms_pfs by hand
```
Simpler: hand-edit `sf-ssls-lfo.stan` RHS back to `sample_ms_pfs`/`sample_ms_right_censored`, recompile, run the test, confirm Scenario A FAILS, then restore the progression-only arrays. Document the observed failure in the commit message.

- [ ] **Step 4: Run the test to verify it PASSES on the fixed model**

Run: `Rscript -e 'testthat::test_file("tests/testthat/test-stan-lfo-ms-all-transition.R")'`
Expected: PASS — both scenarios.

- [ ] **Step 5: Commit**

```bash
git add tests/testthat/test-stan-lfo-ms-all-transition.R tests/testthat/helper-lfo.R
git commit -m "test(lfo): OOS RECIST stamps PD on 0->1 progression, not 0->2 death

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 9: Update the clinical-endpoints specification

**Files:**
- Modify: `quarto/sclc/website/documentation/clinical-endpoints-specification.qmd`

- [ ] **Step 1: Find the OOS RECIST / confusion-matrix section**

Run: `grep -n "confusion\|oos_recist\|OOS RECIST\|out-of-sample RECIST\|operational PFS\|ms_pfs" quarto/sclc/website/documentation/clinical-endpoints-specification.qmd`
Read the relevant section. If no OOS-RECIST subsection exists, add one under the PFS/response section.

- [ ] **Step 2: Document the progression-only convention**

Add/adjust prose stating:

> **Out-of-sample RECIST confusion matrix.** The predicted RECIST label used in the
> OOS confusion matrix marks PD only on **radiographic progression (0→1)** — the first
> of target-lesion PD or the multistate 0→1 progression hazard. Death (0→2) and dropout
> (0→3) do **not** mark PD in this prediction, because the observed RECIST axis is
> scan-based and does not record death as a RECIST category. This is distinct from the
> operational PFS endpoint (`sample_ms_pfs`), which treats death and progression as a
> single composite event. The progression-only quantity used here is `sample_prog_pfs`
> (with `sample_prog_right_censored`), computed in `calculate_all_patients_endpoints_rng`
> from the 0→1 channel before the PFS-from-OS dropout/death graft. PFS and OS endpoints
> are unchanged and continue to treat death as an event.

- [ ] **Step 3: Commit**

```bash
git add quarto/sclc/website/documentation/clinical-endpoints-specification.qmd
git commit -m "docs(spec): OOS RECIST confusion matrix uses progression-only prediction

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 10: Full regression test + parity check

**Files:**
- None (verification only)

- [ ] **Step 1: Run the full Stan LFO test suite**

Run: `Rscript -e 'testthat::test_dir("tests/testthat", filter = "lfo")'`
Expected: PASS for all LFO tests, including the existing #92 parity test (`test-stan-lfo-ms-all-transition.R`). The parity test checks the multistate *likelihood* (`lfo_ll`), which is unaffected by this GQ-only change — it must still report the known value (`lfo_ll = -21.1202170` per project memory). If it changed, STOP — the change leaked into the likelihood.

- [ ] **Step 2: Run the multistate test suite (pioneer/full model regression)**

Run: `Rscript -e 'testthat::test_dir("tests/testthat", filter = "multistate")'`
Expected: PASS. Confirms the tuple-arity changes to shared call sites did not break the full model or pioneer.

- [ ] **Step 3: Compile all three entry models**

```bash
~/.cmdstan/cmdstan-2.38.0/bin/stanc --include-paths=stan --include-paths=stan/tumor stan/tumor/sf-ssm-log-space.stan
~/.cmdstan/cmdstan-2.38.0/bin/stanc --include-paths=stan --include-paths=stan/tumor stan/tumor/sf-ssls-lfo.stan
~/.cmdstan/cmdstan-2.38.0/bin/stanc --include-paths=stan --include-paths=stan/psa stan/psa/pioneer.stan
```
Expected: all PASS (no tuple-arity errors).

- [ ] **Step 4: Final commit (if any test-driven fixes were needed)**

Only if Steps 1–3 required fixes. Otherwise no commit.

---

## Post-Plan: Rerun

After all tasks pass, the SCLC (and CRC) LFO pipeline must be rerun to regenerate `oos_confusion_matrix_sclc` with the corrected prediction. This is a separate launch via the `pioneer-toolkit:start-job` skill (store `lfo`, command `publication_targets.sh -b lfo -G 28 -W 28 -m 'oos_confusion_matrix_sclc'`) — NOT part of this plan. Do not launch without explicit user confirmation.

**Validation after rerun:** the corrected confusion matrix should show predicted-PD marginal much closer to the observed ~19% (down from 53%), and PR recall recover toward the observed ~72%, because death-only forecasts will no longer be mislabeled PD.

---

## Self-Review Notes

- **Spec coverage:** User's two decisions are both covered — progression-only prediction for the confusion matrix (Tasks 1, 6, 7); PFS/OS unchanged (verified by Task 10 Step 1 parity check + the fact that only the GQ stamp source changes); no injected weeks (no task adds visits). Concern 3 needs no code (documented in Background).
- **Type consistency:** `sample_prog_pfs` / `sample_prog_right_censored` are the same names in Task 1 (definition), Task 6 (LFO unpack), and Task 7 (stamp use). Discard slots in Tasks 2–5 do not name them. Sizes: `n_patients` inside the function (Task 1), `n_cutoff_observed_patients` at the LFO call-site receiving arrays (Task 6) — matching the sibling `sample_ms_pfs` arrays in each scope respectively. **Worker must verify both size expressions against the sibling declarations as instructed, because the function-internal size and the call-site receiving-array size legitimately differ.**
- **Placeholder scan:** Test fixtures in Task 8 reference helpers that may not exist; Step 2 explicitly instructs building them from the existing #92 parity fixtures rather than assuming them.
- **Compilation ordering:** Tasks 1–6 intentionally leave the tree non-compiling until Task 6; each task's expected `stanc` result says so.
