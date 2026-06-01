# Spop Dropout PFS-from-OS Convention — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the spop PFS for a 0→3 dropout patient follow SCLC's operational PFS convention — PFS event at the 3→2 death week if death occurs in-window, else censored at the forecast horizon — by reusing the spop OS outcome, while keeping the CIF and sojourn KMs correct.

**Architecture:** Pure generated-quantities + downstream-R change (no MCMC re-fit). In the GQ per-patient loop, reorder the OS draw before the PFS assignment and, for `spop_cause == 3`, graft `spop_os`/`spop_os_censored` into `spop_pfs`/`spop_ms_pfs`. Retain the pre-graft dropout week as a new function output so the 3→2 sojourn KM can still be computed. Reorder the CIF classifier to test `is_dropout` first so grafted dropout-deaths stay in CIF_03. Exclude dropouts from the 1→2 sojourn KM.

**Tech Stack:** Stan (`.stanfunctions`, modular `#include` GQ), R (`targets`, `cmdstanr`), `testthat`.

**Spec:** `docs/superpowers/specs/2026-05-30-spop-dropout-pfs-convention-design.md`

**Source store (read-only, for validation):** `/mnt/data/analysis-results/karim_naguib/publication/main/_targets`

---

## File Structure

| File | Responsibility | Change |
|---|---|---|
| `stan/pfs.stanfunctions` | GQ per-patient endpoint sim + CIF classifier | Add `spop_dropout_week` output; reorder OS-before-PFS; graft for cause 3; reorder `compute_trial_cif` branches |
| `stan/tumor/_tumor_endpoints_generated_quantities.stan` | declares GQ outputs, calls the sim fn, computes sojourn KMs | Declare `spop_dropout_week`; capture it from the call; fix `spop_km_32` and `spop_km_12` sojourn blocks |
| `tests/testthat/test-spop-dropout-pfs.R` | unit coverage for the new derivation rule | Create |
| `quarto/.../documentation/clinical-endpoints-specification.qmd` | model PFS estimand documentation | Add operational-PFS note |

**Stan recompiles** (GQ signature changes). No re-fit: the posterior draws are reused; only generated quantities are recomputed via `generate_quantities` / rebuilding the `tumor_ssls_draws_endpoints` → `km_rvar` → CIF targets.

---

## Task 1: Confirm baseline — Stan compiles and current `died_off_trial` gap

**Files:** none (read-only baseline)

- [ ] **Step 1: Syntax-check the publication model as-is**

Run:
```bash
~/.cmdstan/cmdstan-2.38.0/bin/stanc --include-paths=stan --include-paths=stan/tumor stan/tumor/sf-ssm-log-space.stan
```
Expected: exit 0, no output (clean parse). This is the pre-change baseline.

- [ ] **Step 2: Record the current spop gap for the regression target**

Run (via the targets MCP `eval_expr`, fit already in store):
```r
library(posterior); library(dplyr)
fit <- targets::tar_read(tumor_ssls_res_posterior,
  store = "/mnt/data/analysis-results/karim_naguib/publication/main/_targets")
ad  <- targets::tar_read(all_analysis_data,
  store = "/mnt/data/analysis-results/karim_naguib/publication/main/_targets")
np <- nrow(ad)
m  <- as_draws_matrix(fit$draws(variables = "spop_pfs"))
cols <- match(paste0("spop_pfs[", seq_len(np), "]"), colnames(m))
ad$spop_med_pfs <- apply(m[, cols], 2, median)
ad |> filter(ms_pattern == "died_off_trial") |>
  summarise(obs = median(pfs), spop = median(spop_med_pfs))
```
Expected (baseline, pre-change): `obs = 42`, `spop ≈ 17`. Record these — Task 8 verifies `spop` rises toward 42.

---

## Task 2: Add `spop_dropout_week` output to `calculate_all_patients_endpoints_rng`

**Files:**
- Modify: `stan/pfs.stanfunctions` (return-tuple type, local declaration, return statement)

The spop dropout time (`spop_time_03`) is sampled inside the function but not returned; the 3→2 sojourn KM (Task 6) needs it after the graft overwrites `spop_pfs`. Add it as a new `array[n_patients] int` output, value = the dropout week for cause-3 patients and `0` otherwise.

- [ ] **Step 1: Locate the function's return-type tuple and the matching local-declaration block**

Run:
```bash
grep -n "spop_is_dropout" stan/pfs.stanfunctions
```
Expected: the return-type `tuple(...)` signature line listing `// spop_is_dropout` (~line 1636), the local declarations near ~1696–1742, and the final `return (...)` tuple (~1703–1711).

- [ ] **Step 2: Add `array[] int` to the return-type tuple**

In the `tuple( ... )` return type of `calculate_all_patients_endpoints_rng`, immediately after the `array[] int,  // spop_is_dropout` entry add:
```stan
    array[] int,  // spop_dropout_week (pre-graft cause-3 exit week; 0 otherwise)
```

- [ ] **Step 3: Declare the local array**

After the local declaration `array[n_patients] int spop_is_dropout;` (~line 1696) add:
```stan
  array[n_patients] int spop_dropout_week = zeros_int_array(n_patients);
```

- [ ] **Step 4: Append it to the `return (...)` tuple**

In the final `return (...)` statement, immediately after `spop_is_dropout, sample_is_dropout` add `spop_dropout_week` as the last element (matching the new return-type position — re-check ordering so it is the final element in BOTH the type and the value tuple).

- [ ] **Step 5: Syntax check**

Run:
```bash
~/.cmdstan/cmdstan-2.38.0/bin/stanc --include-paths=stan --include-paths=stan/tumor stan/tumor/sf-ssm-log-space.stan
```
Expected: clean parse (the value is still all-zero; populated in Task 3).

- [ ] **Step 6: Commit**

```bash
git add stan/pfs.stanfunctions
git commit -m "feat(gq): add spop_dropout_week output to endpoints rng (plumbing)"
```

---

## Task 3: Reorder OS-before-PFS and graft for cause 3

**Files:**
- Modify: `stan/pfs.stanfunctions` GQ per-patient loop (~1946–1981)

- [ ] **Step 1: Re-read the current loop region to confirm line numbers**

Run:
```bash
sed -n '1944,1982p' stan/pfs.stanfunctions
```
Expected: `classify_spop_exit` (~1950) → `spop_is_dropout[j] = (spop_cause == 3)` (~1956) → `derive_spop_pfs(...)` assigning `spop_pfs[j]...` (~1958) → OS draw `derive_spop_os_rng(...)` (~1975).

- [ ] **Step 2: Move the OS-draw block above the PFS-derivation block**

Cut the OS block (the `surv_12_s_i`/`surv_12_t_i`/`surv_32_i` guards plus the `(spop_os[j], spop_os_censored[j]) = derive_spop_os_rng(...)` call, ~1965–1980) and paste it so it runs immediately **after** `spop_is_dropout[j] = (spop_cause == 3);` and **before** the `derive_spop_pfs(...)` call. `derive_spop_os_rng` depends only on `spop_cause`, `spop_exit`, and survival rows — all already computed above — so the move is safe.

- [ ] **Step 3: Capture the pre-graft dropout week**

Immediately after `spop_is_dropout[j] = (spop_cause == 3);` add:
```stan
    if (spop_cause == 3) spop_dropout_week[j] = spop_time_03;
```
(`spop_time_03` is the cause-3 exit week computed at ~1942; capture it before the graft overwrites the PFS slot.)

- [ ] **Step 4: Graft OS into PFS for cause 3, after both OS and the (now-following) PFS derivation**

Immediately **after** the `derive_spop_pfs(...)` assignment block, add:
```stan
    // PFS-from-OS convention for dropouts: a 0->3 patient's PFS is its OS
    // outcome — event at the 3->2 death week if death occurs in-window, else
    // censored at the forecast horizon. Reuses spop_os so PFS == OS by
    // construction (preserves PFS <= OS). See spec 2026-05-30.
    if (spop_cause == 3) {
      spop_pfs[j]               = spop_os[j];
      spop_right_censored[j]    = spop_os_censored[j];
      spop_ms_pfs[j]            = spop_os[j];
      spop_ms_right_censored[j] = spop_os_censored[j];
    }
```

- [ ] **Step 5: Syntax check**

Run:
```bash
~/.cmdstan/cmdstan-2.38.0/bin/stanc --include-paths=stan --include-paths=stan/tumor stan/tumor/sf-ssm-log-space.stan
```
Expected: clean parse.

- [ ] **Step 6: Commit**

```bash
git add stan/pfs.stanfunctions
git commit -m "feat(gq): graft spop OS into PFS for 0->3 dropouts (PFS-from-OS convention)"
```

---

## Task 4: Mark the unreachable `derive_spop_pfs` cause-3 inner branch

**Files:**
- Modify: `stan/pfs.stanfunctions` `derive_spop_pfs` cause-3 branch (~1393–1401)

- [ ] **Step 1: Re-read the branch**

Run:
```bash
sed -n '1393,1402p' stan/pfs.stanfunctions
```
Expected: the `} else if (cause == 3) {` block with an inner `if (!c_target && t_target <= t_03)` test.

- [ ] **Step 2: Replace the misleading inner logic with an accurate comment**

Replace the cause-3 body so it returns the censored-at-dropout tuple with a correct note (the inner SLD-PD branch is unreachable: `classify_spop_exit` only returns cause 3 when no uncensored target PD precedes `t_03`):
```stan
  } else if (cause == 3) {
    // Dropout (0->3): the inner "SLD PD before dropout" case is UNREACHABLE here
    // — classify_spop_exit only returns cause 3 when no uncensored target PD
    // precedes t_03 (else progression wins the priority order). For the spop
    // path the caller overwrites this result with the OS outcome (PFS-from-OS
    // graft); this return is the censored-at-dropout fallback.
    return (t_03, 1, t_03, 1);
```

- [ ] **Step 3: Syntax check**

Run:
```bash
~/.cmdstan/cmdstan-2.38.0/bin/stanc --include-paths=stan --include-paths=stan/tumor stan/tumor/sf-ssm-log-space.stan
```
Expected: clean parse.

- [ ] **Step 4: Commit**

```bash
git add stan/pfs.stanfunctions
git commit -m "refactor(gq): document derive_spop_pfs cause-3 branch as unreachable"
```

---

## Task 5: Reorder `compute_trial_cif` to test `is_dropout` first

**Files:**
- Modify: `stan/pfs.stanfunctions` `compute_trial_cif` (~1549–1561) + docstring (~1515–1523)

After the Task 3 graft, a dropout-death has `right_censored == 0` and `pfs == os`, so under the current branch order it is counted as `cnt_02` and never reaches the `is_dropout` branch — emptying CIF_03 of the dropout-death cohort. Test `is_dropout` first.

- [ ] **Step 1: Replace the classification loop body**

Replace lines ~1549–1561 (the `for (j in 1:n) { ... }` body) with:
```stan
  for (j in 1:n) {
    if (is_dropout[j]) {
      cnt_03[pfs[j]] += 1;          // 0->3 dropout (dead or alive) — checked FIRST
    } else if (!right_censored[j]) {
      // PFS event for a non-dropout: distinguish direct death (0->2) from progression (0->1)
      if (!os_censored[j] && pfs[j] == os[j]) {
        cnt_02[pfs[j]] += 1;
      } else {
        cnt_01[pfs[j]] += 1;
      }
    }
    // else: admin-censored non-dropout — no CIF contribution
  }
```

- [ ] **Step 2: Update the docstring classification rule**

In the docstring (~1515–1523) replace the classification-rule lines with:
```stan
 * Classification rule (is_dropout checked FIRST so grafted dropout-deaths,
 * which have right_censored==0 after the PFS-from-OS graft, stay in CIF_03):
 *   0→3 dropout      : is_dropout (dead or alive), regardless of PFS event status
 *   0→2 direct death : NOT dropout AND PFS event AND OS event AND pfs == os
 *   0→1 progression  : NOT dropout AND PFS event AND NOT 0→2
 *   admin censored   : NOT dropout AND right_censored (no CIF contribution)
```

- [ ] **Step 3: Syntax check**

Run:
```bash
~/.cmdstan/cmdstan-2.38.0/bin/stanc --include-paths=stan --include-paths=stan/tumor stan/tumor/sf-ssm-log-space.stan
```
Expected: clean parse.

- [ ] **Step 4: Commit**

```bash
git add stan/pfs.stanfunctions
git commit -m "fix(gq): test is_dropout first in compute_trial_cif so dropout-deaths stay in CIF_03"
```

---

## Task 6: Capture `spop_dropout_week` at the call site and fix both sojourn KMs

**Files:**
- Modify: `stan/tumor/_tumor_endpoints_generated_quantities.stan` (output declaration, call-site tuple, sojourn blocks ~404–418 and ~436–450)

- [ ] **Step 1: Declare the new output array near the other spop endpoint declarations**

After the `spop_is_dropout` declaration (~line 23) add:
```stan
array[n_patients] int<lower = 0> spop_dropout_week;
```

- [ ] **Step 2: Add it to the destructuring of the rng call**

In the `( ... ) = calculate_all_patients_endpoints_rng(` assignment (the LHS tuple ending `... spop_is_dropout, sample_is_dropout`, ~line 234), append `spop_dropout_week` as the final element to match the new return position from Task 2.

- [ ] **Step 3: Fix the 3→2 sojourn block to use the dropout week**

In the `spop_km_32` block (~436–450), replace line ~443:
```stan
            soj[idx]  = max(1, spop_os[i] - spop_pfs[i]);
```
with:
```stan
            soj[idx]  = max(1, spop_os[i] - spop_dropout_week[i]);
```
(After the graft `spop_pfs == spop_os` for dropouts, so the old form collapses to 1; the dropout week recovers the true sojourn.)

- [ ] **Step 4: Exclude dropouts from the 1→2 sojourn block**

In the `spop_km_12` block (~404–418), change the filter at line ~410 from:
```stan
          if (!spop_ms_right_censored[i]) {
```
to:
```stan
          if (!spop_ms_right_censored[i] && !spop_is_dropout[i]) {
```
and correspondingly change the count at line ~405 from
`int n_prog = n_tr - sum(spop_ms_right_censored[tr_start:tr_end]);`
to count only non-dropout progressors:
```stan
        int n_prog = 0;
        for (i in tr_start:tr_end)
          if (!spop_ms_right_censored[i] && !spop_is_dropout[i]) n_prog += 1;
```
(Dropout-deaths now have `spop_ms_right_censored == 0`; without this they leak into the 1→2 progression-sojourn curve.)

- [ ] **Step 5: Syntax check**

Run:
```bash
~/.cmdstan/cmdstan-2.38.0/bin/stanc --include-paths=stan --include-paths=stan/tumor stan/tumor/sf-ssm-log-space.stan
```
Expected: clean parse.

- [ ] **Step 6: Commit**

```bash
git add stan/tumor/_tumor_endpoints_generated_quantities.stan
git commit -m "fix(gq): recover 3->2 sojourn from dropout week; exclude dropouts from 1->2 sojourn"
```

---

## Task 7: Unit test the derivation rule (R, against the regenerated draws)

**Files:**
- Create: `tests/testthat/test-spop-dropout-pfs.R`

This test asserts the invariants the change must satisfy, computed from the regenerated endpoint draws in the store. It is a property test (no hardcoded patient IDs).

- [ ] **Step 1: Write the test**

```r
# Invariants for the spop dropout PFS-from-OS convention (spec 2026-05-30).
# Reads the regenerated endpoint draws from the publication store.
test_that("spop dropout PFS-from-OS invariants hold", {
  store <- file.path("/mnt/data/analysis-results",
                     Sys.getenv("DOMINO_STARTING_USERNAME"),
                     "publication/main/_targets")
  skip_if_not(dir.exists(store), "publication store not present")

  fit <- targets::tar_read(tumor_ssls_res_posterior, store = store)
  ad  <- targets::tar_read(all_analysis_data, store = store)
  np  <- nrow(ad)

  pull <- function(v) {
    m <- posterior::as_draws_matrix(fit$draws(variables = v))
    m[, match(paste0(v, "[", seq_len(np), "]"), colnames(m)), drop = FALSE]
  }
  s_pfs  <- pull("spop_pfs");  s_rc  <- pull("spop_right_censored")
  s_os   <- pull("spop_os");   s_osc <- pull("spop_os_censored")
  drop   <- pull("spop_is_dropout")

  # (1) PFS <= OS pointwise (every draw, every patient)
  expect_true(all(s_pfs <= s_os),
              info = "PFS must never exceed OS")

  # (2) For dropouts: PFS == OS and PFS-censoring == OS-censoring (graft identity)
  dmask <- drop == 1L
  expect_true(all(s_pfs[dmask] == s_os[dmask]),
              info = "dropout PFS must equal OS by construction")
  expect_true(all(s_rc[dmask] == s_osc[dmask]),
              info = "dropout PFS censoring must equal OS censoring")
})
```

- [ ] **Step 2: Run it (will pass only after Task 8 regenerates the draws; before that it reads stale draws)**

Run:
```bash
Rscript -e 'testthat::test_file("tests/testthat/test-spop-dropout-pfs.R")'
```
Expected: FAIL before Task 8 (stale draws have the old censor-at-dropout PFS); PASS after Task 8 regenerates `tumor_ssls_draws_endpoints`.

- [ ] **Step 3: Commit the test**

```bash
git add tests/testthat/test-spop-dropout-pfs.R
git commit -m "test: spop dropout PFS-from-OS invariants (PFS<=OS, PFS==OS for dropouts)"
```

---

## Task 8: Regenerate downstream targets and verify the gap + CIF + sojourn

**Files:** none (pipeline regeneration + validation). No MCMC re-fit.

- [ ] **Step 1: Rebuild the endpoint/KM/CIF targets from the existing fit**

The model recompiles (GQ signature changed) but the posterior draws are reused. Rebuild downstream of the fit via the targets pipeline (use the start-job skill for the Domino run, `-D -m 'tumor_ssls_draws_endpoints_posterior'`, or rebuild locally if resources allow). Confirm `tumor_ssls_km_rvar_posterior`, `all_tumor_ssls_km_rvar`, and the CIF targets rebuild.

- [ ] **Step 2: Verify the `died_off_trial` gap closed**

Run the Task 1 Step 2 snippet again.
Expected: `spop` median PFS for `died_off_trial` rises materially from ~17 toward observed 42.

- [ ] **Step 3: Verify CIF_03 retains the dropout cohort**

```r
# per trial: cnt of is_dropout patients vs CIF_03 plateau * n
ad |> count(trial, wt = NULL)  # n per trial
# inspect spop_cif_03 final value * n_trial ≈ n_dropout per trial
```
Expected: `spop_cif_03` final value × n_trial ≈ number of `spop_is_dropout` patients per trial (CIF_03 did not empty out).

- [ ] **Step 4: Verify the 3→2 sojourn KM is not collapsed**

Inspect `spop_km_32`: it must be a genuine decreasing survival curve, not a step to 0 at t=1.

- [ ] **Step 5: Run the unit test (now passes against regenerated draws)**

```bash
Rscript -e 'testthat::test_file("tests/testthat/test-spop-dropout-pfs.R")'
```
Expected: PASS.

- [ ] **Step 6: Rebuild the website and eyeball the PFS KM + tail**

```bash
rm -rf quarto/publication/website/_freeze/analysis/
quarto render quarto/publication/website
```
Expected: clean render; `survival-analysis.html` spop PFS KM closer to observed; check the tail/censoring per spec Section 4 validation point 5.

---

## Task 9: Document the PFS estimand in the model spec

**Files:**
- Modify: `quarto/publication/website/documentation/clinical-endpoints-specification.qmd` (or the sclc equivalent if that is the canonical model-spec page)

- [ ] **Step 1: Locate the PFS-definition section**

Run:
```bash
grep -rn "progression-free survival\|PFS\|dropout\|off-trial" quarto/publication/website/documentation/clinical-endpoints-specification.qmd | head
```

- [ ] **Step 2: Add the operational-PFS note**

Add a short subsection stating: the model's PFS for off-trial deaths follows SCLC's **operational** PFS definition — death without recorded progression counts as a PFS *event* at the death week (realised via the 0→3→2 path) — and is **not** canonical RECIST/FDA PFS, which would censor a death after a long lost-to-follow-up gap at last assessment. State this is a deliberate estimand choice so the model and the observed comparator share one PFS definition.

- [ ] **Step 3: Commit**

```bash
git add quarto/publication/website/documentation/clinical-endpoints-specification.qmd
git commit -m "docs(spec): state model PFS for off-trial deaths follows operational (not canonical) PFS"
```

---

## Self-Review Notes (spec coverage)

- Spec Section 1 (reorder + graft, retain t_03, degenerate enable_ms_32==0) → Tasks 2, 3.
- Spec Section 2 (CIF branch reorder, no R mirror) → Task 5 (R mirror already deleted, commit 629d245f).
- Spec Section 3 sites 1–4 (graft, compute_trial_cif, spop_km_32, spop_km_12) → Tasks 3, 5, 6.
- Spec Section 3 derive_spop_pfs unreachable note → Task 4.
- Spec Section 4 validation (gap, CIF_03 count, PFS<=OS, sojourn intact, tail) → Tasks 7, 8.
- Spec Section 5 (estimand documentation) → Task 9.
- Sample/conditional path untouched: no task modifies `derive_sample_*` — verified by omission.
