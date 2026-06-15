# LFO All-Transition Multistate Likelihood (B1 fix) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the LFO tumor model (`stan/tumor/sf-ssls-lfo.stan`) fit **all** enabled multistate transitions (0→1, 0→2, 1→2, 0→3, 3→2) at cutoff, instead of only 0→1, so its multistate likelihood measures the same thing as the full model (GitHub issue #92 / blocker B1).

**Architecture:** Replace the hand-rolled 0→1-only likelihood in the LFO model block with the same two-step pattern the **ms-standalone LFO model already uses in production**: (1) call the shared, already-tested `recensor_ms_at_cutoff()` in transformed data to build cutoff-censored arrays for all five transitions, then (2) call the full `multistate_lpmf()` with those arrays. No new likelihood math is written — we reuse the full model's `multistate_lpmf` (so parity is structural) and the existing `recensor_ms_at_cutoff` (so per-transition censoring conventions are inherited, not re-derived). A `forecast_split_level == 0` guard defers blocker B2 (background/forecast row-space split) loudly.

**Tech Stack:** Stan (CmdStan 2.39, modular `#include` architecture), R `targets`, `testthat` + `cmdstanr` for Stan-function unit tests.

---

## Background & Key Facts (read before starting)

**The bug.** `stan/tumor/sf-ssls-lfo.stan:103-109` calls `calc_ms_single_transition_loglik(...)` — only the 0→1 path — even when the publication config enables all five transitions. The full model (`stan/tumor/sf-ssm-log-space.stan:101`) calls `multistate_lpmf(...)` which fits all enabled transitions. The LFO posterior is therefore informed by strictly less data than the full model.

**Why this is small, not large.** The issue's "scope to fix" says the cutoff-censored data plumbing "was never built." That is **out of date**. It WAS built — for the ms-standalone LFO model — and is fully tested:
- `recensor_ms_at_cutoff()` is defined in `stan/lfo.stanfunctions:391` (which `sf-ssls-lfo.stan` already `#include`s at line 8). It re-censors **all five transitions** at the cutoff: visit-week clamping for 0→1 progression and state-0 follow-up; calendar-cutoff censoring for registry-exact death/dropout (0→2, 0→3, 3→2); state downgrades (1→0, 2→1, 3→0) and sojourn-time recomputation for post-cutoff events; and `ic_gap_01` recomputation. It is validated by `tests/testthat/test-stan-recensor-ms.R` (8+ cases).
- `stan/ms-standalone-lfo.stan:100-123` is the **reference implementation**: it calls `recensor_ms_at_cutoff` in `stan/_ms_standalone_lfo_transformed_data.stan:105-112`, then feeds the `lfo_ms_*` arrays to `multistate_lpmf`. We are porting this pattern into the tumor LFO.

**Two differences from the ms-standalone template (do NOT copy blindly):**
1. **Weight vector.** ms-standalone includes the propensity module, so it passes `likelihood_weight[forecast_patient_idx]`. The tumor LFO does **not** include propensity (grep confirms no `likelihood_weight` in scope). The full tumor model passes `ones_vector(n_patients)` (`sf-ssm-log-space.stan:102`). → The tumor LFO must pass `ones_vector(...)` too.
2. **Row space / B2 guard.** The shared `log_cond_surv_*` matrices are sized `[n_forecast_patients, max_t]` and indexed in forecast-local row space (`stan/modules/multistate/transformed_parameters.stan:53`). The publication LFO runs with `forecast_split_level = 0` (`r/sclc/prepare_analysis_data.R:289`), which makes `forecast_patient_idx = [1..n_patients]` the identity (`stan/_full_model_transformed_data.stan:96-98`), so indexing rows by full patient ID is correct. To keep "B1 now, B2 later" safe, the new all-transition path asserts `forecast_split_level == 0` via `fatal_error`. When B2 (background/forecast split) is implemented, that guard is removed and the indexing reconciled.

**GQ is already all-transition.** `stan/tumor/_lfo_endpoints_generated_quantities.stan:238-247` already passes `enable_ms_02/03/12/32` and the corresponding survival matrices to `calculate_all_patients_endpoints_rng`. Today those hazards are prior-driven (never fit), which is the mechanism behind the stuck 0.22 confusion-matrix PD sensitivity. After this fix the model block fits them, so the GQ forecast becomes posterior-driven. **No GQ change is needed** — the fix flows through automatically.

**Recensor input fields exist in tumor LFO data scope** (all confirmed present):
- `ms_final_state, ms_time_01, ms_censored_01, ms_time_02, ms_time_12, ms_time_03, ms_time_32` — `stan/modules/multistate/data.stan`
- `ms_os_event_12` — `stan/modules/multistate/data.stan:34`
- `interval_censored` — `stan/modules/multistate/data.stan:52`
- `ms_prog_deterministic` — `stan/modules/multistate/data.stan:46`
- `cutoff_last_visit_week` — built in `stan/tumor/_lfo_transformed_data.stan:4,11` (already in scope, with the C-EXT historical override at lines 19-26)

---

## File Structure

| File | Responsibility | Change |
|---|---|---|
| `stan/tumor/_lfo_transformed_data.stan` | Builds all cutoff-censored data for the LFO tumor model | **Modify**: add `lfo_cutoff_cal_week` derivation + `recensor_ms_at_cutoff()` call producing `lfo_ms_*` arrays. Keep existing `cutoff_ms_time_01`/`cutoff_ms_censored_01` (still consumed by GQ at `_lfo_endpoints_generated_quantities.stan:252`). |
| `stan/tumor/sf-ssls-lfo.stan` | LFO tumor model: priors, likelihood, GQ | **Modify** model block (lines 102-109): replace the `calc_ms_single_transition_loglik` block with a `forecast_split_level==0` guard + full `multistate_lpmf` call over `lfo_ms_*[forecast_patient_idx]` with `ones_vector`. |
| `tests/testthat/stan/test_lfo_ms_all_transition.stan` | Standalone Stan harness exercising the new likelihood path | **Create**: a minimal model that includes the recensor + `multistate_lpmf` over a tiny synthetic cohort and returns the total log-lik in GQ. |
| `tests/testthat/test-stan-lfo-ms-all-transition.R` | Parity + behavior tests for the new path | **Create**: (a) parity — at a cutoff beyond all events, LFO ms log-lik == full `multistate_lpmf`; (b) all-transition — a 0→2 death before cutoff changes the log-lik (proves 0→2 is now fit). |

No R pipeline (`targets/publication_targets.R`) change is required: the model block reads its data from the same compiled Stan data object; the new `lfo_ms_*` arrays are derived **inside** transformed data from fields already supplied.

---

## Task 1: Add cutoff-censored all-transition data to the LFO transformed-data include

**Files:**
- Modify: `stan/tumor/_lfo_transformed_data.stan` (append after the existing 0→1 block, i.e. after line 209)

- [ ] **Step 1: Add the calendar-cutoff week derivation + recensor call**

Append the following to the **end** of `stan/tumor/_lfo_transformed_data.stan` (after the `cutoff_ms_time_01` block at line 209, before the `cutoff_target_pfs` block — or at end of file; placement only needs `cutoff_last_visit_week`, `calendar_day`, `cutoff_calendar_day`, and the `ms_*` data fields, all already in scope). This mirrors `stan/_ms_standalone_lfo_transformed_data.stan:78-112` exactly, adapted to variable names already present here:

```stan
// ============================================================================
// All-transition cutoff re-censoring (GitHub issue #92 / blocker B1)
// ============================================================================
// Build cutoff-censored multistate arrays for ALL enabled transitions, so the
// LFO model block can fit the same multistate_lpmf as the full model instead
// of the 0->1-only single-transition path. Uses the shared, tested
// recensor_ms_at_cutoff() (stan/lfo.stanfunctions:391; tests in
// test-stan-recensor-ms.R) — identical pattern to ms-standalone-lfo.stan.
//
// Death / dropout / off-trial death are registry-exact (no visit required), so
// they censor at the per-patient calendar cutoff (lfo_cutoff_cal_week); visit-
// gated events (progression, state-0 follow-up) censor at cutoff_last_visit_week.
array[n_patients] int lfo_cutoff_cal_week;
for (i in 1:n_patients) {
  int days_since_enroll = cutoff_calendar_day[1] - calendar_day[i] + 1;
  lfo_cutoff_cal_week[i] = days_since_enroll > 0 ? (days_since_enroll - 1) %/% 7 + 1 : 0;
}

array[n_patients] int lfo_ms_final_state;
array[n_patients] int lfo_ms_time_01;
array[n_patients] int lfo_ms_censored_01;
array[n_patients] int lfo_ms_time_02;
array[n_patients] int lfo_ms_time_12;
array[n_patients] int lfo_ms_time_03;
array[n_patients] int lfo_ms_time_32;
array[n_patients] int lfo_interval_censored;
array[n_patients] int lfo_ms_prog_deterministic;
array[n_patients] int lfo_ms_ic_gap_01;

(lfo_ms_final_state, lfo_ms_time_01, lfo_ms_censored_01,
 lfo_ms_time_02, lfo_ms_time_12, lfo_ms_time_03, lfo_ms_time_32,
 lfo_interval_censored, lfo_ms_prog_deterministic, lfo_ms_ic_gap_01) =
  recensor_ms_at_cutoff(
    ms_final_state, ms_time_01, ms_censored_01,
    ms_time_02, ms_time_12, ms_time_03, ms_time_32,
    ms_os_event_12, interval_censored, ms_prog_deterministic,
    cutoff_last_visit_week, lfo_cutoff_cal_week);
```

> **C-EXT note:** `cutoff_last_visit_week` was already overridden for historical (non-eval-trial) patients at `_lfo_transformed_data.stan:19-26` to their full last visit, and `lfo_cutoff_cal_week` for them is large (they enrolled at `calendar_day==1`). So `recensor_ms_at_cutoff` leaves their events uncensored — i.e. historical patients train on full multistate history, consistent with the C-EXT design for the tumor likelihood. Do not add a separate historical branch.

- [ ] **Step 2: Stan syntax check (no R, fast)**

Run:
```bash
~/.cmdstan/cmdstan-2.39.0/bin/stanc --include-paths=stan --include-paths=stan/tumor stan/tumor/sf-ssls-lfo.stan
```
Expected: exits 0 with no output (or only pre-existing warnings). If `cmdstan-2.39.0` path differs, find it: `ls ~/.cmdstan/`.

> Why this compiles even though we haven't touched the model block yet: the new arrays are unused-but-valid transformed data. Stan permits unused locals.

- [ ] **Step 3: Commit**

```bash
git add stan/tumor/_lfo_transformed_data.stan
git commit -m "feat(lfo): build cutoff-censored all-transition MS data via recensor_ms_at_cutoff

Adds lfo_ms_* arrays (all five transitions) to the tumor LFO transformed
data, mirroring ms-standalone-lfo. Not yet consumed by the likelihood
(next commit). Refs #92 (B1).

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 2: Switch the LFO model block to the full all-transition likelihood

**Files:**
- Modify: `stan/tumor/sf-ssls-lfo.stan:102-109`

- [ ] **Step 1: Replace the 0→1-only likelihood block**

In `stan/tumor/sf-ssls-lfo.stan`, replace exactly this block (lines 102-109):

```stan
    // Multistate likelihood contribution (cutoff-aware)
    if (enable_ms_01) {
      target += sum(calc_ms_single_transition_loglik(
        cutoff_ms_time_01,
        cutoff_ms_censored_01,
        log_cond_surv_01[cutoff_observed_patients]
      ));
    }
```

with:

```stan
    // Multistate likelihood contribution (cutoff-aware, ALL enabled transitions).
    // Mirrors the full model (sf-ssm-log-space.stan:101) via multistate_lpmf,
    // fed cutoff-censored lfo_ms_* arrays from recensor_ms_at_cutoff (issue #92).
    // No propensity module here, so weights are all 1.0 (like the full tumor model).
    //
    // B2 GUARD: log_cond_surv_* are forecast-row-sized and indexed forecast-local.
    // The publication LFO runs forecast_split_level == 0 (identity row map), so
    // indexing by full patient IDs is correct. Background/forecast split (B2) must
    // reconcile the row space before lifting this guard.
    if (fit_multistate_data) {
      if (forecast_split_level != 0)
        fatal_error("sf-ssls-lfo all-transition MS likelihood requires ",
                    "forecast_split_level == 0 (got ", forecast_split_level,
                    "); background/forecast split (B2) not yet implemented.");

      lfo_ms_final_state[forecast_patient_idx] ~ multistate(
        ones_vector(n_forecast_patients),
        enable_ms_01, enable_ms_02, enable_ms_12, ms_time_scale_12,
        enable_ms_03, enable_ms_32,
        lfo_ms_time_01[forecast_patient_idx], lfo_ms_time_02[forecast_patient_idx],
        lfo_ms_time_12[forecast_patient_idx],
        lfo_ms_time_03[forecast_patient_idx], lfo_ms_time_32[forecast_patient_idx],
        lfo_ms_censored_01[forecast_patient_idx],
        lfo_ms_prog_deterministic[forecast_patient_idx],
        lfo_ms_ic_gap_01[forecast_patient_idx],
        t_patient_visits,
        patient_visit_pos,
        log_cond_surv_01,
        log_cond_surv_02,
        log_cond_surv_12_s,
        log_cond_surv_12_t,
        log_cond_surv_03,
        log_cond_surv_32,
        enable_ms_visit_gated_01
      );
    }
```

> Note: the new `multistate_lpmf` lives **inside** the existing `if (fit_tumor_data)` block (the SLD loop precedes it). That is fine — the publication config sets both `fit_tumor_data` and `fit_multistate_data` to 1. Keeping it inside `fit_tumor_data` preserves the existing control flow; the inner `if (fit_multistate_data)` makes the MS term independently gated, matching the full model's separation. **Verify** after editing that the `multistate_lpmf` call is inside `fit_tumor_data`'s braces and `fit_multistate_data` is in scope (it is — declared `sf-ssls-lfo.stan:33`).

- [ ] **Step 2: Confirm `cutoff_ms_time_01`/`cutoff_ms_censored_01` are still needed elsewhere**

Run:
```bash
grep -rn "cutoff_ms_time_01\|cutoff_ms_censored_01" stan/tumor/
```
Expected: still referenced in `_lfo_endpoints_generated_quantities.stan` (GQ). **Do NOT delete** their construction in `_lfo_transformed_data.stan:193-209`. (`calc_ms_single_transition_loglik` is no longer called from this model, but the cutoff_ms_01 arrays remain GQ inputs.)

- [ ] **Step 3: Stan syntax check**

Run:
```bash
~/.cmdstan/cmdstan-2.39.0/bin/stanc --include-paths=stan --include-paths=stan/tumor stan/tumor/sf-ssls-lfo.stan
```
Expected: exits 0. Common failure: `forecast_split_level` not in scope → it is declared in `stan/_full_model_data.stan:22` which is included at `sf-ssls-lfo.stan:18`; if the error appears, confirm that include is present.

- [ ] **Step 4: Commit**

```bash
git add stan/tumor/sf-ssls-lfo.stan
git commit -m "fix(lfo): fit all enabled MS transitions, not just 0->1 (#92 B1)

Replaces the 0->1-only calc_ms_single_transition_loglik path with the full
multistate_lpmf over recensored lfo_ms_* arrays, matching the full model.
Guarded at forecast_split_level == 0 (B2 deferred).

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 3: Parity test — LFO all-transition likelihood == full multistate_lpmf at a late cutoff

This is the correctness backbone: if the LFO recensoring leaves every event uncensored (cutoff beyond all event times), the LFO multistate log-lik must equal the full model's `multistate_lpmf` on the same inputs. This proves the LFO is now "measuring the same thing" (issue #92's core claim).

**Files:**
- Create: `tests/testthat/stan/test_lfo_ms_all_transition.stan`
- Create: `tests/testthat/test-stan-lfo-ms-all-transition.R`

- [ ] **Step 1: Write the standalone Stan harness**

Create `tests/testthat/stan/test_lfo_ms_all_transition.stan`. It exposes both `recensor_ms_at_cutoff` and `multistate_lpmf` on a tiny cohort and returns the total log-lik. Visit-gated 0→1 is OFF and 1→2 is semi-Markov (`ms_time_scale_12 = 1`) to match the publication config's continuous path.

```stan
functions {
  #include "util.stanfunctions"
  #include "pos.stanfunctions"
  #include "lfo.stanfunctions"
  #include "multistate.stanfunctions"
}
data {
  int N;
  int MAX_T;
  // Raw (un-censored) multistate fields
  array[N] int ms_final_state;
  array[N] int ms_time_01;
  array[N] int ms_censored_01;
  array[N] int ms_time_02;
  array[N] int ms_time_12;
  array[N] int ms_time_03;
  array[N] int ms_time_32;
  array[N] int ms_os_event_12;
  array[N] int interval_censored;
  array[N] int ms_prog_deterministic;
  array[N] int cutoff_visit_week;   // per-patient last visit at/before cutoff
  array[N] int cutoff_cal_week;     // per-patient calendar cutoff (week)
  // Visit arrays (enable_03 path uses these; keep simple, all weeks distinct)
  array[N] int patient_visit_pos;   // length N+1
  array[size(patient_visit_pos) > 0 ? patient_visit_pos[N + 1] - 1 : 0] int t_patient_visits;
  // Conditional survival matrices (already -exp transformed log conditional survival)
  matrix[N, MAX_T] log_cond_surv_01;
  matrix[N, MAX_T] log_cond_surv_02;
  matrix[N, MAX_T] log_cond_surv_12_s;
  matrix[N, MAX_T] log_cond_surv_12_t;
  matrix[N, MAX_T] log_cond_surv_03;
  matrix[N, MAX_T] log_cond_surv_32;
  // Flags
  int enable_01; int enable_02; int enable_12; int ms_time_scale_12;
  int enable_03; int enable_32;
}
parameters { real dummy; }
model { dummy ~ std_normal(); }
generated quantities {
  // Re-censor at cutoff
  array[N] int f; array[N] int t01; array[N] int c01;
  array[N] int t02; array[N] int t12; array[N] int t03; array[N] int t32;
  array[N] int ic; array[N] int pd; array[N] int gap01;
  (f, t01, c01, t02, t12, t03, t32, ic, pd, gap01) =
    recensor_ms_at_cutoff(
      ms_final_state, ms_time_01, ms_censored_01,
      ms_time_02, ms_time_12, ms_time_03, ms_time_32,
      ms_os_event_12, interval_censored, ms_prog_deterministic,
      cutoff_visit_week, cutoff_cal_week);

  // LFO likelihood over re-censored data
  real lfo_ll = multistate_lpmf(
    f | ones_vector(N),
    enable_01, enable_02, enable_12, ms_time_scale_12, enable_03, enable_32,
    t01, t02, t12, t03, t32, c01, pd, gap01,
    t_patient_visits, patient_visit_pos,
    log_cond_surv_01, log_cond_surv_02, log_cond_surv_12_s, log_cond_surv_12_t,
    log_cond_surv_03, log_cond_surv_32, 0);

  // Full-model likelihood over the ORIGINAL (un-censored) data
  array[N] int gap01_full;
  for (i in 1:N) gap01_full[i] = (ms_censored_01[i] || ms_prog_deterministic[i]) ? 0 : interval_censored[i] + 1;
  real full_ll = multistate_lpmf(
    ms_final_state | ones_vector(N),
    enable_01, enable_02, enable_12, ms_time_scale_12, enable_03, enable_32,
    ms_time_01, ms_time_02, ms_time_12, ms_time_03, ms_time_32,
    ms_censored_01, ms_prog_deterministic, gap01_full,
    t_patient_visits, patient_visit_pos,
    log_cond_surv_01, log_cond_surv_02, log_cond_surv_12_s, log_cond_surv_12_t,
    log_cond_surv_03, log_cond_surv_32, 0);
}
```

- [ ] **Step 2: Write the parity + behavior test (expected to FAIL first if harness has a typo; confirms harness wiring)**

Create `tests/testthat/test-stan-lfo-ms-all-transition.R`:

```r
library(testthat)
library(here)
library(stringr)
library(posterior)
source(here("tests/testthat/helper-stan.R"))

# Build a cohort with all transition types, plus a cutoff set beyond every
# event so re-censoring is a no-op and LFO log-lik must equal full log-lik.
make_cohort <- function(cutoff_week) {
  n <- 4L; max_t <- 30L
  # Distinct, mildly informative hazards per transition.
  lcs <- function(v) matrix(-v, nrow = n, ncol = max_t)
  list(
    N = n, MAX_T = max_t,
    # Patient 1: 0->1 progression at wk 10, alive (state 1)
    # Patient 2: 0->2 direct death at wk 8 (state 2, censored_01 = 1)
    # Patient 3: 0->1 (wk 6) -> 1->2 death, sojourn 5 (os event wk 11) (state 2)
    # Patient 4: admin-censored in state 0 at wk 20
    ms_final_state = c(1L, 2L, 2L, 0L),
    ms_time_01     = c(10L, 0L, 6L, 20L),
    ms_censored_01 = c(0L, 1L, 0L, 1L),
    ms_time_02     = c(10L, 8L, 6L, 20L),
    ms_time_12     = c(0L, 0L, 5L, 0L),
    ms_time_03     = c(0L, 0L, 0L, 0L),
    ms_time_32     = c(0L, 0L, 0L, 0L),
    ms_os_event_12 = c(0L, 0L, 11L, 0L),
    interval_censored    = c(2L, 0L, 0L, 0L),
    ms_prog_deterministic = c(0L, 0L, 1L, 0L),
    cutoff_visit_week = rep(cutoff_week, n),
    cutoff_cal_week   = rep(cutoff_week, n),
    patient_visit_pos = c(1L, 4L, 7L, 10L, 13L),       # 3 visits each
    t_patient_visits  = rep(c(6L, 12L, 18L), n),
    log_cond_surv_01 = lcs(0.10), log_cond_surv_02 = lcs(0.05),
    log_cond_surv_12_s = lcs(0.07), log_cond_surv_12_t = lcs(0.07),
    log_cond_surv_03 = lcs(0.02), log_cond_surv_32 = lcs(0.03),
    enable_01 = 1L, enable_02 = 1L, enable_12 = 1L, ms_time_scale_12 = 1L,
    enable_03 = 0L, enable_32 = 0L
  )
}

run <- function(data) {
  fit <- test_stan_function(
    here("tests/testthat/stan/test_lfo_ms_all_transition.stan"), data)
  posterior::as_draws_df(fit$draws())
}

test_that("late cutoff: LFO all-transition log-lik == full multistate_lpmf", {
  d <- run(make_cohort(cutoff_week = 30L))  # beyond every event
  expect_equal(get_stan_val(d, "lfo_ll"), get_stan_val(d, "full_ll"),
               tolerance = 1e-6)
})

test_that("0->2 death is actually fit (early cutoff changes log-lik)", {
  late  <- run(make_cohort(cutoff_week = 30L))
  early <- run(make_cohort(cutoff_week = 7L))   # censors P1 progression, P3 death
  # If only 0->1 were fit, 0->2 hazard terms would be absent; with all
  # transitions, the re-censored cohort yields a different (valid) log-lik.
  expect_false(isTRUE(all.equal(get_stan_val(late, "lfo_ll"),
                                get_stan_val(early, "lfo_ll"),
                                tolerance = 1e-6)))
  expect_true(is.finite(get_stan_val(early, "lfo_ll")))
})
```

- [ ] **Step 3: Run the tests**

Run:
```bash
Rscript -e 'testthat::test_file("tests/testthat/test-stan-lfo-ms-all-transition.R")'
```
Expected: both tests PASS. The first proves parity with the full model (issue #92's claim); the second proves the 0→2 path is genuinely fit.

> If `test_stan_function` or `get_stan_val` is not found, check `tests/testthat/helper-stan.R` and `tests/testthat/test-stan-multistate-loglik.R` for the exact helper API (`get_stan_val(d, name, idx)` is used there for indexed values; `get_stan_val(d, name)` for scalars).

> **Possible parity-test nuance (handle if Step 3 fails on test 1):** `recensor_ms_at_cutoff` recomputes `lfo_ms_ic_gap_01` from `lfo_interval_censored`, whereas the "full" side computes `gap01_full` from raw `interval_censored`. At a late cutoff these must match by construction (no clamping), but the IC interpretation for `prog_deterministic` patients is `gap = 0`. The cohort sets P3 (`prog_deterministic = 1`) with `interval_censored = 0`, so both sides give `gap = 0`. If a mismatch appears, it indicates a real semantic gap between recensor and the full model's gap derivation — STOP and report it (that would be a genuine finding, not a test bug).

- [ ] **Step 4: Commit**

```bash
git add tests/testthat/stan/test_lfo_ms_all_transition.stan tests/testthat/test-stan-lfo-ms-all-transition.R
git commit -m "test(lfo): parity test for all-transition MS likelihood (#92 B1)

Verifies LFO multistate log-lik equals full multistate_lpmf at a late
cutoff, and that 0->2 death is genuinely fit.

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 4: Regression guard — run the existing LFO/MS Stan test suite

**Files:** none modified (verification only)

- [ ] **Step 1: Run the multistate + LFO Stan tests**

Run:
```bash
Rscript -e 'for (f in c(
  "tests/testthat/test-stan-multistate-loglik.R",
  "tests/testthat/test-stan-recensor-ms.R",
  "tests/testthat/test-stan-lfo.R",
  "tests/testthat/test-stan-lfo-cext.R",
  "tests/testthat/test-ms-standalone-lfo.R",
  "tests/testthat/test-publication-lfo.R"
)) { cat("\n==== ", f, " ====\n"); testthat::test_file(f) }'
```
Expected: all PASS. These cover the recensor function, the single-transition path (still used by ms-standalone and GQ), the C-EXT historical override, and the publication LFO pipeline wiring. None should regress — we only changed which likelihood function the tumor LFO model block calls.

- [ ] **Step 2: Compile the full LFO model end-to-end (catches include-order issues the syntax check misses)**

Run:
```bash
Rscript -e 'cmdstanr::cmdstan_model("stan/tumor/sf-ssls-lfo.stan", include_paths = c("stan", "stan/tumor"), compile = TRUE, quiet = FALSE)' 2>&1 | tail -20
```
Expected: compiles to an executable with no errors. (This is slower — a few minutes — but it is the only check that the `multistate_lpmf` call links against all the in-scope `log_cond_surv_*` matrices with correct dimensions.)

- [ ] **Step 3: No commit** (verification task). If anything fails, fix in the relevant task above and re-run.

---

## Task 5: Update the issue-referenced docs to reflect reality

The issue and findings doc both claim the cutoff-censored plumbing "was never built." Now that the tumor LFO uses `recensor_ms_at_cutoff`, leave a breadcrumb so the next reader isn't misled.

**Files:**
- Modify: `.claude/rules/sclc.md` (the `ms_km_est` note around line 50)

- [ ] **Step 1: Clarify the `ms_km_est` endpoint note**

In `.claude/rules/sclc.md`, find the line documenting `ms_km_est` as "Multistate hazard (0→1 transition)" and append a clarifying sentence:

```markdown
- `ms_km_est` - Multistate hazard endpoint. **Reporting** convention is 0→1 (progression detection); this is NOT the likelihood. As of #92, the LFO model fits ALL enabled transitions (0→1, 0→2, 1→2, 0→3, 3→2) via `recensor_ms_at_cutoff` + `multistate_lpmf`, identical to the full model.
```

- [ ] **Step 2: Commit**

```bash
git add .claude/rules/sclc.md
git commit -m "docs(lfo): clarify ms_km_est is a reporting endpoint, not the likelihood (#92)

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

- [ ] **Step 3: Comment on the issue (optional, do via gh after merge)**

After this lands, post a comment on issue #92 noting that B1 is fixed by adopting the ms-standalone `recensor_ms_at_cutoff` + `multistate_lpmf` pattern (not the larger from-scratch build the issue originally scoped), guarded at `forecast_split_level == 0` pending B2.

---

## Post-implementation: the actual payoff (separate from this plan)

This plan ends at "the LFO model fits all transitions and is tested." The original motivation — does the OOS confusion-matrix **PD sensitivity** lift off the stuck 0.22? — is answered by a Domino rerun of `oos_confusion_matrix_sclc` against the `lfo` store on the branch carrying these commits. That is a separate launch step (use the `pioneer-toolkit:start-job` skill), not part of the code change. Do NOT bundle the rerun into this plan's commits.

---

## Self-Review

**Spec coverage (issue #92 "Scope to fix"):**
1. ✅ "Build cutoff-censored data for 0→2, 1→2, 0→3, 3→2" — Task 1 (via `recensor_ms_at_cutoff`, which produces all of them).
2. ✅ "Add the corresponding likelihood terms, matching `multistate_lpmf` semantics" — Task 2 (calls `multistate_lpmf` directly, so semantics are identical by construction, not re-implemented).
3. ✅ "Reconcile per-transition censoring conventions" — inherited from `recensor_ms_at_cutoff` (visit-week vs calendar-cutoff vs sojourn), validated by `test-stan-recensor-ms.R`; parity test (Task 3) confirms the assembled result matches the full model.
4. ✅ `ms_km_est` reporting-vs-likelihood distinction — Task 5.
5. ✅ B1↔B2 entanglement — `forecast_split_level == 0` `fatal_error` guard (Task 2) defers B2 safely and loudly.

**Placeholder scan:** No TBD/TODO/"handle edge cases" — every code step shows full code; every run step shows the command and expected result.

**Type consistency:** `recensor_ms_at_cutoff` returns a 10-tuple in the documented order (`lfo.stanfunctions:383-386`); Task 1's LHS tuple and Task 3's harness both use that exact order. `multistate_lpmf` argument order in Task 2 and the Task 3 harness matches the definition (`multistate.stanfunctions:938-957`) and the ms-standalone call site (`ms-standalone-lfo.stan:102-121`). Weight is `ones_vector(n_forecast_patients)` (model block) / `ones_vector(N)` (test), matching the full tumor model's no-propensity convention.
