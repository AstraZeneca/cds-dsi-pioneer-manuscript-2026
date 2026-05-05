# Flatten GQ Trial Dimension in Pioneer

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Remove the `array[n_trials]` outer dimension from all GQ endpoint output arrays in the pioneer full model, replacing trial-level aggregation with flat aggregation over all `n_forecast_patients`. Subgroup conditioning continues to use the existing `cond_*` arrays.

**Architecture:** Shared Stan files (`burden_endpoints.stan`, `sojourn_km_generated_quantities.stan`) stay unchanged since they're also used by the tumor/sclc pipeline. New PSA-specific flat versions are created in `stan/psa/` and substituted via `#include` changes in `pioneer.stan` only. The ms-standalone model is out of scope for this PR (n_trials=1 there, so the existing trial-indexed output is benign; addressed in a follow-up).

**Tech Stack:** Stan (stanc syntax checks), R/tidyverse, targets, tidybayes `spread_rvars`

---

## Context

GQ output arrays like `sample_pfs_km_est`, `sample_target_km_est`, `sample_os_km_est`, `sample_km_12`, `sample_rpfs_km_est`, and all associated quantile/n-at-risk arrays are currently typed `array[n_trials] vector[...]`. For pioneer, forecasting only targets one trial, so the trial dimension is meaningless—it produces wasted output slots (zeros for Flatiron in the Laplace variant) and awkward `[trial, t]` indexing in R. The `cond_*` arrays already handle all subgroup conditioning. Issue #68 in `cds-dsi-pioneer-pioneer-2026`.

---

## Critical Files

| File | Role | Change |
|------|------|--------|
| `stan/psa/pioneer.stan` | Full PSA model | Flatten ~15 GQ declarations; swap 2 includes |
| `stan/psa/_psa_rpfs_generated_quantities.stan` | rPFS GQ | Flatten `array[n_trials]` + collapse trial loop |
| `stan/psa/_psa_burden_endpoints.stan` | **New** | Flat aggregation replacing trial-indexed `aggregate_trial_metrics` |
| `stan/psa/_psa_sojourn_km_generated_quantities.stan` | **New** | Flat sojourn KMs (no trial loop) |
| `stan/modules/state_space/burden_endpoints.stan` | Shared — **do not modify** | Used by tumor + ms-standalone |
| `stan/modules/multistate/sojourn_km_generated_quantities.stan` | Shared — **do not modify** | Used by ms-standalone (sclc) |
| `stan/pfs.stanfunctions` | **Do not modify** | `aggregate_trial_metrics` stays for tumor model |
| `targets/pioneer_targets.R` | R pipeline | Remove `[trial, t]` indexing from KM/endpoint targets |
| `r/plot_functions.R` | Plot utilities | Add `flat_cif = FALSE` param to `plot_competing_risks_cif` |

---

## Task 1: Create `stan/psa/_psa_burden_endpoints.stan`

**Files:**
- Create: `stan/psa/_psa_burden_endpoints.stan`

This is a PSA-specific version of `stan/modules/state_space/burden_endpoints.stan`. It keeps sections 1–2 and the conditional group aggregation identical, but replaces the `aggregate_trial_metrics` call (section 3) and the CIF loop (section 5) with flat equivalents.

- [ ] **Step 1: Create `stan/psa/_psa_burden_endpoints.stan`**

```stan
// ============================================================================
// BURDEN ENDPOINTS — PSA FLAT AGGREGATION
// ============================================================================
// Flat (no trial dimension) version of modules/state_space/burden_endpoints.stan.
// Used by pioneer.stan only. All other models use burden_endpoints.stan.
//
// Contract interface: identical to burden_endpoints.stan.
// GQ output variables in scope must be flat (not array[n_trials]):
//   sample_target_orr, spop_target_orr                 — real
//   sample_target_km_est, spop_target_km_est, ...      — vector[max_all_t + 1]
//   sample_target_pfs_quant, ...                       — vector[n_pfs_quantiles]
//   sample_target_pfs_quant_exceeds_max, ...           — array[n_pfs_quantiles] int
//   sample_target_pfs_n, ...                           — vector[n_pfs_timepoints]
//   spop_cif_01, spop_cif_02, spop_cif_03              — vector[max_all_t + 1]
//   sample_cif_01, sample_cif_02, sample_cif_03        — vector[max_all_t + 1]

profile("burden_endpoints") {
  (sample_target_pfs, sample_target_right_censored,
   spop_target_pfs, spop_target_right_censored,
   spop_target_obs_cens_pfs, spop_target_obs_cens_right_censored,
   sample_ms_pfs, sample_ms_right_censored,
   spop_ms_pfs, spop_ms_right_censored,
   sample_pfs, sample_right_censored,
   spop_pfs, spop_right_censored,
   sample_target_confirmed_response, sample_target_unconfirmed_response,
   spop_target_confirmed_response, spop_target_unconfirmed_response,
   forecast_target_pfs, forecast_target_right_censored,
   sample_os, sample_os_censored,
   spop_os, spop_os_censored,
   spop_is_dropout, sample_is_dropout
  ) = calculate_all_patients_endpoints_rng(
    forecast_patient_idx,
    obs_biomarker_cat,
    rep_biomarker_cat,
    forecast_obs_biomarker_cat,
    forecast_obs_visits_pos,
    forecast_observation_interval,
    log_cond_surv_01,
    enable_ms_02,
    enable_ms_03,
    enable_ms_32,
    enable_ms_12,
    ms_time_scale_12,
    log_cond_surv_02,
    log_cond_surv_03,
    log_cond_surv_32,
    log_cond_surv_12_s,
    log_cond_surv_12_t,
    burden_pfs,
    burden_right_censored,
    burden_target_pfs,
    burden_target_right_censored,
    ms_censored_01,
    ms_final_state,
    ms_time_02,
    ms_censored_02,
    ms_censored_12,
    ms_time_12,
    ms_os_event_12,
    ms_time_03,
    ms_time_32,
    ms_censored_32,
    patient_visit_pos,
    forecast_visits_pos,
    patient_last_obs_visit,
    last_predict_visit,
    t_patient_visits,
    max_all_t,
    n_patient_screening_visits,
    burden_enable_ms_visit_gated_01,
    burden_tv_coef_01_val,
    burden_forecast_obs_log_psa,
    burden_median_log_psa_obs,
    burden_iqr_log_psa_obs
  );
}

array[n_forecast_patients] int orr_sample_response =
  burden_orr_use_confirmed_response
    ? sample_target_confirmed_response
    : sample_target_unconfirmed_response;
array[n_forecast_patients] int orr_spop_response =
  burden_orr_use_confirmed_response
    ? spop_target_confirmed_response
    : spop_target_unconfirmed_response;

// ── Flat aggregation over all forecast patients ──────────────────────────────
profile("flat_aggregate_metrics") {
  sample_target_orr = mean(orr_sample_response);
  spop_target_orr   = mean(orr_spop_response);

  sample_target_km_est        = estimate_kaplan_meier(sample_target_pfs, sample_target_right_censored, max_all_t).1;
  spop_target_km_est          = estimate_kaplan_meier(spop_target_pfs, spop_target_right_censored, max_all_t, 0).1;
  spop_target_obs_cens_km_est = estimate_kaplan_meier(spop_target_obs_cens_pfs, spop_target_obs_cens_right_censored, max_all_t, 0).1;
  sample_ms_pfs_km_est        = estimate_kaplan_meier(sample_ms_pfs, sample_ms_right_censored, max_all_t, 0).1;
  spop_ms_pfs_km_est          = estimate_kaplan_meier(spop_ms_pfs, spop_ms_right_censored, max_all_t, 0).1;
  sample_pfs_km_est           = estimate_kaplan_meier(sample_pfs, sample_right_censored, max_all_t, 0).1;
  spop_pfs_km_est             = estimate_kaplan_meier(spop_pfs, spop_right_censored, max_all_t, 0).1;
  sample_os_km_est            = estimate_kaplan_meier(sample_os, sample_os_censored, max_all_t, 0).1;
  spop_os_km_est              = estimate_kaplan_meier(spop_os, spop_os_censored, max_all_t, 0).1;

  (sample_target_pfs_quant, sample_target_pfs_quant_exceeds_max) = km_quantiles(sample_target_km_est, pfs_quantiles);
  (spop_target_pfs_quant,   spop_target_pfs_quant_exceeds_max)   = km_quantiles(spop_target_km_est,   pfs_quantiles);
  (sample_ms_pfs_quant,     sample_ms_pfs_quant_exceeds_max)     = km_quantiles(sample_ms_pfs_km_est, pfs_quantiles);
  (spop_ms_pfs_quant,       spop_ms_pfs_quant_exceeds_max)       = km_quantiles(spop_ms_pfs_km_est,   pfs_quantiles);
  (sample_pfs_quant,        sample_pfs_quant_exceeds_max)        = km_quantiles(sample_pfs_km_est,    pfs_quantiles);
  (spop_pfs_quant,          spop_pfs_quant_exceeds_max)          = km_quantiles(spop_pfs_km_est,      pfs_quantiles);
  (sample_os_quant,         sample_os_quant_exceeds_max)         = km_quantiles(sample_os_km_est,     pfs_quantiles);
  (spop_os_quant,           spop_os_quant_exceeds_max)           = km_quantiles(spop_os_km_est,       pfs_quantiles);

  for (n in 1:n_pfs_timepoints) {
    sample_target_pfs_n[n] = calc_km_pfs_n(sample_target_km_est, months_to_weeks(pfs_timepoints[n]));
    spop_target_pfs_n[n]   = calc_km_pfs_n(spop_target_km_est,   months_to_weeks(pfs_timepoints[n]));
    sample_ms_pfs_n[n]     = calc_km_pfs_n(sample_ms_pfs_km_est, months_to_weeks(pfs_timepoints[n]));
    spop_ms_pfs_n[n]       = calc_km_pfs_n(spop_ms_pfs_km_est,   months_to_weeks(pfs_timepoints[n]));
    sample_pfs_n[n]        = calc_km_pfs_n(sample_pfs_km_est,    months_to_weeks(pfs_timepoints[n]));
    spop_pfs_n[n]          = calc_km_pfs_n(spop_pfs_km_est,      months_to_weeks(pfs_timepoints[n]));
    sample_os_n[n]         = calc_km_pfs_n(sample_os_km_est,     months_to_weeks(pfs_timepoints[n]));
    spop_os_n[n]           = calc_km_pfs_n(spop_os_km_est,       months_to_weeks(pfs_timepoints[n]));
  }
}

// ── Conditional group aggregation (unchanged from burden_endpoints.stan) ─────
(cond_sample_target_orr, cond_spop_target_orr,
 cond_sample_target_km_est, cond_spop_target_km_est, cond_spop_target_obs_cens_km_est,
 cond_sample_ms_pfs_km_est, cond_spop_ms_pfs_km_est,
 cond_sample_pfs_km_est, cond_spop_pfs_km_est,
 cond_sample_target_pfs_quant, cond_spop_target_pfs_quant,
 cond_sample_target_pfs_quant_exceeds_max, cond_spop_target_pfs_quant_exceeds_max,
 cond_sample_ms_pfs_quant, cond_spop_ms_pfs_quant,
 cond_sample_ms_pfs_quant_exceeds_max, cond_spop_ms_pfs_quant_exceeds_max,
 cond_sample_pfs_quant, cond_spop_pfs_quant,
 cond_sample_pfs_quant_exceeds_max, cond_spop_pfs_quant_exceeds_max,
 cond_sample_target_pfs_n, cond_spop_target_pfs_n,
 cond_sample_ms_pfs_n, cond_spop_ms_pfs_n,
 cond_sample_pfs_n, cond_spop_pfs_n,
 cond_sample_os_km_est, cond_spop_os_km_est,
 cond_sample_os_quant, cond_spop_os_quant,
 cond_sample_os_quant_exceeds_max, cond_spop_os_quant_exceeds_max,
 cond_sample_os_n, cond_spop_os_n) =
  aggregate_conditional_group_metrics(
    orr_sample_response,
    orr_spop_response,
    sample_target_pfs,
    sample_target_right_censored,
    spop_target_pfs,
    spop_target_right_censored,
    spop_target_obs_cens_pfs,
    spop_target_obs_cens_right_censored,
    sample_ms_pfs,
    sample_ms_right_censored,
    spop_ms_pfs,
    spop_ms_right_censored,
    sample_pfs,
    sample_right_censored,
    spop_pfs,
    spop_right_censored,
    sample_os,
    sample_os_censored,
    spop_os,
    spop_os_censored,
    cond_group,
    cond_group_pos,
    max_all_t,
    pfs_quantiles,
    pfs_timepoints
  );

// ── Flat CIF (empirical subdistribution for all forecast patients) ────────────
(spop_cif_01,   spop_cif_02,   spop_cif_03) = compute_trial_cif(
  spop_pfs,   spop_right_censored,   spop_is_dropout,
  spop_os,   spop_os_censored,   max_all_t);
(sample_cif_01, sample_cif_02, sample_cif_03) = compute_trial_cif(
  sample_pfs, sample_right_censored, sample_is_dropout,
  sample_os, sample_os_censored, max_all_t);
```

- [ ] **Step 2: Commit**

```bash
git add stan/psa/_psa_burden_endpoints.stan
git commit -m "feat(pioneer): add flat burden_endpoints for PSA GQ"
```

---

## Task 2: Create `stan/psa/_psa_sojourn_km_generated_quantities.stan`

**Files:**
- Create: `stan/psa/_psa_sojourn_km_generated_quantities.stan`

Flat version of `stan/modules/multistate/sojourn_km_generated_quantities.stan`. The loop over trials is replaced by a single computation over all `n_forecast_patients`.

- [ ] **Step 1: Create `stan/psa/_psa_sojourn_km_generated_quantities.stan`**

```stan
// ============================================================================
// SOJOURN KM CURVES — PSA FLAT VERSION (no trial dimension)
// ============================================================================
// Flat version of modules/multistate/sojourn_km_generated_quantities.stan.
// Used by pioneer.stan only (ms-standalone still uses the trial-indexed version).
//
// Requires same contract as the trial-indexed version:
//   sample_ms_pfs, sample_ms_right_censored, sample_os, sample_os_censored,
//   sample_pfs, sample_is_dropout  (and spop_ variants)
//   ms_max_sojourn_t, ms_max_sojourn_t_32, estimate_kaplan_meier()

vector<lower=0, upper=1>[ms_max_sojourn_t + 1]
  sample_km_12 = ones_vector(ms_max_sojourn_t + 1),
  spop_km_12   = ones_vector(ms_max_sojourn_t + 1);
vector<lower=0, upper=1>[ms_max_sojourn_t_32 + 1]
  sample_km_32 = ones_vector(ms_max_sojourn_t_32 + 1),
  spop_km_32   = ones_vector(ms_max_sojourn_t_32 + 1);

// ── 1→2: post-progression sojourn KM ─────────────────────────────────────────
{
  int n_prog = n_forecast_patients - sum(sample_ms_right_censored);
  if (n_prog > 0) {
    array[n_prog] int soj; array[n_prog] int cens;
    int idx = 1;
    for (i in 1:n_forecast_patients) {
      if (!sample_ms_right_censored[i]) {
        soj[idx]  = max(1, sample_os[i] - sample_ms_pfs[i]);
        cens[idx] = sample_os_censored[i];
        idx += 1;
      }
    }
    sample_km_12 = estimate_kaplan_meier(soj, cens, ms_max_sojourn_t, 0).1;
  }
}
{
  int n_prog = n_forecast_patients - sum(spop_ms_right_censored);
  if (n_prog > 0) {
    array[n_prog] int soj; array[n_prog] int cens;
    int idx = 1;
    for (i in 1:n_forecast_patients) {
      if (!spop_ms_right_censored[i]) {
        soj[idx]  = max(1, spop_os[i] - spop_ms_pfs[i]);
        cens[idx] = spop_os_censored[i];
        idx += 1;
      }
    }
    spop_km_12 = estimate_kaplan_meier(soj, cens, ms_max_sojourn_t, 0).1;
  }
}

// ── 3→2: off-trial sojourn KM ─────────────────────────────────────────────────
{
  int n_drop = sum(to_array_1d(sample_is_dropout));
  if (n_drop > 0) {
    array[n_drop] int soj; array[n_drop] int cens;
    int idx = 1;
    for (i in 1:n_forecast_patients) {
      if (sample_is_dropout[i]) {
        soj[idx]  = max(1, sample_os[i] - sample_pfs[i]);
        cens[idx] = sample_os_censored[i];
        idx += 1;
      }
    }
    sample_km_32 = estimate_kaplan_meier(soj, cens, ms_max_sojourn_t_32, 0).1;
  }
}
{
  int n_drop = sum(spop_is_dropout);
  if (n_drop > 0) {
    array[n_drop] int soj; array[n_drop] int cens;
    int idx = 1;
    for (i in 1:n_forecast_patients) {
      if (spop_is_dropout[i]) {
        soj[idx]  = max(1, spop_os[i] - spop_pfs[i]);
        cens[idx] = spop_os_censored[i];
        idx += 1;
      }
    }
    spop_km_32 = estimate_kaplan_meier(soj, cens, ms_max_sojourn_t_32, 0).1;
  }
}
```

- [ ] **Step 2: Commit**

```bash
git add stan/psa/_psa_sojourn_km_generated_quantities.stan
git commit -m "feat(pioneer): add flat sojourn KM GQ for PSA model"
```

---

## Task 3: Flatten `stan/psa/_psa_rpfs_generated_quantities.stan`

**Files:**
- Modify: `stan/psa/_psa_rpfs_generated_quantities.stan`

Replace the `array[n_trials]` declarations and `for (s in 1:n_trials)` loop with flat declarations and a direct computation.

- [ ] **Step 1: Replace trial-indexed declarations and loop**

Full replacement of the file from line 18 onwards (lines 1–16 — patient-level declarations — stay identical):

```stan
// Trial-level rPFS KM curves (flat — no trial dimension)
vector<lower=0, upper=1>[max_all_t + 1] sample_rpfs_km_est, spop_rpfs_km_est;
array[n_cond_group] vector<lower=0, upper=1>[max_all_t + 1]
  cond_sample_rpfs_km_est, cond_spop_rpfs_km_est;

// rPFS quantiles — flat
vector<lower=0>[n_pfs_quantiles]
  sample_rpfs_quant = zeros_vector(n_pfs_quantiles),
  spop_rpfs_quant   = zeros_vector(n_pfs_quantiles);
array[n_pfs_quantiles] int
  sample_rpfs_quant_exceeds_max = zeros_int_array(n_pfs_quantiles),
  spop_rpfs_quant_exceeds_max   = zeros_int_array(n_pfs_quantiles);

// rPFS quantiles — conditional group-level (unchanged)
array[n_cond_group] vector<lower=0>[n_pfs_quantiles]
  cond_sample_rpfs_quant = rep_array(zeros_vector(n_pfs_quantiles), n_cond_group),
  cond_spop_rpfs_quant   = rep_array(zeros_vector(n_pfs_quantiles), n_cond_group);
array[n_cond_group, n_pfs_quantiles] int
  cond_sample_rpfs_quant_exceeds_max = rep_array(zeros_int_array(n_pfs_quantiles), n_cond_group),
  cond_spop_rpfs_quant_exceeds_max   = rep_array(zeros_int_array(n_pfs_quantiles), n_cond_group);

// rPFS at fixed timepoints
vector<lower=0, upper=1>[n_pfs_timepoints] sample_rpfs_n, spop_rpfs_n;
array[n_cond_group] vector<lower=0, upper=1>[n_pfs_timepoints] cond_sample_rpfs_n, cond_spop_rpfs_n;

// Patient-level: rPFS = min(radiographic PFS, OS), censored otherwise
for (i in 1:n_forecast_patients) {
  if (!spop_ms_right_censored[i]) {
    spop_rpfs[i] = spop_ms_pfs[i]; spop_rpfs_censored[i] = 0;
  } else if (!spop_os_censored[i]) {
    spop_rpfs[i] = spop_os[i]; spop_rpfs_censored[i] = 0;
  } else {
    spop_rpfs[i] = spop_ms_pfs[i]; spop_rpfs_censored[i] = 1;
  }

  if (!sample_ms_right_censored[i]) {
    sample_rpfs[i] = sample_ms_pfs[i]; sample_rpfs_censored[i] = 0;
  } else if (!sample_os_censored[i]) {
    sample_rpfs[i] = sample_os[i]; sample_rpfs_censored[i] = 0;
  } else {
    sample_rpfs[i] = sample_ms_pfs[i]; sample_rpfs_censored[i] = 1;
  }
}

// Flat aggregation
sample_rpfs_km_est = estimate_kaplan_meier(sample_rpfs, sample_rpfs_censored, max_all_t, 0).1;
spop_rpfs_km_est   = estimate_kaplan_meier(spop_rpfs,   spop_rpfs_censored,   max_all_t, 0).1;

(sample_rpfs_quant, sample_rpfs_quant_exceeds_max) = km_quantiles(sample_rpfs_km_est, pfs_quantiles);
(spop_rpfs_quant,   spop_rpfs_quant_exceeds_max)   = km_quantiles(spop_rpfs_km_est,   pfs_quantiles);

for (n in 1:n_pfs_timepoints) {
  sample_rpfs_n[n] = calc_km_pfs_n(sample_rpfs_km_est, months_to_weeks(pfs_timepoints[n]));
  spop_rpfs_n[n]   = calc_km_pfs_n(spop_rpfs_km_est,   months_to_weeks(pfs_timepoints[n]));
}

// Conditional group aggregation (unchanged)
for (c in 1:n_cond_group) {
  int curr_group_size = get_pos_size(cond_group_pos, c);
  if (curr_group_size > 0) {
    array[curr_group_size] int curr_group_patients =
      get_int_sub_array(cond_group, cond_group_pos, c);

    cond_sample_rpfs_km_est[c] = estimate_kaplan_meier(
      sample_rpfs[curr_group_patients], sample_rpfs_censored[curr_group_patients],
      max_all_t, 0).1;
    cond_spop_rpfs_km_est[c] = estimate_kaplan_meier(
      spop_rpfs[curr_group_patients], spop_rpfs_censored[curr_group_patients],
      max_all_t, 0).1;

    (cond_sample_rpfs_quant[c], cond_sample_rpfs_quant_exceeds_max[c]) =
      km_quantiles(cond_sample_rpfs_km_est[c], pfs_quantiles);
    (cond_spop_rpfs_quant[c], cond_spop_rpfs_quant_exceeds_max[c]) =
      km_quantiles(cond_spop_rpfs_km_est[c], pfs_quantiles);

    for (n in 1:n_pfs_timepoints) {
      cond_sample_rpfs_n[c, n] = calc_km_pfs_n(cond_sample_rpfs_km_est[c], months_to_weeks(pfs_timepoints[n]));
      cond_spop_rpfs_n[c, n]   = calc_km_pfs_n(cond_spop_rpfs_km_est[c],   months_to_weeks(pfs_timepoints[n]));
    }
  } else {
    cond_sample_rpfs_km_est[c] = zeros_vector(max_all_t + 1);
    cond_spop_rpfs_km_est[c]   = zeros_vector(max_all_t + 1);
    cond_sample_rpfs_n[c]      = zeros_vector(n_pfs_timepoints);
    cond_spop_rpfs_n[c]        = zeros_vector(n_pfs_timepoints);
  }
}
```

- [ ] **Step 2: Commit**

```bash
git add stan/psa/_psa_rpfs_generated_quantities.stan
git commit -m "feat(pioneer): flatten rPFS GQ — remove n_trials dimension"
```

---

## Task 4: Flatten GQ declarations in `stan/psa/pioneer.stan`

**Files:**
- Modify: `stan/psa/pioneer.stan`

Two changes: (a) flatten all `array[n_trials]` GQ output declarations, (b) swap the two `#include` lines.

- [ ] **Step 1: Replace GQ declarations (lines 227–287 in current file)**

Replace the existing trial-indexed block. The `array[n_forecast_patients]` patient-level declarations stay unchanged. Change only the aggregated arrays:

```stan
  // ── Flat aggregated endpoint outputs ──────────────────────────────────────
  real<lower = 0, upper = 1> sample_target_orr, spop_target_orr;
  vector<lower = 0, upper = 1>[n_cond_group] cond_sample_target_orr = zeros_vector(n_cond_group),
                                              cond_spop_target_orr   = zeros_vector(n_cond_group);

  vector<lower = 0, upper = 1>[max_all_t + 1]
    sample_target_km_est, spop_target_km_est, spop_target_obs_cens_km_est,
    sample_ms_pfs_km_est, spop_ms_pfs_km_est,
    sample_pfs_km_est, spop_pfs_km_est;
  array[n_cond_group] vector<lower = 0, upper = 1>[max_all_t + 1]
    cond_sample_target_km_est, cond_spop_target_km_est, cond_spop_target_obs_cens_km_est,
    cond_sample_ms_pfs_km_est, cond_spop_ms_pfs_km_est,
    cond_sample_pfs_km_est, cond_spop_pfs_km_est;

  vector<lower = 0, upper = 1>[n_pfs_timepoints]
    sample_target_pfs_n, spop_target_pfs_n,
    sample_ms_pfs_n, spop_ms_pfs_n,
    sample_pfs_n, spop_pfs_n;
  array[n_cond_group] vector<lower = 0, upper = 1>[n_pfs_timepoints]
    cond_sample_target_pfs_n, cond_spop_target_pfs_n,
    cond_sample_ms_pfs_n, cond_spop_ms_pfs_n,
    cond_sample_pfs_n, cond_spop_pfs_n;

  vector<lower = 0>[n_pfs_quantiles]
    sample_target_pfs_quant = zeros_vector(n_pfs_quantiles),
    spop_target_pfs_quant   = zeros_vector(n_pfs_quantiles),
    sample_ms_pfs_quant     = zeros_vector(n_pfs_quantiles),
    spop_ms_pfs_quant       = zeros_vector(n_pfs_quantiles),
    sample_pfs_quant        = zeros_vector(n_pfs_quantiles),
    spop_pfs_quant          = zeros_vector(n_pfs_quantiles);
  array[n_pfs_quantiles] int
    sample_target_pfs_quant_exceeds_max = zeros_int_array(n_pfs_quantiles),
    spop_target_pfs_quant_exceeds_max   = zeros_int_array(n_pfs_quantiles),
    sample_ms_pfs_quant_exceeds_max     = zeros_int_array(n_pfs_quantiles),
    spop_ms_pfs_quant_exceeds_max       = zeros_int_array(n_pfs_quantiles),
    sample_pfs_quant_exceeds_max        = zeros_int_array(n_pfs_quantiles),
    spop_pfs_quant_exceeds_max          = zeros_int_array(n_pfs_quantiles);
  array[n_cond_group] vector<lower = 0>[n_pfs_quantiles]
    cond_sample_target_pfs_quant = rep_array(zeros_vector(n_pfs_quantiles), n_cond_group),
    cond_spop_target_pfs_quant   = rep_array(zeros_vector(n_pfs_quantiles), n_cond_group),
    cond_sample_ms_pfs_quant     = rep_array(zeros_vector(n_pfs_quantiles), n_cond_group),
    cond_spop_ms_pfs_quant       = rep_array(zeros_vector(n_pfs_quantiles), n_cond_group),
    cond_sample_pfs_quant        = rep_array(zeros_vector(n_pfs_quantiles), n_cond_group),
    cond_spop_pfs_quant          = rep_array(zeros_vector(n_pfs_quantiles), n_cond_group);
  array[n_cond_group, n_pfs_quantiles] int
    cond_sample_target_pfs_quant_exceeds_max = rep_array(zeros_int_array(n_pfs_quantiles), n_cond_group),
    cond_spop_target_pfs_quant_exceeds_max   = rep_array(zeros_int_array(n_pfs_quantiles), n_cond_group),
    cond_sample_ms_pfs_quant_exceeds_max     = rep_array(zeros_int_array(n_pfs_quantiles), n_cond_group),
    cond_spop_ms_pfs_quant_exceeds_max       = rep_array(zeros_int_array(n_pfs_quantiles), n_cond_group),
    cond_sample_pfs_quant_exceeds_max        = rep_array(zeros_int_array(n_pfs_quantiles), n_cond_group),
    cond_spop_pfs_quant_exceeds_max          = rep_array(zeros_int_array(n_pfs_quantiles), n_cond_group);

  vector<lower = 0, upper = 1>[max_all_t + 1] sample_os_km_est, spop_os_km_est;
  array[n_cond_group] vector<lower = 0, upper = 1>[max_all_t + 1]
    cond_sample_os_km_est, cond_spop_os_km_est;
  vector<lower = 0>[n_pfs_quantiles]
    sample_os_quant = zeros_vector(n_pfs_quantiles),
    spop_os_quant   = zeros_vector(n_pfs_quantiles);
  array[n_pfs_quantiles] int
    sample_os_quant_exceeds_max = zeros_int_array(n_pfs_quantiles),
    spop_os_quant_exceeds_max   = zeros_int_array(n_pfs_quantiles);
  array[n_cond_group] vector<lower = 0>[n_pfs_quantiles]
    cond_sample_os_quant = rep_array(zeros_vector(n_pfs_quantiles), n_cond_group),
    cond_spop_os_quant   = rep_array(zeros_vector(n_pfs_quantiles), n_cond_group);
  array[n_cond_group, n_pfs_quantiles] int
    cond_sample_os_quant_exceeds_max = rep_array(zeros_int_array(n_pfs_quantiles), n_cond_group),
    cond_spop_os_quant_exceeds_max   = rep_array(zeros_int_array(n_pfs_quantiles), n_cond_group);
  vector<lower = 0, upper = 1>[n_pfs_timepoints] sample_os_n, spop_os_n;
  array[n_cond_group] vector<lower = 0, upper = 1>[n_pfs_timepoints]
    cond_sample_os_n, cond_spop_os_n;

  // CIF (flat — no trial dimension)
  vector<lower = 0, upper = 1>[max_all_t + 1]
    spop_cif_01   = zeros_vector(max_all_t + 1),
    spop_cif_02   = zeros_vector(max_all_t + 1),
    spop_cif_03   = zeros_vector(max_all_t + 1),
    sample_cif_01 = zeros_vector(max_all_t + 1),
    sample_cif_02 = zeros_vector(max_all_t + 1),
    sample_cif_03 = zeros_vector(max_all_t + 1);
```

- [ ] **Step 2: Swap the two `#include` lines in the GQ computation block**

Find and replace:
```stan
    #include "modules/state_space/burden_endpoints.stan"
```
→
```stan
    #include "_psa_burden_endpoints.stan"
```

Find and replace:
```stan
  #include "modules/multistate/sojourn_km_generated_quantities.stan"
```
→
```stan
  #include "_psa_sojourn_km_generated_quantities.stan"
```

- [ ] **Step 3: Commit**

```bash
git add stan/psa/pioneer.stan
git commit -m "feat(pioneer): flatten GQ declarations — remove n_trials dimension"
```

---

## Task 5: Stanc syntax check

**Files:** (read-only)

- [ ] **Step 1: Run stanc on pioneer.stan**

```bash
~/.cmdstan/cmdstan-2.38.0/bin/stanc \
  --include-paths=stan \
  --include-paths=stan/psa \
  stan/psa/pioneer.stan
```

Expected: `Model name=pioneer_model` with no errors.

- [ ] **Step 2: Fix any errors, re-run until clean, commit fix if needed**

Common issues:
- Type mismatch: `real` vs `vector` assignment — check that ORR assignments in `_psa_burden_endpoints.stan` use `mean(array[] int)` which returns `real`
- Undeclared variable: check that all variables referenced in the new includes match the new declarations in `pioneer.stan`

---

## Task 6: Update `targets/pioneer_targets.R` — full model KM targets

**Files:**
- Modify: `targets/pioneer_targets.R`

Remove the `[trial, t]` / `[trial, q]` / `[trial, n]` indices from all trial-level spread targets and remove the `recover_types(select(all_analysis_data, trial)) |>` pipe step from those same targets. The resulting tibbles no longer have a `trial` column.

The `cond_*` targets are unchanged — they keep `[cond_group, t]` indexing.

- [ ] **Step 1: Update KM curve targets (pioneer full model)**

Targets to update — `pioneer_pfs_km_rvar`, `pioneer_ms_pfs_km_rvar`, `pioneer_km_os_rvar`, `pioneer_rpfs_km_rvar`, `pioneer_km_12_rvar`, `pioneer_km_32_rvar`. Pattern: remove `recover_types(...) |>` and change `[trial, t]` → `[t]`.

Example — `pioneer_pfs_km_rvar`:
```r
tar_target(
  pioneer_pfs_km_rvar,
  pioneer_draws_km |>
    spread_rvars(
      sample_pfs_km_est[t],
      spop_pfs_km_est[t]
    ) |>
    mutate(fit_type = type)
),
```

Apply the same pattern to `sample_ms_pfs_km_est[t]`, `sample_os_km_est[t]`, `sample_rpfs_km_est[t]`, `sample_km_12[t]`, `sample_km_32[t]`.

- [ ] **Step 2: Update endpoint targets (quantiles, n-at-risk, ORR)**

Targets to update: `pioneer_rpfs_quant`, `pioneer_rpfs_n_rvar`, `pioneer_pfs_quant`, `pioneer_quant_os`, `pioneer_pfs_n_rvar`, `pioneer_os_n_rvar`, `pioneer_target_orr`.

For quantile targets: `[trial, q]` → `[q]`, remove `recover_types`.

Example — `pioneer_pfs_quant`:
```r
tar_target(
  pioneer_pfs_quant,
  pioneer_draws_endpoints |>
    spread_rvars(sample_pfs_quant[q], spop_pfs_quant[q]) |>
    left_join(
      enframe(all_stan_data$pfs_quantiles, name = "q", value = "quantile"),
      by = "q"
    ) |>
    mutate(fit_type = type)
),
```

For n-at-risk targets: `[trial, n]` → `[n]`, remove `recover_types`.

Example — `pioneer_pfs_n_rvar`:
```r
tar_target(
  pioneer_pfs_n_rvar,
  pioneer_draws_endpoints |>
    spread_rvars(sample_pfs_n[n], spop_pfs_n[n]) |>
    left_join(pfs_timepoints, by = "n") |>
    mutate(fit_type = type)
),
```

For ORR — `pioneer_target_orr`:
```r
tar_target(
  pioneer_target_orr,
  pioneer_draws_endpoints |>
    spread_rvars(sample_target_orr, spop_target_orr) |>
    mutate(fit_type = type)
),
```

- [ ] **Step 3: Commit**

```bash
git add targets/pioneer_targets.R
git commit -m "feat(pioneer): update R targets — remove trial dimension from KM/endpoint spread"
```

---

## Task 7: Update `r/plot_functions.R` — CIF plot for flat arrays

**Files:**
- Modify: `r/plot_functions.R:1178`

`plot_competing_risks_cif` currently builds column names as `spop_cif_01[tr,t_idx]`. After flattening, the pioneer CIF columns are `spop_cif_01[t_idx]`. The sclc version still uses `[tr,t_idx]`. Add a `flat_cif = FALSE` parameter to handle both.

- [ ] **Step 1: Add `flat_cif` parameter to `plot_competing_risks_cif`**

In the function signature, add `flat_cif = FALSE`:
```r
plot_competing_risks_cif <- function(
  stan_data,
  draws_cif,
  cif_prefix = c("spop", "sample"),
  time_step = 4L,
  x_breaks_months = seq(0, 48, by = 6),
  trials = NULL,
  causes = 1:3,
  flat_cif = FALSE
) {
```

Then replace the column-name construction inside the function. Find the line:
```r
col <- str_c(cif_prefix, "_cif_0", cause_id, "[", tr, ",", t_idx, "]")
```
Replace with:
```r
col <- if (flat_cif) {
  str_c(cif_prefix, "_cif_0", cause_id, "[", t_idx, "]")
} else {
  str_c(cif_prefix, "_cif_0", cause_id, "[", tr, ",", t_idx, "]")
}
```

When `flat_cif = TRUE`, the `trials` loop should be skipped since there's only one cohort. Wrap the `map_dfr(trials, ...)` for `model_cif` with a guard:
```r
model_cif <- if (flat_cif) {
  map_dfr(time_grid, function(t) {
    t_idx <- t + 1L
    extract_cif <- function(cause_id) {
      col <- str_c(cif_prefix, "_cif_0", cause_id, "[", t_idx, "]")
      if (!col %in% colnames(draws_mat)) return(tibble(med = NA_real_, lo = NA_real_, hi = NA_real_))
      x <- draws_mat[, col]
      q <- quantile(x, c(0.10, 0.50, 0.90))
      tibble(med = q[[2]], lo = q[[1]], hi = q[[3]])
    }
    bind_cols(
      tibble(t = t, trial = 1L),  # single pseudo-trial for downstream compat
      map(1:3, function(cause_id) {
        rename_with(extract_cif(cause_id), ~str_c("cif_0", cause_id, "_", .))
      }) |> bind_cols()
    )
  })
} else {
  map_dfr(trials, function(tr) {
    map_dfr(time_grid, function(t) {
      # ... existing code unchanged ...
    })
  })
}
```

Update `quarto/pioneer/website/analysis/competing-risks.qmd` to pass `flat_cif = TRUE`:
```r
plot_competing_risks_cif(stan_data, cif_draws, cif_prefix = "sample", causes = 1:2, flat_cif = TRUE)
```

- [ ] **Step 2: Commit**

```bash
git add r/plot_functions.R quarto/pioneer/website/analysis/competing-risks.qmd
git commit -m "feat(pioneer): support flat CIF in plot_competing_risks_cif"
```

---

## Task 8: R syntax check and tests

**Files:** (read-only)

- [ ] **Step 1: Check R syntax on modified files**

```bash
Rscript -e 'source("targets/pioneer_targets.R")'
Rscript -e 'source("r/plot_functions.R")'
```

Expected: no errors (sourcing targets files may warn about unresolved targets symbols — that's fine).

- [ ] **Step 2: Run testthat**

```bash
Rscript -e 'testthat::test_dir("tests/testthat")'
```

Expected: all tests pass. The existing hierarchy and Stan tests don't cover GQ output shape, so no test failures are expected. If there are failures, investigate before proceeding.

- [ ] **Step 3: Final commit if any fixes needed**

```bash
git add -p
git commit -m "fix(pioneer): address R syntax issues from flat GQ refactor"
```

---

## Out of Scope (Follow-up)

The `ms-standalone` model (`stan/ms-standalone.stan` → `stan/tumor/_ms_standalone_generated_quantities.stan`) is shared with the sclc pipeline. Its GQ still uses `array[n_trials]`. Since pioneer ms_standalone always has `n_trials = 1`, the existing trial-indexed output works correctly for now. Flattening it would require a PSA-specific fork of `ms-standalone.stan` (changing only the final `#include` line) plus corresponding updates to the `ms_standalone_*_rvar` targets in `pioneer_targets.R`. File a follow-up issue after this PR merges.

---

## Verification

1. `stanc` exits cleanly on `stan/psa/pioneer.stan`
2. `testthat::test_dir("tests/testthat")` — all green
3. Manual spot-check: after a model run, confirm that `posterior::variables(fit)` for `pioneer` shows `sample_pfs_km_est[t]` (1D, `t` from 1 to `max_all_t + 1`) rather than `sample_pfs_km_est[trial,t]` (2D)
4. `spread_rvars(sample_pfs_km_est[t])` in R produces a tibble with columns `t`, `sample_pfs_km_est` — no `trial` column
