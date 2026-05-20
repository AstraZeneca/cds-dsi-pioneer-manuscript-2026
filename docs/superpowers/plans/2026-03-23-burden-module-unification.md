# Burden Module Unification Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

## Background and Motivation

### What are "burden endpoints"?

The PIONEER models (tumor `sf-ssm-log-space.stan` and PSA `pioneer.stan`) share the same endpoint generation logic: given a simulated biomarker trajectory, classify the response (RECIST for SLD/tumor, PCWG3 for PSA), then compute PFS/OS/KM curves, aggregation, and competing-risks CIF. This logic currently lives in two near-identical files:

- `stan/tumor/_tumor_endpoints_generated_quantities.stan` (~404 lines) — tumor/SLD version
- `stan/psa/_psa_ms_endpoints_generated_quantities.stan` (~253 lines) — PSA version

The PSA file was written for an older function signature and **does not compile** on this branch (missing `forecast_obs_*` args added by the `karim/data-in` merge).

### How we got here: the merge situation

This branch (`karim/burden-endpoints`, branched from `karim/trial-data`) is the result of integrating two feature branches:

1. **`karim/data-in` (core)** — brought: classify-first multistate architecture, assessment-visit RECIST (noise only at 6-week assessment visits, not weekly), `forecast_obs_visits_pos` / `forecast_observation_interval` data fields, `ms_os_event_12` for OS endpoint, CIF refactor
2. **`karim/laplace-for-rwd`** — brought: dual-indexing pattern (`j` = HMC loop index, `p = patient_idx[j]` = unified patient data index), Laplace IC marginalization module, pioneer website, updated PSA endpoint caller

The auto-merge broke `calculate_all_patients_endpoints_rng` (body used `i` directly instead of the `j`/`p` dual-index pattern). That was fixed in commit `064015f`. The PSA endpoint caller still doesn't compile because it was never updated for the data-in function signature.

### Key architectural concepts for this codebase

**Stan `#include` is text substitution**: There is no namespacing. A file included in the GQ block sees all variables declared before it in the same block. This is why the "contract" pattern works — the model file declares interface variables, then includes the shared computation file.

**Dual-indexing (j/p pattern)**: In generated quantities loops, `j` is the HMC-local loop counter (1 to n_hmc_patients), `p = patient_idx[j]` is the unified patient data index. Output arrays use `j`; data arrays use `p`. This pattern lets the model run HMC on a subset of patients while maintaining correct data lookups.

**Assessment-visit noise**: Measurement noise only occurs at actual assessment visits (every 6 weeks). The weekly forecast trajectory (`forecast_patient_log_*`) is deterministic mean — only `forecast_obs_*` arrays (subsampled at assessment intervals) get fresh Student-t noise. This is critical: endpoint computation must use `forecast_obs_*` not the weekly trajectory.

**RECIST vs PCWG3**: Same 1–4 encoding (CR/PR/SD/PD = Undetectable/PSA50/Stable/PSA-PD), but completely different algorithms. RECIST is a single vectorized call (`calculate_all_patients_recist`). PCWG3 requires per-patient loops with absolute PSA values, study days, and nadir tracking — this is why the categorization files cannot be merged into one.

---

**Goal:** Eliminate duplication between tumor and PSA endpoint generated quantities by creating a shared burden endpoints computation file that both models include.

**Architecture:** Stan `#include` is text substitution. Each model's GQ block declares its own output variables, computes biomarker-specific categories (RECIST/PCWG3) via a categorization include, then opens a local `{ }` block that sets up contract interface variables and includes the shared computation file. Contract variables are local (not saved to CSV). Output variables are at GQ scope (saved). The visits module (assessment-visit scheduling) moves from tumor-only to shared.

**Tech Stack:** Stan 2.38, CmdStan stanc compiler

**Branch:** `karim/burden-endpoints` (branched from `karim/trial-data` at commit `614003c`)

---

## Current State (Problems)

1. `_tumor_endpoints_generated_quantities.stan` (404 lines) contains trajectory generation that duplicates `modules/state_space/generated_quantities.stan`
2. `_psa_ms_endpoints_generated_quantities.stan` (253 lines) is ~95% identical to the tumor file
3. PSA caller is missing 3 args (`forecast_obs_*`) and `ms_os_event_12` — doesn't compile
4. PSA caller has wrong arg alignment (written for laplace-for-rwd's signature which had `interval_censored` as separate param; data-in removed it)
5. Visits module (`forecast_obs_visits_pos`, `forecast_observation_interval`) only included by tumor models, not pioneer

## Target State

Each model's GQ block follows this pattern:

```
  #include "modules/state_space/generated_quantities.stan"   // shared trajectory generation
  #include "_MODEL_categorization.stan"                       // biomarker categorization (GQ scope)

  <endpoint output variable declarations>                     // GQ scope — saved to CSV

  {
    <contract interface variables = model-specific names>     // local — NOT saved to CSV

    #include "modules/state_space/burden_endpoints.stan"      // shared computation
  }
```

Concrete:

```
tumor/sf-ssm-log-space.stan GQ block:
  #include "modules/state_space/generated_quantities.stan"
  #include "_tumor_categorization.stan"
  <declare endpoint output vars>
  { <contract: recist, pfs, ...>  #include "modules/state_space/burden_endpoints.stan" }
  #include "modules/tumor/generated_quantities.stan"

psa/pioneer.stan GQ block:
  #include "modules/state_space/generated_quantities.stan"
  #include "_psa_categorization.stan"
  <declare endpoint output vars>
  { <contract: pcwg3, psa_pfs, ...>  #include "modules/state_space/burden_endpoints.stan" }
```

## Contract Interface

The shared `burden_endpoints.stan` expects these variables in scope (set in the local `{ }` block by each model before the include):

```stan
// Biomarker categories (int arrays, 1-4 encoding: CR/PR/SD/PD or Undetectable/PSA50/Stable/PSA-PD)
array[n_total_visits] int obs_biomarker_cat;
array[n_total_visits] int rep_biomarker_cat;
array[n_total_forecast_obs_visits] int forecast_obs_biomarker_cat;

// Burden-specific PFS data (FULL unified arrays — function uses patient_idx[j] to index)
array[n_patients] int burden_pfs;
array[n_patients] int burden_right_censored;
array[n_patients] int burden_target_pfs;
array[n_patients] int burden_target_right_censored;
```

Everything else the shared file needs (`log_cond_surv_*`, `ms_*` fields, `patient_visit_pos`, etc.) comes from modules already in scope for both models.

## File Map

| Action | File | Purpose |
|--------|------|---------|
| Create | `stan/modules/state_space/burden_endpoints.stan` | Shared computation: endpoint RNG + aggregation + CIF |
| Create | `stan/tumor/_tumor_categorization.stan` | SLD aliases, weekly RECIST, assessment-visit RECIST |
| Create | `stan/psa/_psa_categorization.stan` | PSA aliases, weekly PCWG3, assessment-visit PCWG3 |
| Modify | `stan/tumor/sf-ssm-log-space.stan` | GQ: endpoint decls + contract block + new includes |
| Modify | `stan/psa/pioneer.stan` | Add visits module, GQ: endpoint decls + contract block + new includes |
| Delete | `stan/psa/_psa_ms_endpoints_generated_quantities.stan` | Replaced by shared burden_endpoints.stan |
| Keep   | `stan/tumor/_tumor_endpoints_generated_quantities.stan` | Keep for now (LFO still references patterns from it); delete later |
| Keep   | `stan/psa/_psa_endpoints_generated_quantities.stan` | Keep for now; `_psa_categorization.stan` subsumes it |
| No change | `stan/tumor/sf-ssls-lfo.stan` | LFO stays as-is (different sizing, separate concern) |

---

### Task 1: Create `_tumor_categorization.stan`

Extract from `_tumor_endpoints_generated_quantities.stan`: SLD aliases, weekly RECIST, assessment-visit RECIST computation. No contract variables here — those go in the model file.

**Files:**
- Create: `stan/tumor/_tumor_categorization.stan`
- Reference: `stan/tumor/_tumor_endpoints_generated_quantities.stan` (lines 3-7, 176-221)

- [ ] **Step 1: Create the file**

```stan
// ============================================================================
// TUMOR/SLD CATEGORIZATION
// ============================================================================
// Converts generic state-space trajectories to SLD-specific quantities and
// RECIST categories. Included at GQ scope — all variables here are saved.
//
// Requires (from state_space/generated_quantities.stan):
//   rep_patient_log_obs, rep_mean_patient_log_obs
//   forecast_patient_log_obs, forecast_mean_patient_log_obs

// --- SLD aliases (for R targets backward compatibility) ---
vector[n_total_visits] rep_patient_log_sld = rep_patient_log_obs;
vector[n_total_visits] rep_mean_patient_log_sld = rep_mean_patient_log_obs;
vector[n_total_forecast_visits] forecast_patient_log_sld = forecast_patient_log_obs;
vector[n_total_forecast_visits] forecast_mean_patient_log_sld = forecast_mean_patient_log_obs;

// --- Weekly RECIST from noisy SLD (for visualization/diagnostics) ---
array[n_total_visits] int<lower = CR, upper = PD + 1> rep_recist = rep_array(PD + 1, n_total_visits);
array[n_total_forecast_visits] int<lower = CR, upper = PD> forecast_recist;

(rep_recist, forecast_recist) = calculate_all_patients_recist(
  rep_patient_log_sld,
  forecast_patient_log_sld,
  patient_visit_pos,
  forecast_visits_pos,
  n_patient_screening_visits
);

// --- Assessment-visit noisy SLD and RECIST (for endpoint computation) ---
// These are intermediates used to build the contract variable forecast_obs_biomarker_cat.
// forecast_obs_log_sld is at GQ scope but not extracted by R code.
vector[n_total_forecast_obs_visits] forecast_obs_log_sld;
array[n_total_forecast_obs_visits] int forecast_obs_recist;

for (i in 1:n_patients) {
  int n_assessment = get_pos_size(forecast_obs_visits_pos, i);
  if (n_assessment > 0) {
    int assess_start, assess_end;
    (assess_start, assess_end) = get_pos(forecast_obs_visits_pos, i);

    int forecast_visit_start, forecast_visit_end;
    (forecast_visit_start, forecast_visit_end) = get_pos(forecast_visits_pos, i);
    int forecast_size = forecast_visit_end - forecast_visit_start + 1;

    // Subsample mean SLD at assessment weeks
    vector[n_assessment] obs_visit_mean;
    for (a in 1:n_assessment) {
      int forecast_idx = min(a * forecast_observation_interval, forecast_size);
      obs_visit_mean[a] = forecast_mean_patient_log_sld[forecast_visit_start + forecast_idx - 1];
    }

    // Apply measurement noise at assessment visits only
    forecast_obs_log_sld[assess_start:assess_end] =
      to_vector(student_t_rng(measure_nu_obs, obs_visit_mean, measure_sd_obs));
  }
}

// Assessment-visit RECIST from noisy SLD
{
  array[n_total_visits] int unused_rep_recist;
  (unused_rep_recist, forecast_obs_recist) = calculate_all_patients_recist(
    rep_patient_log_sld,
    forecast_obs_log_sld,
    patient_visit_pos,
    forecast_obs_visits_pos,
    n_patient_screening_visits
  );
}
```

- [ ] **Step 2: Verify stanc compiles tumor model with this file**

This step happens after Task 3 and Task 4 (need the shared endpoints file and model update first).

---

### Task 2: Create `_psa_categorization.stan`

Parallel to tumor: PSA aliases, weekly PCWG3, assessment-visit PCWG3. No contract variables — those go in the model file.

**Files:**
- Create: `stan/psa/_psa_categorization.stan`
- Reference: `stan/psa/_psa_endpoints_generated_quantities.stan` (all), `stan/psa/_psa_ms_endpoints_generated_quantities.stan` (assessment-visit section)

- [ ] **Step 1: Create the file**

```stan
// ============================================================================
// PSA/PCWG3 CATEGORIZATION
// ============================================================================
// Converts generic state-space trajectories to PSA-specific quantities and
// PCWG3 categories. Included at GQ scope — all variables here are saved.
//
// Requires (from state_space/generated_quantities.stan):
//   rep_patient_log_obs, rep_mean_patient_log_obs
//   forecast_patient_log_obs, forecast_mean_patient_log_obs

// --- PSA aliases (for R targets backward compatibility) ---
vector[n_total_visits] rep_patient_log_psa = rep_patient_log_obs;
vector[n_total_visits] rep_mean_patient_log_psa = rep_mean_patient_log_obs;
vector[n_total_forecast_visits] forecast_patient_log_psa = forecast_patient_log_obs;
vector[n_total_forecast_visits] forecast_mean_patient_log_psa = forecast_mean_patient_log_obs;

// --- Weekly PCWG3 from replicated PSA (for visualization/diagnostics) ---
array[n_total_visits] int<lower=1, upper=5> rep_pcwg3 = rep_array(5, n_total_visits);
array[n_total_forecast_visits] int<lower=1, upper=4> forecast_pcwg3;

// --- Assessment-visit noisy PSA and PCWG3 (for endpoint computation) ---
array[n_total_forecast_obs_visits] int forecast_obs_pcwg3;

profile("psa_categorization") {
  for (i in 1:n_patients) {
    int visit_start, visit_end;
    (visit_start, visit_end) = get_pos(patient_visit_pos, i);
    int visit_size = get_pos_size(patient_visit_pos, i);
    int n_screening = n_patient_screening_visits[i];
    int n_treat = visit_size - n_screening;

    // Convert replicated log(PSA) means to absolute PSA for PCWG3 evaluation
    vector[visit_size] rep_psa_abs = exp(rep_mean_patient_log_obs[visit_start:visit_end]);

    // Replicated PCWG3 (observed visit weeks)
    if (n_treat > 0) {
      array[n_treat] int rep_cats = calculate_psa_category(
        rep_psa_abs,
        t_patient_visits_day[visit_start:visit_end],
        n_screening,
        psa_undetectable_threshold
      );
      rep_pcwg3[(visit_start + n_screening):visit_end] = rep_cats;
    }

    // Weekly forecast PCWG3 (for visualization)
    int forecast_start, forecast_end;
    (forecast_start, forecast_end) = get_pos(forecast_visits_pos, i);
    int forecast_size = get_pos_size(forecast_visits_pos, i);

    if (forecast_size > 0) {
      // Concatenate last observed + forecast for nadir tracking continuity
      vector[1 + forecast_size] combined_psa;
      combined_psa[1] = rep_psa_abs[visit_size];
      combined_psa[2:] = exp(forecast_mean_patient_log_obs[forecast_start:forecast_end]);

      array[1 + forecast_size] int combined_days = linspaced_int_array(
        1 + forecast_size,
        (patient_last_obs_visit[i] - 1) * 7 + 4,
        (last_predict_visit - 1) * 7 + 4);

      real obs_nadir = n_treat > 0
        ? min(rep_psa_abs[(n_screening + 1):visit_size])
        : rep_psa_abs[visit_size];

      array[forecast_size] int forecast_cats = calculate_psa_category(
        combined_psa, combined_days, obs_nadir, 1, psa_undetectable_threshold
      );
      forecast_pcwg3[forecast_start:forecast_end] = forecast_cats;
    }

    // Assessment-visit forecast PCWG3 (for endpoint computation)
    int n_assessment = get_pos_size(forecast_obs_visits_pos, i);
    if (n_assessment > 0) {
      int assess_start, assess_end;
      (assess_start, assess_end) = get_pos(forecast_obs_visits_pos, i);

      // Subsample mean log PSA at assessment weeks, apply fresh noise
      vector[n_assessment] obs_visit_mean;
      for (a in 1:n_assessment) {
        int forecast_idx = min(a * forecast_observation_interval, forecast_size);
        obs_visit_mean[a] = forecast_mean_patient_log_obs[forecast_start + forecast_idx - 1];
      }
      vector[n_assessment] noisy_log_psa =
        to_vector(student_t_rng(measure_nu_obs, obs_visit_mean, measure_sd_obs));

      // Convert noisy log PSA to absolute for PCWG3
      vector[1 + n_assessment] combined_obs_psa;
      combined_obs_psa[1] = rep_psa_abs[visit_size];
      combined_obs_psa[2:] = exp(noisy_log_psa);

      array[1 + n_assessment] int combined_obs_days;
      combined_obs_days[1] = (patient_last_obs_visit[i] - 1) * 7 + 4;
      for (a in 1:n_assessment) {
        int week = min(patient_last_obs_visit[i] + a * forecast_observation_interval, last_predict_visit);
        combined_obs_days[a + 1] = (week - 1) * 7 + 4;
      }

      real obs_nadir_for_assess = n_treat > 0
        ? min(rep_psa_abs[(n_screening + 1):visit_size])
        : rep_psa_abs[visit_size];

      array[n_assessment] int obs_cats = calculate_psa_category(
        combined_obs_psa, combined_obs_days, obs_nadir_for_assess, 1, psa_undetectable_threshold
      );
      forecast_obs_pcwg3[assess_start:assess_end] = obs_cats;
    }
  }
}
```

Note: PSA does not have separate `target_pfs` vs `pfs` — for PSA, target=overall. The `burden_target_*` aliases point to the same arrays as `burden_*`. This is handled in the model file's contract block (Task 5).

---

### Task 3: Create `modules/state_space/burden_endpoints.stan`

The shared computation file. Contains ONLY the endpoint function call, aggregation calls, and CIF computation. No variable declarations — those are in each model's GQ block. Included inside a local `{ }` block after contract variables are set.

**Files:**
- Create: `stan/modules/state_space/burden_endpoints.stan`
- Reference: `stan/tumor/_tumor_endpoints_generated_quantities.stan` (lines 223-403)

- [ ] **Step 1: Create the shared file**

The file assumes in scope:
- Contract variables (set in the local block before this include): `obs_biomarker_cat`, `rep_biomarker_cat`, `forecast_obs_biomarker_cat`, `burden_pfs`, `burden_right_censored`, `burden_target_pfs`, `burden_target_right_censored`
- Endpoint output variables (declared at GQ scope before the block): `sample_target_pfs`, `spop_km_est`, etc.
- All multistate fields (`ms_*`, `log_cond_surv_*`) from multistate module
- `hmc_patient_idx`, `patient_visit_pos`, `forecast_visits_pos`, `forecast_obs_visits_pos`, etc.

```stan
// ============================================================================
// BURDEN ENDPOINTS — SHARED COMPUTATION
// ============================================================================
// Included inside a local { } block. Contract interface variables and endpoint
// output variables must already be in scope.
//
// Contract interface (local, set by model before this include):
//   obs_biomarker_cat, rep_biomarker_cat, forecast_obs_biomarker_cat
//   burden_pfs, burden_right_censored, burden_target_pfs, burden_target_right_censored
//
// Output variables (GQ scope, declared by model before the { } block):
//   sample_target_pfs, spop_target_pfs, ..., spop_cif_01, etc.

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
    hmc_patient_idx,
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
    n_patient_screening_visits
  );
}

// Aggregate to trial-level metrics
(sample_target_orr, spop_target_orr,
 sample_target_km_est, spop_target_km_est, spop_target_obs_cens_km_est,
 sample_ms_km_est, spop_ms_km_est,
 sample_km_est, spop_km_est,
 sample_target_quant_pfs, spop_target_quant_pfs,
 sample_target_quant_pfs_exceeds_max, spop_target_quant_pfs_exceeds_max,
 sample_ms_quant_pfs, spop_ms_quant_pfs,
 sample_ms_quant_pfs_exceeds_max, spop_ms_quant_pfs_exceeds_max,
 sample_quant_pfs, spop_quant_pfs,
 sample_quant_pfs_exceeds_max, spop_quant_pfs_exceeds_max,
 sample_target_pfs_n, spop_target_pfs_n,
 sample_ms_pfs_n, spop_ms_pfs_n,
 sample_pfs_n, spop_pfs_n,
 sample_os_km_est, spop_os_km_est,
 sample_os_quant, spop_os_quant,
 sample_os_quant_exceeds_max, spop_os_quant_exceeds_max,
 sample_os_n, spop_os_n) =
  aggregate_trial_metrics(
    sample_target_confirmed_response,
    spop_target_confirmed_response,
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
    hmc_trial_patient_pos,
    max_all_t,
    pfs_quantiles,
    pfs_timepoints
  );

// Aggregate to conditional group-level metrics
(cond_sample_target_orr, cond_spop_target_orr,
 cond_sample_target_km_est, cond_spop_target_km_est, cond_spop_target_obs_cens_km_est,
 cond_sample_ms_km_est, cond_spop_ms_km_est,
 cond_sample_km_est, cond_spop_km_est,
 cond_sample_target_quant_pfs, cond_spop_target_quant_pfs,
 cond_sample_target_quant_pfs_exceeds_max, cond_spop_target_quant_pfs_exceeds_max,
 cond_sample_ms_quant_pfs, cond_spop_ms_quant_pfs,
 cond_sample_ms_quant_pfs_exceeds_max, cond_spop_ms_quant_pfs_exceeds_max,
 cond_sample_quant_pfs, cond_spop_quant_pfs,
 cond_sample_quant_pfs_exceeds_max, cond_spop_quant_pfs_exceeds_max,
 cond_sample_target_pfs_n, cond_spop_target_pfs_n,
 cond_sample_ms_pfs_n, cond_spop_ms_pfs_n,
 cond_sample_pfs_n, cond_spop_pfs_n,
 cond_sample_os_km_est, cond_spop_os_km_est,
 cond_sample_os_quant, cond_spop_os_quant,
 cond_sample_os_quant_exceeds_max, cond_spop_os_quant_exceeds_max,
 cond_sample_os_n, cond_spop_os_n) =
  aggregate_conditional_group_metrics(
    sample_target_confirmed_response,
    spop_target_confirmed_response,
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

// Competing Risks CIF (per-trial, empirical subdistribution)
for (s in 1:n_trials) {
  int n_tr = get_pos_size(trial_patient_pos, s);
  if (n_tr > 0) {
    int tr_start; int tr_end;
    (tr_start, tr_end) = get_pos(trial_patient_pos, s);

    (spop_cif_01[s], spop_cif_02[s], spop_cif_03[s]) = compute_trial_cif(
      spop_pfs[tr_start:tr_end], spop_right_censored[tr_start:tr_end],
      spop_is_dropout[tr_start:tr_end],
      spop_os[tr_start:tr_end],  spop_os_censored[tr_start:tr_end],
      max_all_t
    );

    (sample_cif_01[s], sample_cif_02[s], sample_cif_03[s]) = compute_trial_cif(
      sample_pfs[tr_start:tr_end], sample_right_censored[tr_start:tr_end],
      sample_is_dropout[tr_start:tr_end],
      sample_os[tr_start:tr_end],  sample_os_censored[tr_start:tr_end],
      max_all_t
    );
  }
}
```

- [ ] **Step 2: Remove unused `forecast_biomarker_category` param from `calculate_all_patients_endpoints_rng`**

In `stan/modules/state_space/sf.stanfunctions`, remove the `forecast_biomarker_category` parameter from the function signature and the dead `patient_forecast_cat` variable from the body (line 1388). Update all callers (tumor endpoints, LFO endpoints).

---

### Task 4: Update `sf-ssm-log-space.stan` (tumor model)

Replace the old includes with the new structure: categorization → output declarations → local contract block → shared computation.

**Files:**
- Modify: `stan/tumor/sf-ssm-log-space.stan` (generated quantities block, ~lines 180-188)

- [ ] **Step 1: Update the GQ block**

Change from:
```stan
  // Biomarker-agnostic trajectory generation
  #include "modules/state_space/generated_quantities.stan"

  // SLD-specific: RECIST classification + PFS endpoints
  #include "_tumor_endpoints_generated_quantities.stan"

  // RECIST accuracy metrics
  #include "modules/tumor/generated_quantities.stan"
```

To:
```stan
  // Biomarker-agnostic trajectory generation
  #include "modules/state_space/generated_quantities.stan"

  // SLD-specific: aliases + RECIST categorization + assessment-visit RECIST
  #include "_tumor_categorization.stan"

  // Endpoint output declarations (GQ scope — saved to CSV)
  array[n_hmc_patients] int<lower = 0> sample_target_pfs, spop_target_pfs, sample_ms_pfs, spop_ms_pfs, spop_target_obs_cens_pfs,
                                   sample_pfs, spop_pfs;
  array[n_hmc_patients] int<lower = 0, upper = 1>
    sample_target_right_censored, spop_target_right_censored, spop_target_obs_cens_right_censored,
    sample_ms_right_censored, spop_ms_right_censored,
    sample_right_censored, spop_right_censored;
  array[n_hmc_patients] int<lower = 0> sample_os, spop_os;
  array[n_hmc_patients] int<lower = 0, upper = 1> sample_os_censored, spop_os_censored;
  array[n_patients] int<lower = 0, upper = 1> spop_is_dropout, sample_is_dropout;
  array[sum(target_right_censored)] int<lower = 0> forecast_target_pfs;
  array[sum(target_right_censored)] int<lower = 0, upper = 1> forecast_target_right_censored;
  array[n_hmc_patients] int<lower = 0, upper = 1> sample_target_confirmed_response, spop_target_confirmed_response;
  array[n_hmc_patients] int<lower = 0, upper = 1> sample_target_unconfirmed_response, spop_target_unconfirmed_response;
  vector<lower = 0, upper = 1>[n_trials] sample_target_orr, spop_target_orr;
  vector<lower = 0, upper = 1>[n_cond_group] cond_sample_target_orr = zeros_vector(n_cond_group), cond_spop_target_orr = zeros_vector(n_cond_group);
  // KM estimates, quantiles, PFS-at-timepoints, OS — copy from _tumor_endpoints_generated_quantities.stan lines 29-101
  // CIF arrays — copy from _tumor_endpoints_generated_quantities.stan lines 375-381
  // (Full declarations omitted for brevity — copy verbatim from current file)

  // Shared endpoint computation (contract vars are local — not saved to CSV)
  {
    array[n_total_visits] int obs_biomarker_cat = recist;
    array[n_total_visits] int rep_biomarker_cat = rep_recist;
    array[n_total_forecast_obs_visits] int forecast_obs_biomarker_cat = forecast_obs_recist;
    array[n_patients] int burden_pfs = pfs;
    array[n_patients] int burden_right_censored = right_censored;
    array[n_patients] int burden_target_pfs = target_pfs;
    array[n_patients] int burden_target_right_censored = target_right_censored;

    #include "modules/state_space/burden_endpoints.stan"
  }

  // RECIST accuracy metrics
  #include "modules/tumor/generated_quantities.stan"
```

- [ ] **Step 2: Compile tumor model**

```bash
~/.cmdstan/cmdstan-2.38.0/bin/stanc --include-paths=stan --include-paths=stan/tumor stan/tumor/sf-ssm-log-space.stan
```

---

### Task 5: Update `pioneer.stan` (PSA model)

Add visits module, replace PSA endpoint includes with new structure.

**Files:**
- Modify: `stan/psa/pioneer.stan` (data, transformed_data, generated quantities blocks)

- [ ] **Step 1: Add visits module to data and transformed_data blocks**

In the `data` block, add:
```stan
  #include "modules/visits/data.stan"
```

In the `transformed data` block, add:
```stan
  #include "modules/visits/transformed_data.stan"
```

- [ ] **Step 2: Update the GQ block**

Change from:
```stan
  // Biomarker-agnostic trajectory generation
  #include "modules/state_space/generated_quantities.stan"

  // PSA-specific endpoints (PCWG3 classification)
  #include "_psa_endpoints_generated_quantities.stan"

  // Multistate endpoints (PFS, KM curves, OS) using PCWG3 categories
  #include "_psa_ms_endpoints_generated_quantities.stan"
```

To:
```stan
  // Biomarker-agnostic trajectory generation
  #include "modules/state_space/generated_quantities.stan"

  // PSA-specific: aliases + PCWG3 categorization + assessment-visit PCWG3
  #include "_psa_categorization.stan"

  // Endpoint output declarations (GQ scope — saved to CSV)
  // Same declarations as tumor model (copy from Task 4) except:
  //   forecast_target_pfs sized by sum(psa_right_censored) instead of sum(target_right_censored)
  array[n_hmc_patients] int<lower = 0> sample_target_pfs, spop_target_pfs, ...;
  // ... (full declarations — same pattern as tumor) ...
  array[sum(psa_right_censored)] int<lower = 0> forecast_target_pfs;
  array[sum(psa_right_censored)] int<lower = 0, upper = 1> forecast_target_right_censored;
  // ... KM, quantiles, OS, CIF arrays ...

  // Shared endpoint computation (contract vars are local — not saved to CSV)
  {
    array[n_total_visits] int obs_biomarker_cat = pcwg3_category;
    array[n_total_visits] int rep_biomarker_cat = rep_pcwg3;
    array[n_total_forecast_obs_visits] int forecast_obs_biomarker_cat = forecast_obs_pcwg3;
    // PSA has no target vs non-target distinction: target = overall
    array[n_patients] int burden_pfs = psa_pfs;
    array[n_patients] int burden_right_censored = psa_right_censored;
    array[n_patients] int burden_target_pfs = psa_pfs;
    array[n_patients] int burden_target_right_censored = psa_right_censored;

    #include "modules/state_space/burden_endpoints.stan"
  }
```

- [ ] **Step 3: Compile pioneer model**

```bash
~/.cmdstan/cmdstan-2.38.0/bin/stanc --include-paths=stan --include-paths=stan/psa stan/psa/pioneer.stan
```

---

### Task 6: Update R pipeline for visits module in pioneer

The visits module requires `forecast_observation_interval` in the Stan data. This must be provided by the R pipeline for pioneer.

**Files:**
- Modify: `r/priors.R` or the pioneer data preparation function

- [ ] **Step 1: Check if `forecast_observation_interval` is already in pioneer stan data**

Search for where `forecast_observation_interval` is set in R code for tumor, and ensure it's also set for pioneer.

- [ ] **Step 2: Add if missing**

The visits module expects `forecast_observation_interval` (typically 6 for 6-week assessment cycles). Add to the pioneer data preparation if not already present.

---

### Task 7: Update LFO model caller

The LFO model (`sf-ssls-lfo.stan`) uses `_lfo_endpoints_generated_quantities.stan` which calls `calculate_all_patients_endpoints_rng` with a different patient set (cutoff patients). This is NOT part of the unification — LFO stays as-is. But if we remove `forecast_biomarker_category` from the function signature (Task 3 Step 2), we must also update the LFO caller.

**Files:**
- Modify: `stan/tumor/_lfo_endpoints_generated_quantities.stan`

- [ ] **Step 1: Remove `forecast_biomarker_category` arg from LFO caller if removed from function**

- [ ] **Step 2: Compile LFO model**

```bash
~/.cmdstan/cmdstan-2.38.0/bin/stanc --include-paths=stan --include-paths=stan/tumor stan/tumor/sf-ssls-lfo.stan
```

---

### Task 8: Verify ms_os_event_12 availability in pioneer

The endpoint function takes `ms_os_event_12` (exact death calendar week for observed 1→2 deaths). Verify this field is defined in the PSA/pioneer multistate data.

**Files:**
- Check: `stan/modules/multistate/data.stan` for `ms_os_event_12`
- Check: `r/pioneer/` data prep for `ms_os_event_12` population

- [ ] **Step 1: Grep for ms_os_event_12 in multistate module and R code**

If missing from multistate module data.stan, it was added by data-in specifically for tumor. Need to add it to the shared multistate module or handle it in the burden endpoints file.

---

### Task 9: Final compilation + tests

- [ ] **Step 1: Compile all three models**

```bash
~/.cmdstan/cmdstan-2.38.0/bin/stanc --include-paths=stan --include-paths=stan/tumor stan/tumor/sf-ssm-log-space.stan
~/.cmdstan/cmdstan-2.38.0/bin/stanc --include-paths=stan --include-paths=stan/psa stan/psa/pioneer.stan
~/.cmdstan/cmdstan-2.38.0/bin/stanc --include-paths=stan --include-paths=stan/tumor stan/tumor/sf-ssls-lfo.stan
```

- [ ] **Step 2: Run R tests**

```bash
Rscript -e 'testthat::test_dir("tests/testthat")'
```

- [ ] **Step 3: Commit**

```bash
git add stan/modules/state_space/burden_endpoints.stan \
        stan/tumor/_tumor_categorization.stan \
        stan/psa/_psa_categorization.stan \
        stan/tumor/sf-ssm-log-space.stan \
        stan/psa/pioneer.stan \
        stan/modules/state_space/sf.stanfunctions \
        stan/tumor/_lfo_endpoints_generated_quantities.stan
git commit -m "Unify tumor and PSA endpoint logic into shared burden module

Create modules/state_space/burden_endpoints.stan with shared endpoint
computation, aggregation, and CIF logic. Each model provides a
categorization file (_tumor_categorization.stan, _psa_categorization.stan)
that computes biomarker-specific categories at GQ scope. Contract
interface variables are set in a local { } block in each model file
(not saved to CSV), with the shared computation included inside.

- Add visits module to pioneer for assessment-visit scheduling
- Remove unused forecast_biomarker_category from endpoint function
- Eliminate ~250 lines of duplication between tumor and PSA endpoints
- No new GQ output variables — contract vars are local scope"
```

---

## Risks and Mitigations

1. **RNG sequence changes**: The tumor model currently draws trajectory noise twice (state_space GQ + tumor endpoints GQ). Removing the redundant draw changes the RNG sequence, which changes all downstream samples. This is correct but breaks bitwise reproducibility. **Mitigation**: Document this in the commit message.

2. **LFO model divergence**: LFO uses different sizing and patient sets. It stays on its own endpoint file. **Mitigation**: Only update LFO for function signature changes, not for the burden unification.

3. **PSA assessment-visit PCWG3 correctness**: The PCWG3 computation requires absolute PSA, study days, and nadir tracking — more complex than RECIST. **Mitigation**: The code in Task 2 follows the same pattern as the existing weekly PCWG3 computation in `_psa_endpoints_generated_quantities.stan`, adapted for assessment-visit subsampling.

4. **Endpoint output declarations duplicated across models**: Each model declares ~50 lines of endpoint output variables in its GQ block. These are nearly identical (only `forecast_target_pfs` sizing differs). This is intentional — each model explicitly owns its output interface. **Mitigation**: Accept the duplication; it's declarations only, no logic.

---

## Review Notes

Plan reviewed by subagent. Issues triaged and resolved:

| # | Reviewer Issue | Status | Resolution |
|---|---------------|--------|------------|
| 1 | `measure_nu_obs`/`measure_sd_obs` wrong | **False positive** | Both models define these in GQ before include (pioneer.stan:219-220, sf-ssm-log-space.stan:177-178) |
| 3 | Shared endpoint call lacks full arg list | **Fixed** | Added explicit 26-value LHS + full RHS argument list in Task 3 |
| 5 | Missing `spop_is_dropout`/`sample_is_dropout` | **Fixed** | Included in full LHS capture |
| 6 | CIF requires `trial_patient_pos` in PSA | **False positive** | Available via `_base_hierarchy_transformed_data.stan:55`, included by both models |
| 8 | Missing `burden_*` in tumor categorization | **Fixed** | Contract vars now in model file's local block, not categorization file |
| 9 | `forecast_target_pfs` sizing | **Fixed** | Each model uses its own data array for sizing (`target_right_censored` / `psa_right_censored`) |
| — | Contract vars saved to CSV (doubles output) | **Fixed** | Contract vars declared in local `{ }` block — not saved. Output decls at GQ scope in model file |
