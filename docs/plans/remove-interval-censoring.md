# Plan: Remove Interval Censoring — Make `pfs` Mean Detection Week

## Context

The model currently splits PFS into two variables:
- `pfs` = last clean assessment week (1-indexed)
- `interval_censored` = gap to detection: `progress_week - pfs - 1`

This forces `pfs + interval_censored + 1` scattered throughout Stan and R to reconstruct the actual detection week. It caused bugs: the multistate likelihood placed events at the wrong week, the GP index for post-progression death was inconsistent between likelihood and GQ, and sojourn times (`ms_time_12`) were inflated by the IC gap.

**Goal**: `pfs` means "the week progression was detected" everywhere. Convert in R data prep. Remove `interval_censored` from Stan entirely.

**GitHub Issue**: #126 (continuation — Layer 2 assessment gating is done, this is Layer 3)

## Partially-Applied Changes to Revert/Complete

These changes were started in the current session and must be handled:

| File | Current State | Action |
|------|--------------|--------|
| `stan/modules/multistate/transformed_data.stan` | Has 50 new lines computing `detection_week`, `ms_detection_time_01`, `ms_detection_time_12` | **Remove** — conversion moves to R |
| `stan/multistate.stanfunctions` | `calc_ms_single_transition_loglik` updated with `-1` for event patients | **Keep** — correct for detection-week convention |
| `stan/sf-ssm-log-space.stan` | Model block passes `ms_detection_time_01` | **Revert** to `ms_time_01` (now detection-week from R) |
| `stan/_lfo_transformed_data.stan` | References `ms_detection_time_01` | **Revert** to `ms_time_01` |
| `stan/modules/state_space/sf.stanfunctions` | Signature has `detection_week`, `interval_censored`, `ms_detection_time_12` | **Simplify** — use `pfs` directly |

## Implementation

### Step 1: R Data Prep — Convert to Detection-Week Convention

**File: `r/sclc/prepare_analysis_data.R`**

Add a second `mutate()` block after the existing one (after line ~130) that converts all "last clean week" variables to detection-week:

```r
|> mutate(
  # Convert to detection-week convention:
  # For event patients: detection_week = old_pfs + IC + 1
  # For censored patients: pfs already = patient_max_t (unchanged)
  pfs = if_else(!right_censored, pfs + interval_censored + 1L, pfs),
  ms_time_01_full = if_else(!ms_censored_01_full, ms_time_01_full + interval_censored + 1L, ms_time_01_full),
  ms_time_01_sld = if_else(!ms_censored_01_sld, ms_time_01_sld + interval_censored + 1L, ms_time_01_sld),
  ms_time_12 = case_when(
    ms_time_12 == 0L ~ 0L,
    !ms_censored_01_full ~ pmax(1L, ms_time_12 - interval_censored - 1L),
    TRUE ~ ms_time_12
  ),
)
```

In stan_data assembly (~line 296):
- **Remove** `interval_censored = analysis_data$interval_censored`
- **Change** `target_pfs` to use detection-week:
  ```r
  target_pfs = if_else(
    analysis_data$det_right_censored,
    analysis_data$det_pfs,
    analysis_data$det_pfs + analysis_data$det_interval_censored + 1L
  ),
  ```

**File: `r/sclc/prepare_ms_standalone_data.R`** (~line 90)
- Remove `interval_censored = analysis_data$interval_censored`

### Step 2: Stan — Remove `interval_censored` Declaration

**File: `stan/modules/tumor/data.stan`** (line 28)
- Remove: `array[n_patients] int<lower = 0> interval_censored;`

**File: `stan/ms-standalone.stan`** (lines 43-46)
- Remove: `interval_censored` declaration and comment

### Step 3: Stan — Remove Transformed Data Conversion Block

**File: `stan/modules/multistate/transformed_data.stan`**
- Remove the entire block added this session (lines 14-50): `ms_detection_time_01`, `ms_detection_time_12`, `detection_week` declarations and loops
- Keep original content starting from "Validate Final State Upper Bound"

### Step 4: Stan — Revert Model Block Variable Names

**File: `stan/sf-ssm-log-space.stan`** (line 84)
- Change `ms_detection_time_01, ms_time_02, ms_detection_time_12` back to `ms_time_01, ms_time_02, ms_time_12`
- (These now contain detection-week values from R)

### Step 5: Stan — Update GQ Function Signature and Body

**File: `stan/modules/state_space/sf.stanfunctions`**

Function signature changes (remove 2 params, rename 1):
- Remove `array[] int detection_week` parameter
- Remove `array[] int interval_censored` parameter
- Rename `array[] int ms_detection_time_12` → `array[] int ms_time_12`

Update doc comment to reflect detection-week convention.

Body changes:
- Line 1476: `sample_target_pfs[i] = detection_week[i]` → `sample_target_pfs[i] = target_pfs[i]`
  (target_pfs is now detection-week from R)
- Line 1505: `assessment_gated_survival_time_rng(..., pfs[i], 1)` — **keep as-is** (pfs is now detection-week, conditioning starts from detection week)
- Line 1508: `sample_ms_pfs[i] = pfs[i] + interval_censored[i] + 1` → `sample_ms_pfs[i] = pfs[i]`
- Line 1630: `survival_time_rng(log_cond_surv_02[i], pfs[i], 1, 0)` — **keep** (pfs is detection-week, correct conditioning point)
- Line 1670: `sample_os[i] = pfs[i] + ms_time_12[i]` — **keep** (detection_week + corrected_sojourn = death_week ✓)
- Line 1679: `pfs[i], ms_time_12[i]` in `sample_post_progression_death_rng` — **keep** (pfs_week=detection_week, sojourn=corrected)
- All `ms_detection_time_12` → `ms_time_12` throughout function body

### Step 6: Stan — Update GQ Callers

**File: `stan/_endpoints_generated_quantities.stan`** (~line 252-253)
- Remove `interval_censored,` from args to `calculate_all_patients_endpoints_rng`

**File: `stan/_ms_standalone_generated_quantities.stan`** (line 97)
- Change `sample_ms_pfs[i] = ms_time_01[i] + interval_censored[i] + 1` → `sample_ms_pfs[i] = ms_time_01[i]`
- Change `spop_ms_pfs[i] += 1` (line 87) → remove the `+= 1` (spop sampling via `survival_time_rng` returns 0-indexed, needs +1... actually keep this, it's not affected by IC removal. The standalone model still uses `survival_time_rng`, not the assessment-gated version)

### Step 7: Stan — Fix LFO Model

**File: `stan/_lfo_transformed_data.stan`**
- Remove `cutoff_interval_censored` declaration (line 133) and all assignments (lines 143, 148)
- Revert `ms_detection_time_01[i]` back to `ms_time_01[i]` (lines 160, 165)

**File: `stan/_lfo_endpoints_generated_quantities.stan`** (line 271)
- Remove `cutoff_interval_censored,` from args to `calculate_all_patients_endpoints_rng`

### Step 8: R Plotting/KM — Simplify

**File: `r/util.R`**
- `km_to_tibble()` (line 556): Change `pfs + interval_censored + 1 - censored` → `pfs + 1 - censored` (lb and ub converge)
- `apply_calendar_cutoff()` (lines 713-733): Remove all `interval_censored` recomputation logic

**File: `r/plot_functions.R`** (line 866)
- Change `aes(pfs + interval_censored + 1)` → `aes(pfs)`

**File: `r/sclc/plot_functions.R`** (line 927)
- Remove `interval_censored = if (endpoint == "os") 0L else interval_censored`

**File: `r/sclc/prepare_sclc_forecast.R`** (line 23)
- Remove `mutate(pfs = pfs + interval_censored + 1)`

**File: `targets/sclc_targets.R`** (lines 773, 783, 794, 805, 816)
- Remove all `interval_censored = 0L` lines from OS KM targets

### Step 9: Other Projects (Lower Priority)

**Files: `r/breast/prepare_analysis_data.R`, `r/endometrial-to-lung/prepare_analysis_data.R`**
- Same pattern as Step 1: add detection-week conversion, remove IC from stan_data

**Files: `targets/endometrial_to_lung_targets*.R`**
- Change `max(pfs + interval_censored)` → `max(pfs)`

### Step 10: Documentation

**File: `quarto/website/documentation/index-conventions.qmd`**
- Update `pfs` definition: "week progression was detected" (not "last clean assessment")
- Remove `detection week = pfs + IC + 1` formula
- Add note about `assessment_gated_survival_time_rng` return convention
- Update conversion table: remove IC columns
- Note that Pitfall 2 (GP index mismatch) is now structurally impossible

## Files NOT Changed

- `r/data_preparation_pipeline/` — upstream pipeline still computes IC (user instruction)
- `r/util.R` `determine_pfs()` — still returns `det_interval_censored` (used in Step 1 conversion)
- `stan/pfs.stanfunctions` `calc_pch_loglik` — keeps IC support internally for legacy callers
- `stan/legacy/` — all legacy files untouched
- `stan/pfs.stanfunctions` `survival_time_rng` — keeps IC parameter (callers will pass 0)

## Key Correctness Arguments

1. **Likelihood consistency**: `multistate_lpmf` treats `time_01` as the event week (hazard at `time_01`). `calc_ms_single_transition_loglik` now subtracts 1 for events before passing to `calc_pch_loglik` (which adds 1 back). Both place the event at detection_week. ✓

2. **GP index alignment**: Likelihood queries post-progression GP starting at `time_01 + 1 = detection_week + 1`. GQ's `sample_post_progression_death_rng` receives `pfs_week = detection_week` and queries `log_cond_surv_12_t[detection_week + s]`. Identical indexing. ✓

3. **OS computation**: `pfs[i] + ms_time_12[i]` = `detection_week + (death_week - detection_week)` = `death_week`. ✓

4. **Sojourn time**: Old: `death_week - last_clean` included IC gap. New: `death_week - detection_week` is correct sojourn. ✓

## Verification

1. **Stan syntax check**: `stanc --include-paths=stan stan/sf-ssm-log-space.stan && stanc --include-paths=stan stan/sf-ssls-lfo.stan && stanc --include-paths=stan stan/ms-standalone.stan`

2. **R round-trip test**: Before/after conversion, verify `pfs_new == pfs_old + interval_censored + 1` for all event patients, and `pfs_new == pfs_old` for all censored patients

3. **Run 200-patient test**: `sclc_targets.sh -b ic-removal-test -t 200`, compare KM curves against previous run with old convention

4. **Edge cases to verify**:
   - Patients with IC=0 (adjacent visits): pfs changes by +1
   - Death without progression: pfs unchanged (death patients have IC=0 by construction)
   - Censored patients: pfs unchanged (= patient_max_t)
   - ms_time_12 with IC > 0: verify sojourn decreases by IC+1, clamped to min 1
