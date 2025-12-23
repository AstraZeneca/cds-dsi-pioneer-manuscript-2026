# Moving from Competing Risks to Independent Events Model

**Date**: November 12, 2025  
**Issue**: Patient data violates competing risks assumption - both target and other-events can be observed  
**Solution**: Model target progression and other-events as independent processes with precedence rule

## Background

### Problem Discovery

Found patient with:
- Non-target progression at week 6
- Continued study participation  
- Target progression at week 18

This violates the competing risks framework which assumes only ONE event type can occur (they compete, and first event ends follow-up).

**Root cause**: Some "other events" (e.g., non-target PD) are **non-terminal** - they don't end study participation. Patient can have both non-target AND target progression.

### Model Framework: Independent Events with Precedence

**Key Assumptions**:
1. **Independence**: Target progression and other-events are modeled as independent processes
2. **Both Can Occur**: A patient can have both non-target PD and target PD
3. **Precedence Rule**: If target PD and other-event occur **in the same week**, target takes precedence
   - This is for clinical reporting/decision purposes
   - Non-target event is treated as NOT occurring (target is the primary event)

### Previous Model (Competing Risks)

```stan
// Other events censored if target progressed first
other_events_right_censored[i] = right_censored[i] || (!target_right_censored[i] && pfs[i] >= target_pfs[i]);

// Assertion: at most one event type can be observed
if (!target_right_censored[i] && !other_events_right_censored[i]) {
  reject("Both events uncensored - violates competing risks");
}
```

**Problem**: This assertion fails when both events are observed.

### New Model (Independent Events)

Target progression (from tumor dynamics) and other-events (death, non-target PD, new lesions, dropout) are **independent processes**:
- Both can be observed for the same patient
- Each modeled separately with its own survival function
- No mutual exclusivity constraint

```stan
// Other events PFS and censoring are INPUT DATA
array[n_patients] int<lower=0> other_events_pfs;
array[n_patients] int<lower=0, upper=1> other_events_right_censored;

// No competing risks logic - just interval censoring
other_events_interval_censored[i] = other_events_right_censored[i] ? 0 : interval_censored[i];
```

## Changes Made

### 1. Stan Data (`stan/ssls/modules/other_events/data.stan`)

**Added** input data for other-events:
```stan
// Other events progression-free survival (weeks from baseline)
array[n_patients] int<lower=0> other_events_pfs;

// Whether other events were right-censored
array[n_patients] int<lower=0, upper=1> other_events_right_censored;
```

### 2. Stan Transformed Data (`stan/ssls/modules/other_events/transformed_data.stan`)

**Removed**:
- Competing risks logic computing `other_events_pfs` and `other_events_right_censored`
- Assertion that fails when both events observed

**Kept**:
- Interval censoring logic (unchanged)
- QR decomposition for covariates (unchanged)

### 3. R Data Preparation (`r/sclc/prepare_analysis_data.R`)

**Added** other-events variables to Stan data:
```r
other_events_pfs = analysis_data$pfs,
other_events_right_censored = analysis_data$right_censored,
```

**IMPORTANT TODO**: Current implementation uses overall `pfs` which may include target PD. Need to:
1. Track non-target PD separately from target PD in data wrangling
2. Compute `other_events_pfs` as time to FIRST of (non-target PD, death, dropout, new lesions)
3. This should be INDEPENDENT from target PD

## Interpretation Changes

### Before (Competing Risks)

- `pfs` = minimum of (target_pfs, other_events_pfs)
- Only one event type can occur
- Other events "compete out" target progression

### After (Independent Events)

- `target_pfs` = time to target lesion progression (from SLD model)
- `other_events_pfs` = time to non-target events (independent)
- Both can be observed for same patient
- No minimum constraint

### Example Patient

**Data**:
- Non-target PD at week 6
- Target PD at week 18
- Patient continued between week 6 and 18

**Variables**:
- `target_pfs = 18`, `target_right_censored = 0` (target PD observed)
- `other_events_pfs = 6`, `other_events_right_censored = 0` (non-target PD observed at week 6)
- Both events modeled independently ✓

**Example with Simultaneous Events**:
- Both non-target and target PD detected at week 12

**Variables**:
- `target_pfs = 12`, `target_right_censored = 0` (target PD at week 12)
- `non_target_pd = FALSE` (precedence rule: target takes priority)
- `other_events_right_censored = TRUE` (treated as no other-event)
- Only target progression is recorded ✓

## Data Requirements

### Current Limitation

The R code currently sets:
```r
other_events_pfs = analysis_data$pfs
```

where `pfs` is the overall PFS from the raw data, which may be the MINIMUM of target and other events.

### Required Fix

Need to modify data wrangling to track:

1. **Target-only events**: Time to target lesion PD (already have this as `det_pfs`)

2. **Other events separately**:
   - Non-target progression
   - New lesions
   - Death
   - Study dropout
   
3. **Compute `other_events_pfs`** as minimum of these, **independent** of target PD

### Where to Fix

Likely in `/mnt/code/r/data_preparation_pipeline/sclc/wrangle_data_sclc.R`:
- Currently computes overall `pfs` from `fjd$lbpfs`
- Need to separate target-specific vs other events
- May need to go back to raw RECIST data to distinguish target vs non-target PD

## Testing

### Before Deployment

1. **Verify data separation**: Check that `other_events_pfs` excludes target-only PD
2. **Test compilation**: Ensure Stan model compiles
3. **Validate results**: Compare with previous competing risks model
4. **Check edge cases**: Patients with only target PD, only other events, or both

### Expected Differences

- Patients with both events will now be modeled (previously would have failed assertion)
- Other-events hazard estimates may change (no longer censored by target PD)
- Predictions should be more realistic for sequential events

## References

- Original bug report: `/mnt/code/fix_summary_sample_target_pfs.md`
- Variable reference: `/mnt/code/pfs_censoring_variables_reference.md`
- Git commits: a553b2b, 988ec3b (introduced PFS variable bugs)

---

**Status**: Code changes complete, but data wrangling needs update to properly separate other-events from target PD.
