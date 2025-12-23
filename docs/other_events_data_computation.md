# Computing Other Events Data for Independent Events Model

**Date**: November 12, 2025  
**Purpose**: Document how `other_events_pfs` and `other_events_right_censored` are computed

## Overview

For the independent events model, we need to separate:
- **Target progression**: From tumor SLD dynamics (already computed as `det_pfs`, `det_right_censored`)
- **Other events**: Non-target PD, death, dropout (independent from target)

**Precedence Rule**: If target and other-event occur in the **same week**, target takes precedence. This means:
- Non-target PD at same visit as target PD is ignored
- Only target PD is recorded for that week
- Clinically, target progression is the primary driver for treatment decisions

## Key Data Sources

### From Trial Data

1. **`response`**: Overall RECIST response (includes target + non-target lesions)
2. **`det_response`**: Target-only RECIST response (from SLD trajectory)
3. **`progression_before_death`**: Flag indicating PD occurred before death
4. **`death`**: Whether patient died
5. **`trial_death`**: Death within 6 weeks of last assessment

## Computation Logic

### Step 1: Detect Non-Target PD

Non-target PD occurs when overall PD is detected but it's NOT due to target lesions:

```r
non_target_pd = map_lgl(visit_data, \(d) {
  # Find first overall PD
  overall_pd_idx <- which(fct_match(d$response, "PD"))
  if (length(overall_pd_idx) == 0) return(FALSE)  # No PD at all
  
  # Find first target PD
  target_pd_idx <- which(fct_match(d$det_response, "PD"))
  if (length(target_pd_idx) == 0) return(TRUE)   # Overall PD but no target PD
  
  # Non-target PD only if it occurred STRICTLY BEFORE target PD
  # PRECEDENCE RULE: If same visit, target takes precedence
  min(overall_pd_idx) < min(target_pd_idx)
})
```

**Logic**:
- If overall `response == "PD"` but `det_response != "PD"`: Non-target PD
- If both are PD, check which occurred first (by visit index)
- **If they occur at same visit index**: Target takes precedence, `non_target_pd = FALSE`
- Non-target PD includes: non-target lesion progression, new lesions

### Step 2: Flag Progression Events

```r
progression_event = !right_censored & progression_before_death
```

This indicates whether an uncensored event was due to PD (vs death without PD).

### Step 3: Compute Other Events Variables

```r
# Other events are censored if NEITHER non-target PD NOR death occurred
other_events_right_censored = !(non_target_pd | trial_death)

# If other event occurred, use pfs; otherwise use last visit
other_events_pfs = if_else(other_events_right_censored, patient_max_t, pfs)
```

**Important**: `pfs` here is the overall PFS which may be:
- Time to non-target PD (if that occurred)
- Time to death (if that occurred without non-target PD)
- Last visit (if censored)

## Example Scenarios

### Scenario 1: Non-Target PD Only

**Data**:
- Non-target PD at week 6
- No target PD observed
- Patient continued to week 18 (last visit)

**Results**:
- `non_target_pd = TRUE`
- `trial_death = FALSE`
- `other_events_right_censored = FALSE` (non-target event occurred)
- `other_events_pfs = 6`
- `target_right_censored = TRUE` (no target PD)
- `target_pfs = 18`

### Scenario 2: Both Non-Target and Target PD

**Data**:
- Non-target PD at week 6
- Target PD at week 18

**Results**:
- `non_target_pd = TRUE` (occurred before target)
- `trial_death = FALSE`
- `other_events_right_censored = FALSE`
- `other_events_pfs = 6` (first event)
- `target_right_censored = FALSE`
- `target_pfs = 18`

### Scenario 2b: Simultaneous Target and Non-Target PD

**Data**:
- Both non-target and target PD detected at week 12 (same visit)

**Results**:
- `non_target_pd = FALSE` (precedence rule: target takes priority)
- `trial_death = FALSE`
- `other_events_right_censored = TRUE` (no other-event recorded)
- `other_events_pfs = patient_max_t`
- `target_right_censored = FALSE`
- `target_pfs = 12`

### Scenario 3: Target PD Only

**Data**:
- Target PD at week 12
- No non-target PD
- No death

**Results**:
- `non_target_pd = FALSE`
- `trial_death = FALSE`
- `other_events_right_censored = TRUE` (no other-event)
- `other_events_pfs = patient_max_t`
- `target_right_censored = FALSE`
- `target_pfs = 12`

### Scenario 4: Death Without PD

**Data**:
- Death at week 8
- No PD of any kind

**Results**:
- `non_target_pd = FALSE`
- `trial_death = TRUE`
- `other_events_right_censored = FALSE` (death is an other-event)
- `other_events_pfs = 8`
- `target_right_censored = TRUE`
- `target_pfs = patient_max_t`

### Scenario 5: Late Death (Censored)

**Data**:
- Death at week 50
- Last assessment at week 40
- No PD observed

**Results**:
- `trial_death = FALSE` (death > 6 weeks after last visit)
- `non_target_pd = FALSE`
- `other_events_right_censored = TRUE` (death too late, treated as censored)
- `other_events_pfs = 40`
- `target_right_censored = TRUE`
- `target_pfs = 40`

## Important Considerations

### Timing of Overall PFS

The `pfs` variable from raw data represents the time to first event of ANY type. When `other_events_right_censored = FALSE`, this `pfs` should represent the time to the other-event (non-target PD or death).

### Interval Censoring

Both target and other-events use the same `interval_censored` value:
- Shared because assessments occur on the same schedule
- If both events observed, we assume same assessment interval uncertainty applies to whichever occurred

### Edge Cases

1. **Simultaneous events (PRECEDENCE RULE)**: If non-target PD and target PD occur at same visit:
   - Target takes precedence
   - `non_target_pd = FALSE`
   - Only target PD is recorded
   - `other_events_right_censored = TRUE` (unless death also occurred)

2. **Death same week as PD**: 
   - `progression_before_death` flag determines precedence
   - If `progression_before_death = FALSE`, death is the event
   - PFS time is same for both

3. **Non-target PD followed by target PD**:
   - Both can be observed (independent events)
   - `other_events_pfs` = time to non-target event
   - `target_pfs` = time to target event
   - Different times allowed

## Files Modified

- `r/sclc/prepare_analysis_data.R`: Added `non_target_pd`, `progression_event`, `other_events_pfs`, `other_events_right_censored` computation

## References

- Independent events model change: `docs/independent_events_model_change.md`
- PFS variables reference: `pfs_censoring_variables_reference.md`
