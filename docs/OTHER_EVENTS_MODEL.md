# Other Events Model: Complete Guide

**Last Updated:** November 14, 2025

## Table of Contents

1. [Overview](#overview)
2. [Model Evolution: From Competing Risks to Independent Events](#model-evolution)
3. [Data Computation](#data-computation)
4. [Mechanistic Tumor Covariates](#mechanistic-tumor-covariates)
5. [Implementation Details](#implementation-details)

---

## Overview

The **other events model** captures progression-free survival events that are **not** due to target lesion progression. This includes:

- Non-target lesion progression
- New lesions
- Death
- Study dropout

This model works in conjunction with the target lesion (SLD) model to provide a complete picture of disease progression.

---

## Model Evolution

### From Competing Risks to Independent Events

**Date of Change:** November 12, 2025  
**Branch:** karim/non-target

#### Problem Discovery

Found patient cases that violated the competing risks assumption:

**Example:**
- Non-target PD at week 6
- Patient continued study participation
- Target PD at week 18

This revealed that some "other events" (e.g., non-target PD) are **non-terminal** - they don't end study participation. Patients can experience both non-target AND target progression.

#### Previous Model: Competing Risks

```stan
// Computed other events censoring based on target PFS
other_events_right_censored[i] = right_censored[i] || 
  (!target_right_censored[i] && pfs[i] >= target_pfs[i]);

// Assertion: at most one event type
if (!target_right_censored[i] && !other_events_right_censored[i]) {
  reject("Both events uncensored - violates competing risks");
}
```

**Problem:** This assertion failed when both events were observed for the same patient.

#### Current Model: Independent Events

Target progression and other events are now modeled as **independent processes**:

- Both can be observed for the same patient
- Each has its own survival function
- No mutual exclusivity constraint

```stan
// Other events data passed as INPUT
array[n_patients] int<lower=0> other_events_pfs;
array[n_patients] int<lower=0, upper=1> other_events_right_censored;

// Simple interval censoring (no competing risks logic)
other_events_interval_censored[i] = 
  other_events_right_censored[i] ? 0 : interval_censored[i];
```

#### Precedence Rule

**Important:** If target PD and other-event occur **in the same week**, target takes precedence.

**Rationale:**
- Clinical reporting treats target progression as the primary event
- Simplifies interpretation when events coincide
- Reflects treatment decision-making (target progression drives decisions)

**Implementation:**
- Non-target event at same visit as target PD is treated as **not occurring**
- Only target PD is recorded
- `non_target_pd = FALSE` in this scenario

---

## Data Computation

### Key Variables

#### From Trial Data

1. **`response`** - Overall RECIST response (target + non-target)
2. **`det_response`** - Target-only RECIST response (from SLD model)
3. **`progression_before_death`** - PD occurred before death
4. **`death`** - Patient died
5. **`trial_death`** - Death within 6 weeks of last assessment

#### Computed Variables

1. **`non_target_pd`** - Non-target progression occurred
2. **`other_events_pfs`** - Time to other event (weeks)
3. **`other_events_right_censored`** - Whether other events were censored

### Computation Logic

#### Step 1: Detect Non-Target PD

```r
non_target_pd = map_lgl(visit_data, \(d) {
  # Find first overall PD
  overall_pd_idx <- which(fct_match(d$response, "PD"))
  if (length(overall_pd_idx) == 0) return(FALSE)
  
  # Find first target PD
  target_pd_idx <- which(fct_match(d$det_response, "PD"))
  if (length(target_pd_idx) == 0) return(TRUE)
  
  # Non-target PD only if STRICTLY BEFORE target PD
  # PRECEDENCE RULE: Same visit → target takes precedence
  min(overall_pd_idx) < min(target_pd_idx)
})
```

**Logic:**
- Overall PD but no target PD → Non-target PD
- Both PD types → Check which occurred first
- **Same visit index** → Target precedence, `non_target_pd = FALSE`

#### Step 2: Compute Other Events Variables

```r
# Other events censored if NEITHER non-target PD NOR death occurred
other_events_right_censored = !(non_target_pd | trial_death)

# If other event occurred, use pfs; otherwise use last visit
other_events_pfs = if_else(other_events_right_censored, patient_max_t, pfs)
```

### Example Scenarios

#### Scenario 1: Non-Target PD Only

**Data:**
- Non-target PD at week 6
- No target PD
- Last visit week 18

**Results:**
```r
non_target_pd = TRUE
trial_death = FALSE
other_events_right_censored = FALSE
other_events_pfs = 6

target_right_censored = TRUE
target_pfs = 18
```

#### Scenario 2: Both Event Types

**Data:**
- Non-target PD at week 6
- Target PD at week 18

**Results:**
```r
non_target_pd = TRUE
other_events_right_censored = FALSE
other_events_pfs = 6

target_right_censored = FALSE
target_pfs = 18
```

Both events are independently modeled ✓

#### Scenario 3: Simultaneous Events (Precedence Rule)

**Data:**
- Both non-target and target PD at week 12 (same visit)

**Results:**
```r
non_target_pd = FALSE  # Precedence rule applies
other_events_right_censored = TRUE
other_events_pfs = patient_max_t

target_right_censored = FALSE
target_pfs = 12
```

Only target progression recorded ✓

#### Scenario 4: Death Without PD

**Data:**
- Death at week 8
- No PD

**Results:**
```r
non_target_pd = FALSE
trial_death = TRUE
other_events_right_censored = FALSE
other_events_pfs = 8

target_right_censored = TRUE
target_pfs = patient_max_t
```

#### Scenario 5: Late Death (Censored)

**Data:**
- Death at week 50
- Last assessment week 40

**Results:**
```r
trial_death = FALSE  # Death > 6 weeks after last visit
other_events_right_censored = TRUE
other_events_pfs = 40
```

---

## Mechanistic Tumor Covariates

**Date Implemented:** November 13, 2025  
**Updated:** November 14, 2025 (Normalization approach)

### Motivation

The mechanistic state-space model provides rich patient-specific information beyond current tumor burden:

1. **Decrease rate** - How fast tumor regresses (treatment effect)
2. **Growth rate** - How fast tumor grows (disease aggressiveness)  
3. **SLD velocity** - Rate of tumor size change (derivative) **[DISABLED]**

These capture underlying biology predictive of other events independent of current tumor size.

### Three Tumor Covariates (Current)

The model estimates **3 population-level coefficients per cause**:

| Index | Covariate | Type | Normalization | Interpretation | Expected Sign |
|-------|-----------|------|---------------|----------------|---------------|
| 1 | log(SLD) | Time-varying | Median/IQR from observed data | Current tumor burden | Positive (+) |
| 2 | log(decrease rate) | Time-invariant | None (patient-normalized) | Treatment response rate | Negative (-) |
| 3 | log(growth rate) | Time-invariant | None (patient-normalized) | Disease aggressiveness | Positive (+) |
| ~~4~~ | ~~SLD velocity~~ | ~~Time-varying~~ | ~~DISABLED~~ | ~~Rate of tumor change~~ | ~~Positive (+)~~ |

**Note:** SLD velocity (covariate 4) is disabled (`n_tumor_covar = 3`) as it may be redundant with rates.

### Implementation

Variables are accessed directly from transformed parameters:

```stan
// In other_events/transformed_parameters.stan

// Patient-level mechanistic rates (already patient-normalized)
patient_log_decrease_rate[i]  // Time-invariant, no normalization needed
patient_log_growth_rate[i]    // Time-invariant, no normalization needed

// Time-varying SLD from states_full_grid
// states_full_grid[1] = log(regression component)
// states_full_grid[2] = log(growth component)
// log(SLD) = log_sum_exp(log_regression, log_growth)
```

### Normalization Strategy

#### 1. log(SLD) - Median/IQR Normalization

**Purpose:** Capture absolute tumor burden effect

**Method:** 
- Compute median and IQR from **observed SLD data** (in `transformed_data`)
- Filters out zeros to avoid `-Inf`
- Fixed across MCMC iterations for stable prior interpretation

```stan
// In other_events/transformed_data.stan
vector[sum(n_patient_visits)] log_sld_all_obs;
int n_positive = 0;

for (i in 1:sum(n_patient_visits)) {
  if (sum_tumor_size[i] > 0) {
    n_positive += 1;
    log_sld_all_obs[n_positive] = log(sum_tumor_size[i]);
  }
}

array[3] real quantiles_obs = quantile(log_sld_all_obs[1:n_positive], {0.25, 0.5, 0.75});
median_log_sld_obs = quantiles_obs[2];
iqr_log_sld_obs = quantiles_obs[3] - quantiles_obs[1];
```

**Application:**
```stan
// In transformed_parameters
log_sld_standardized = (log_sld_absolute - median_log_sld_obs) / iqr_log_sld_obs;
```

**Interpretation:** "log HR per IQR increase in log(SLD)"

**Why median/IQR instead of mean/SD?**
- Robust to outliers and extreme latent state values
- Avoids `-Inf` from zero observed SLD
- More appropriate for skewed distributions

#### 2. Patient Rates - No Normalization

**Purpose:** Capture patient-specific rate of change

**Method:** 
- Rates are already **patient-normalized** (relative to each patient's baseline)
- Used directly without population-level normalization
- Captures how fast tumor changes relative to patient's own baseline

```stan
// In transformed_parameters - used raw
oe_time_varying_log_hazard_ratio[k, i] += oe_tumor_coef_pop[k][2] * patient_log_decrease_rate[i];
oe_time_varying_log_hazard_ratio[k, i] += oe_tumor_coef_pop[k][3] * patient_log_growth_rate[i];
```

**Interpretation:**
- "log HR per unit increase in patient's log(decrease rate)"
- "log HR per unit increase in patient's log(growth rate)"

**Why no normalization?**
- Rates measure **relative change**, not absolute burden
- Already normalized per patient (inherent in state-space model)
- Population normalization would obscure patient-specific dynamics

#### 3. SLD Velocity - DISABLED

Previously computed as first derivative of log(SLD):
```stan
sld_velocity[t] = log_sld[t] - log_sld[t-1]
```

**Why disabled?**
- Potentially redundant with patient-level rates
- Rates already capture tumor dynamics characteristics
- Velocity may just add noise without new information

Can be re-enabled by setting `n_tumor_covar = 4` in R.

### Prior Specifications

Weakly informative priors (per cause):

```r
# n_tumor_covar = 3
oe_tumor_coef_pop_mean = array(rep(0, 3), dim = c(n_causes, 3))
oe_tumor_coef_pop_sd   = array(rep(0.5, 3), dim = c(n_causes, 3))
```

Expected effects:
- **log(SLD):** Positive - Higher burden (relative to median) → Higher risk
- **log(decrease rate):** Negative - Faster regression → Lower risk
- **log(growth rate):** Positive - Faster growth → Higher risk

### Configuration

Set in R data preparation (`r/sclc/prepare_analysis_data.R`):

```r
# Use SLD + rates (current default)
n_tumor_covar <- 3

# Or use only log(SLD)
n_tumor_covar <- 1

# Or enable velocity (all covariates)
n_tumor_covar <- 4
```

---

## Implementation Details

### Files Modified

#### Stan Files

1. **`stan/ssls/modules/other_events/data.stan`**
   - Added `other_events_pfs` and `other_events_right_censored` as input data
   - Removed `n_tumor_covar` from data block (moved to transformed_data)

2. **`stan/ssls/modules/other_events/transformed_data.stan`**
   - Set `n_tumor_covar = 4` (hardcoded)
   - Removed competing risks logic
   - Removed assertion about mutual exclusivity
   - Kept interval censoring logic

3. **`stan/ssls/modules/other_events/transformed_parameters.stan`**
   - Access patient-level rates directly
   - Compute SLD velocity
   - Z-score normalize all covariates
   - Broadcast time-invariant rates to all timepoints

4. **`stan/ssls/modules/other_events/parameters.stan`**
   - `oe_tumor_coef_pop` sized based on `n_tumor_covar`

#### R Files

5. **`r/sclc/prepare_analysis_data.R`**
   - Compute `non_target_pd` variable
   - Compute `other_events_pfs` and `other_events_right_censored`
   - Pass these to Stan data
   - Set `n_tumor_covar = 4`

6. **`r/priors.R`**
   - `oe_tumor_coef_pop_mean/sd` use `stan_data$n_tumor_covar`
   - Provides priors for all 4 coefficients

### Data Flow

```
Raw trial data
    ↓
r/sclc/prepare_analysis_data.R
    ↓ (compute non_target_pd)
    ↓ (compute other_events_pfs/right_censored)
    ↓
Stan data list
    ↓
stan/ssls/modules/other_events/
    ↓ (access patient rates)
    ↓ (compute SLD velocity)
    ↓ (z-score normalize)
    ↓
Hazard model with 4 tumor covariates
```

### Backward Compatibility

The implementation maintains compatibility:

1. If `n_tumor_covar = 1`: Only log(SLD) coefficient (previous behavior)
2. If `n_tumor_covar = 4`: All covariates including rates and velocity
3. Feature can be toggled via R configuration

---

## Testing Recommendations

### Correctness Validation

1. **Verify non-target PD detection** - Check patients with known non-target events
2. **Test precedence rule** - Ensure simultaneous events handled correctly
3. **Compare against target PFS** - Both events should be independent

### Model Comparison

1. **Fit with only log(SLD)** (`n_tumor_covar = 1`) for baseline
2. **Fit with all covariates** (`n_tumor_covar = 4`) for full model
3. **Compare via LOO-CV** to quantify predictive improvement
4. **Check coefficient signs** match clinical expectations

### Clinical Plausibility

1. **Examine patient-level predictions** - Do hazards make sense?
2. **Validate against observed events** - Are high-risk patients progressing?
3. **Check rate effects** - Do they align with treatment response?

---

## Future Extensions

Potential additional tumor covariates:

- **SLD acceleration**: Second derivative `d²(log_SLD)/dt²`
- **Growth fraction**: Proportion in growing state
- **Resistance score**: Increase in growth rate after nadir
- **Time since nadir**: Duration at minimum SLD
- **Nadir depth**: How low SLD went

These could be added to `states_full_grid` and accessed similarly.

---

## References

### Related Documentation

- **Architecture:** [`docs/ARCHITECTURE.md`](ARCHITECTURE.md)
- **State space model:** `sld_state_space_model.md`
- **Code ownership:** [`docs/CODEOWNERS`](CODEOWNERS)

### Git History

- Independent events change: Commits a553b2b, 988ec3b
- Mechanistic covariates: November 13, 2025

---

**Document Version:** 2.0 (Consolidated)  
**Status:** Active Reference  
**Maintainer:** Update when other events model changes
