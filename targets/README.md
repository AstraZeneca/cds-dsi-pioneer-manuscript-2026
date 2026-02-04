# Targets Pipeline Structure Guide - SCLC-01 Analysis

This README explains how the nested `tar_map()` structures in the SCLC-01 targets files generate target names.

## Overview

The targets pipeline uses nested `tar_map()` calls to create combinations of parameters. Each nesting level adds a suffix to the target name, separated by underscores.

## Naming Convention

Target names follow this pattern:
```
<base_name>_<level1>_<level2>_<level3>
```

Where each level corresponds to a `tar_map()` nesting level.

## Nesting Structure Overview

The pipeline uses **4 nesting levels** (not 5!):

1. **Level 1**: DCO (Data Cut-Off) names - outermost
2. **Level 2**: Covariate groups  
3. **Level 3**: Prior vs Posterior
4. **Level 4**: Two parallel branches (siblings, not nested in each other):
   - **Level 4a**: Patient/trial level (for binned parameter intervals)
   - **Level 4b**: Event type (for patient states and SLD values)

These two Level 4 branches create **different sets of targets** - they don't both apply to the same target.

## Main Target Structure

### Level 1: DCO (Data Cut-Off) Names

The outermost `tar_map()` creates targets for different data cut-offs:

```r
tar_map(
  tibble(
    dco_name = c("apr25", "aug25"),
    ...
  ),
  names = "dco_name",
  ...
)
```

**Suffixes:** `_apr`, `_aug`

**Example targets:**
- `sclc_analysis_data_apr`
- `sclc_analysis_data_aug`

---

### Level 2: Covariate Groups

The second `tar_map()` creates targets for different covariate configurations:

```r
tar_map(
  tribble(
    ~model,  ~model_formula,  ~limited, ~hist, 
    "ctdna",       covar_formula,         TRUE,     TRUE,
    "no_ctdna",    covar_formula_no_ctdna, TRUE,    TRUE,
    "ctdna_only",  covar_formula_ctdna_only, TRUE,  TRUE,
    "all_trials",  covar_formula_no_ctdna, FALSE,   TRUE,
    "no_covar",    NULL,                   TRUE,     TRUE,
    "no_hist",     covar_formula,          TRUE,     FALSE,
  ),
  names = "model",
  ...
)
```

**Suffixes:** `_ctdna`, `_no_ctdna`, `_ctdna_only`, `_all_trials`, `_no_covar`, `_no_hist`

**Example targets:**
- `all_analysis_data_ctdna_apr`
- `all_analysis_data_no_ctdna_aug`
- `tumor_priors_no_hist_apr`

---

### Level 3: Prior vs Posterior

The third `tar_map()` creates targets for both prior and posterior analysis:

```r
tar_map(
  tibble(
    type = c("prior", "posterior"),
    fit_data = c(FALSE, TRUE),
    base_name = c("prior_tumor_ssls", "tumor_ssls"),
    ...
  ),
  names = "type",
  ...
)
```

**Suffixes:** `_prior`, `_posterior`

**Note:** The `base_name` column actually changes the base name itself, not just adds a suffix:
- When `type = "prior"`: targets start with `prior_tumor_ssls_`
- When `type = "posterior"`: targets start with `tumor_ssls_`

**Example targets:**
- `prior_tumor_ssls_stan_data_ctdna_apr`
- `tumor_ssls_stan_data_ctdna_apr`
- `prior_tumor_ssls_res_no_covar_aug`
- `tumor_ssls_res_no_covar_aug`

---

### Level 4a: Patient/Trial Level (Parallel to Level 4b)

A fourth `tar_map()` creates targets for different analysis levels (this is a **sibling** to Level 4b below, not nested):

```r
tar_map(
  tibble(level = c("patient")),
  names = "level",
  ...
)
```

**Suffixes:** `_patient`

**Example targets:**
- `tumor_ssls_rates_bpi_patient_prior_ctdna_apr`
- `tumor_ssls_rates_bpi_patient_posterior_no_hist_aug`
- `tumor_ssls_logits_bpi_patient_prior_no_ctdna_aug`

**Note:** This branch only generates targets for binned parameter intervals (`*_bpi` targets).

---

### Level 4b: Event Type (Parallel to Level 4a)

Another fourth-level `tar_map()` creates targets for different event types (this is a **sibling** to Level 4a above, not nested):

```r
tar_map(
  tibble(
    event_type = c("right_censored", "uncensored"),
    event_slicer = c(\(d, n) d, \(d, n) slice_sample(d, n = n)),
  ),
  names = "event_type",
  ...
)
```

**Suffixes:** `_right_censored`, `_uncensored`

**Example targets:**
- `tumor_ssls_patient_states_rvar_right_censored_prior_ctdna_apr`
- `tumor_ssls_patient_states_rvar_uncensored_posterior_no_hist_aug`
- `tumor_ssls_rep_sld_rvar_right_censored_posterior_no_ctdna_aug`

**Note:** This branch generates targets for patient-level states, SLD values, and RECIST assessments.

---

## Complete Examples

### Example 1: Level 4a branch (patient-level parameters)

```
tumor_ssls_rates_bpi_patient_posterior_no_ctdna_aug
```

Breaking this down:
- `tumor_ssls_rates_bpi` - base target name
- `patient` - Level 4a: analysis level (from first sibling tar_map)
- `posterior` - Level 3: type (posterior fit)
- `no_ctdna` - Level 2: covariate group
- `aug25` - Level 1: DCO name

### Example 2: Level 4b branch (event type states)

```
tumor_ssls_patient_states_rvar_right_censored_posterior_no_hist_apr
```

Breaking this down:
- `tumor_ssls_patient_states_rvar` - base target name
- `right_censored` - Level 4b: event type (from second sibling tar_map)
- `posterior` - Level 3: type (posterior fit)
- `no_hist` - Level 2: covariate group
- `apr25` - Level 1: DCO name

## Aggregated Targets

Some targets aggregate across prior/posterior types and are prefixed with `all_`:

**Examples:**
- `all_tumor_ssls_orr_rvar` - combines `tumor_ssls_orr_rvar_prior` and `tumor_ssls_orr_rvar_posterior`
- `all_tumor_ssls_km_rvar` - combines KM estimates from both prior and posterior
- `all_tumor_ssls_patient_states_rvar` - combines patient states across event types and fit types

These targets include both:
1. Data from prior and posterior runs (via `bind_rows()`)
2. A `fit_type` column indicating which run
3. A `dco` column indicating the data cut-off

## Finding Targets

### Using `tar_manifest()`

To see all target names in a pipeline:

```r
tar_manifest() |> pull(name)
```

### Using `tar_visnetwork()`

To visualize the target dependency graph:

```r
tar_visnetwork()
```

### Filtering by Pattern

To find all targets for a specific configuration:

```r
# All August DCO targets
tar_manifest() |> filter(str_detect(name, "_aug$"))

# All no_ctdna covariate group targets
tar_manifest() |> filter(str_detect(name, "_no_ctdna_"))

# All posterior fit targets
tar_manifest() |> filter(str_detect(name, "tumor_ssls_.*_posterior_"))
```

## Common Target Patterns

### Data Preparation Targets
- `*_patient_data_file_*` - Raw patient data files
- `*_visit_data_file_*` - Raw visit data files
- `*_analysis_data_*` - Cleaned and prepared analysis datasets
- `all_analysis_data_*` - Combined data across trials

### Model Targets
- `*_stan_data_*` - Prepared Stan input data
- `*_res_*` - Model fit results
- `*_draws_*` - MCMC draws from the model
- `*_nuts_param_*` - NUTS diagnostic parameters

### Analysis Targets
- `*_rvar_*` - Random variable summaries
- `*_bpi_*` - Binned posterior intervals
- `*_km_*` - Kaplan-Meier estimates
- `*_orr_*` - Objective response rate
- `*_pfs_*` - Progression-free survival

### LFO (Leave-Future-Out) Targets
- `lfo_cutoffs_*` - Time points for cross-validation
- `tumor_ssls_lfo_*` - LFO cross-validation results
- `oos_confusion_matrix_*` - Out-of-sample confusion matrices

## Tips

1. **Start from the end:** If you know what you want (e.g., "posterior KM estimates for no_ctdna group in August"), work backwards through the suffixes.

2. **Use wildcards:** When loading targets, you can use patterns:
   ```r
   tar_read(starts_with("tumor_ssls_km_rvar_posterior"))
   ```

3. **Check dependencies:** Use `tar_deps()` to see what a target depends on:
   ```r
   tar_deps(tumor_ssls_res_posterior_no_ctdna_aug)
   ```

4. **Branch patterns:** Targets with `pattern = map()` create one target per row in the mapped data frame.

## Configuration Variables

Key variables that control the pipeline behavior:

- `covar_formula` - Full covariate formula with ctDNA
- `covar_formula_no_ctdna` - Covariate formula excluding ctDNA
- `covar_formula_ctdna_only` - Only ctDNA as covariate
- `lfo_groups` - Number of groups for parallel LFO computation (default: 5)
- `lfo_step` - Days between LFO cutoffs (default: 30)
- `extend_max_all_t` - Maximum time for forecasting (default: 200)
