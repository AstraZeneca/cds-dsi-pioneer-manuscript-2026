---
name: stan-flags
description: Display enable_* flags from Stan data targets. Use to check model configuration flags for SCLC pipeline runs. Supports comparing across branches and models.
user-invocable: true
allowed-tools: [mcp__targets__tar_read, mcp__targets__tar_progress, mcp__targets__eval_expr, mcp__targets__list_pipelines, Bash]
---

# Stan Flags Skill

Display `enable_*` flags from Stan data targets to quickly check model configuration for SCLC pipeline runs.

## Usage

```
/stan-flags [target1] [target2] ...
```

### Target Notation

Use `branch::model` notation to specify targets:
- **branch**: Pipeline store name (e.g., `test2`, `test3`, `main`)
- **model**: Model type from `sclc_targets.R` (e.g., `ctdna`, `no_oe`, `no_ctdna`)

### Show Hints

When invoked with no arguments or invalid format, display usage hints:

```
Usage: /stan-flags [target1] [target2] ...

Format: branch::model

Available models:
  ctdna, no_oe, no_ctdna, ctdna_only, all_trials, no_covar, no_hist, no_pl

Available branches: [detect from /mnt/data/analysis-results/$DOMINO_STARTING_USERNAME/sclc/]
  test2 (2 models), test3 (1 model), main (3 models)

Examples:
  /stan-flags test2::ctdna              # Single target
  /stan-flags test2                     # All models in test2
  /stan-flags test2::ctdna test3::ctdna # Compare across branches
```

### Examples

```bash
# Check specific branch and model
/stan-flags test2::ctdna

# Check all built models in a branch
/stan-flags test2

# Compare same model across branches
/stan-flags test2::ctdna test3::ctdna

# Compare different models in same branch
/stan-flags test2::ctdna test2::no_oe

# Compare multiple configurations
/stan-flags test2::ctdna test3::ctdna test2::no_oe

# Show available targets (no arguments)
/stan-flags
```

## Model Types

From `targets/sclc_targets.R`, available models are:
- `ctdna` - With covariates and other events
- `no_oe` - Without other events
- `no_ctdna` - Without ctDNA covariates
- `ctdna_only` - Only ctDNA covariates
- `all_trials` - All trials included
- `no_covar` - Without covariates
- `no_hist` - Without historical controls
- `no_pl` - Without patient level effects

## Implementation

### 1. Parse Arguments and Show Hints

**If no arguments provided:**
- Show usage hints with available models and branches
- List branches by scanning `/mnt/data/analysis-results/$DOMINO_STARTING_USERNAME/sclc/`
- For each branch, use `tar_progress` to count completed stan_data targets
- Display examples

**If arguments provided:**
Parse `branch::model` notation:
```
test2::ctdna → branch="test2", model="ctdna"
test2 → branch="test2", model=null
```

### 2. Construct Store Path

Store path pattern:
```
/mnt/data/analysis-results/$DOMINO_STARTING_USERNAME/sclc/{branch}/_targets
```

### 4. Determine Targets

**If model specified:**
- Target name: `tumor_ssls_stan_data_posterior_{model}_aug`

**If model not specified:**
- Use `mcp__targets__tar_progress(store)` to list all targets
- Filter for targets matching `tumor_ssls_stan_data_posterior_*_aug`
- Extract model name from target name pattern
- Only include targets with `progress="completed"`

### 5. Read Stan Data

For each target:
```
mcp__targets__tar_read(
  name = "tumor_ssls_stan_data_posterior_{model}_aug",
  store = store_path
)
```

### 6. Extract Enable Flags

Use `mcp__targets__eval_expr` to extract flags:
```r
stan_data <- targets::tar_read(
  tumor_ssls_stan_data_posterior_{model}_aug,
  store = "{store_path}"
)
stan_data[grep("^enable_", names(stan_data))]
```

### 7. Format Output

**Single Target:**
```
## Stan Flags: test2::ctdna

Tumor Regression (tr):
  enable_patient_intercept_tr: TRUE
  enable_patient_process_noise_tr: TRUE
  enable_trial_intercept_tr: FALSE
  ...

Growth Fraction (frac):
  enable_patient_intercept_frac: TRUE
  enable_pop_cov_frac: TRUE
  ...

Initial State (init):
  enable_patient_intercept_init: TRUE
  enable_pop_cov_init: TRUE
  ...
```

**Multiple Targets (Comparison):**
```
## Stan Flags Comparison

| Flag | test2::ctdna | test3::ctdna |
|------|--------------|--------------|
| **Tumor Regression (tr)** | | |
| enable_patient_intercept_tr | TRUE | TRUE |
| enable_patient_process_noise_tr | TRUE | FALSE ⚠️ |
| enable_trial_intercept_tr | FALSE | FALSE |
| **Growth Fraction (frac)** | | |
| enable_patient_intercept_frac | TRUE | TRUE |
| enable_pop_cov_frac | TRUE | TRUE |
...

⚠️ = Differences highlighted
```

**No Arguments (Show Hints):**
```
Usage: /stan-flags [target1] [target2] ...

Format: branch::model

Available models:
  ctdna, no_oe, no_ctdna, ctdna_only, all_trials, no_covar, no_hist, no_pl

Available branches:
  test2 (2 models: ctdna, no_oe)
  test3 (1 model: ctdna)
  main (2 models: ctdna, no_ctdna)

Examples:
  /stan-flags test2::ctdna              # Single target
  /stan-flags test2                     # All models in test2
  /stan-flags test2::ctdna test3::ctdna # Compare across branches
```

## Error Handling

**Branch not found:**
```
Error: Branch 'test99' not found
Available branches: test2, test3, main
```

**Model not built:**
```
Error: Model 'no_ctdna' not found in branch 'test2'
Available models in test2: ctdna, no_oe
```

**Target not completed:**
```
Error: Target 'tumor_ssls_stan_data_posterior_ctdna_aug' in test2 is 'dispatched' (not completed)
```

**No stan data targets found:**
```
No stan data targets found in branch 'test2'
Run tar_progress to check pipeline status.
```

## Flag Categories

Group flags by module for better readability:

**Tumor Regression (tr):**
- `enable_patient_intercept_tr`
- `enable_trial_intercept_tr`
- `enable_pop_cov_tr`
- `enable_trial_cov_tr`
- `enable_patient_cov_tr`
- `enable_patient_process_noise_tr`
- `enable_patient_process_noise_sd_tr`
- `enable_patient_process_noise_phi_tr`

**Growth Fraction (frac):**
- `enable_patient_intercept_frac`
- `enable_trial_intercept_frac`
- `enable_pop_cov_frac`
- `enable_trial_cov_frac`
- `enable_patient_cov_frac`

**Initial State (init):**
- `enable_patient_intercept_init`
- `enable_trial_intercept_init`
- `enable_pop_cov_init`
- `enable_trial_cov_init`
- `enable_patient_cov_init`

## Related Files

- `targets/sclc_targets.R` - Model definitions and configurations
- `r/priors.R` - Prior specifications
- `stan/modules/*/flags.stan` - Stan flag definitions
