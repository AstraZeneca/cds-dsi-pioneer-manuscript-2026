# Fixed Effects Level Intercepts for Multistate Module

**Date**: 2026-03-25
**Branch**: `karim/burden-endpoints`
**Scope**: Multistate module only (`stan/modules/multistate/`)

## Problem

The multistate module's level baseline hazard uses random effects (RE) with an estimated
hierarchical SD parameter (`level_intercept_sd[lv]`). With only 2 groups at the trial level
(n=107 trial vs n=1,666 RWD), the estimated SD shrinks toward zero, pulling the trial
intercept toward the RWD-dominated population mean. This causes the superpopulation CIF for
the trial arm to inherit the RWD's event patterns (dropout > progression) even though the
trial data shows the opposite (progression 35% > dropout 21%).

Fixed effects (FE) would give each group an independent intercept with a data-specified prior
SD, eliminating pooling and letting the trial arm reflect its own data.

## Design

### Flag Encoding

`enable_ms_level_baseline_hazard` remains `array[n_levels] int` but the encoding changes:

| Value | Mode | Intercept SD source | Pooling | GP kernel |
|-------|------|---------------------|---------|-----------|
| 0 | Off | — | — | No |
| 1 | FE intercept (new) | Data hyperparameter | None | No |
| 2 | RE intercept | Estimated parameter | Yes | No |
| 3 | RE GP | Estimated parameter | Yes | Yes |

Upper bound in `flags.stan` changes from `2` to `3`. Update the comment at line 21 of
`flags.stan` to reflect the new 4-state encoding.

### Encoding Migration

All existing configs must update:

- Old `1` (RE intercept) → new `2`
- Old `2` (RE GP) → new `3`
- `0` stays `0`

Affected files:
- `targets/sclc_targets.R` (two sites: lines ~977 and ~2320)
- `targets/pioneer_targets.R` (line ~412)

### New Data Hyperparameters

Per-transition FE prior SD in `hyperparams.stan`:

```stan
array[n_levels] real<lower=0> fe_log_lambda_gp_01_level_intercept_sd;
array[n_levels] real<lower=0> fe_log_lambda_gp_02_level_intercept_sd;
array[n_levels] real<lower=0> fe_log_lambda_gp_12_s_level_intercept_sd;
array[n_levels] real<lower=0> fe_log_lambda_gp_12_t_level_intercept_sd;
array[n_levels] real<lower=0> fe_log_lambda_gp_03_level_intercept_sd;
array[n_levels] real<lower=0> fe_log_lambda_gp_32_s_level_intercept_sd;
```

These are only read when the corresponding level's flag == 1. Default values set in
`r/priors.R`.

### Transformed Data

Add a per-transition boolean `any_re_level_XX` computed as:

```stan
int any_re_level_01 = enable_ms_01
  ? (max(to_array_1d(enable_ms_level_baseline_hazard)) >= 2 ? 1 : 0)
  : 0;
```

This gates the conditional sizing of RE parameters. Note: the ternary guards
against calling `max()` on a zero-length array when the transition is disabled.

**`compute_n_enabled_groups` and `create_enabled_pos`**: No changes needed. These use
`if (enabled[i])` which is truthy for any nonzero value (1, 2, or 3), so FE groups are
already included in the enabled set.

### GP Threshold Migration (`== 2` → `== 3`)

The GP detection mask at `transformed_data.stan` line 62 must change:

```stan
// Before:
ms_level_baseline_is_gp[lv] = (enable_ms_level_baseline_hazard[lv] == 2) ? 1 : 0;
// After:
ms_level_baseline_is_gp[lv] = (enable_ms_level_baseline_hazard[lv] == 3) ? 1 : 0;
```

This array gates all `n_gp_groups_ms_baseline_*` counts and `gp_level_pos_ms_baseline`.
Without this change, flag=3 would produce zero-sized GP eta matrices.

**Complete list of `== 2` → `== 3` threshold changes:**

| File | Lines | What |
|------|-------|------|
| `multistate/transformed_data.stan` | ~62 | `ms_level_baseline_is_gp` |
| `multistate/transformed_parameters.stan` | ~44, ~181, ~295, ~384, ~480, ~546 | GP branch in each transition block |
| `r/initializers_ms.R` | ~30 | `n_gp_groups_ms_baseline` computation |
| `r/initializers_fixed.R` | ~27 | Same |

### Parameters

The RE SD parameter becomes conditionally sized:

```stan
// Only declare when at least one level uses RE mode (flag >= 2)
array[enable_ms_01 && any_re_level_01 ? n_levels : 0]
  real<lower=0> log_lambda_gp_01_level_intercept_sd;
```

The raw intercept vector is unchanged — used by both FE and RE:

```stan
vector[n_enabled_groups_ms_baseline_01] raw_log_lambda_gp_01_level_intercept;
```

### Transformed Parameters

The scaling step becomes conditional on the flag value:

```stan
for (lv in 1:n_levels) {
  if (enable_ms_level_baseline_hazard[lv] >= 1) {
    int lv_start, lv_end;
    (lv_start, lv_end) = get_pos(enabled_level_pos_ms_baseline, lv);

    if (enable_ms_level_baseline_hazard[lv] == 1) {
      // FE: scale by fixed data hyperparameter
      log_lambda_gp_01_level_intercept[lv_start:lv_end] =
        raw_log_lambda_gp_01_level_intercept[lv_start:lv_end] *
        fe_log_lambda_gp_01_level_intercept_sd[lv];
    } else {
      // RE (2 or 3): scale by estimated parameter
      log_lambda_gp_01_level_intercept[lv_start:lv_end] =
        raw_log_lambda_gp_01_level_intercept[lv_start:lv_end] *
        log_lambda_gp_01_level_intercept_sd[lv];
    }

    if (enable_ms_level_baseline_hazard[lv] == 3) {
      // GP residual (unchanged from current flag==2 logic)
      ...
    } else {
      // Intercept-only (FE=1 or RE=2): constant shift
      ...
    }
  }
}
```

This pattern repeats for all 6 transition blocks (01, 02, 12_s, 12_t, 03, 32_s).

### Priors

```stan
for (lv in 1:n_levels) {
  if (enable_ms_level_baseline_hazard[lv] >= 2) {
    // RE prior on estimated SD (only for RE modes 2 and 3)
    log_lambda_gp_01_level_intercept_sd[lv] ~ normal(
      0, log_lambda_gp_01_level_intercept_sd_sd[lv]);
  }
}
// Raw prior unchanged — applies to both FE and RE
raw_log_lambda_gp_01_level_intercept ~ std_normal();
```

### R Pipeline Changes

#### `r/priors.R`

Add default FE prior SDs per transition. Reasonable default: 1.0 (weakly informative on
log-hazard scale — allows ~2.7× hazard ratio per group).

```r
fe_log_lambda_gp_01_level_intercept_sd = rep(1.0, n_levels),
fe_log_lambda_gp_02_level_intercept_sd = rep(1.0, n_levels),
# ... etc for all transitions
```

#### `r/initializers_ms.R` and `r/initializers_fixed.R`

Conditionally skip `level_intercept_sd` initialization when no RE levels are active:

```r
if (any(enable_ms_level_baseline_hazard >= 2L)) {
  init$log_lambda_gp_01_level_intercept_sd <- ...
}
```

Update GP group count threshold from `== 2L` to `== 3L`:

```r
# Before:
n_gp_groups_ms_baseline <- sum(n_groups_per_level[enable_ms_level_baseline_hazard == 2L])
# After:
n_gp_groups_ms_baseline <- sum(n_groups_per_level[enable_ms_level_baseline_hazard == 3L])
```

#### `targets/pioneer_targets.R`

Migrate encoding:
```r
# Before: c(trial = TRUE, arm = FALSE, patient = FALSE)  →  c(1, 0, 0)
# After:  c(trial = 2L, arm = 0L, patient = 0L)  for RE intercept (same behavior)
# Or:     c(trial = 1L, arm = 0L, patient = 0L)  for FE intercept (new)
enable_ms_level_baseline_hazard = as.array(c(trial = 1L, arm = 0L, patient = 0L)),
```

#### `targets/sclc_targets.R`

Migrate encoding at both sites:
```r
# Line ~977: c(trial = 2L, patient = 0L) → c(trial = 3L, patient = 0L)  (was GP)
# Line ~2320: c(trial = 1L, patient = 0L) → c(trial = 2L, patient = 0L)  (was RE intercept)
```

### Laplace Module

`stan/modules/laplace/model.stan` references `enable_ms_level_baseline_hazard` at lines
~77, ~120, ~171. These check `if (enable_ms_level_baseline_hazard[lv])` which remains
true for all values ≥ 1. The Laplace code uses `patient_ms_baseline_flat_idx` for routing,
which is computed in `laplace/transformed_data.stan` based on the same flag. Both handle
FE correctly because the group indexing is unchanged — only the scaling source differs,
which happens in `multistate/transformed_parameters.stan` before Laplace sees the values.

No Laplace module changes needed.

### Stan Syntax Check

After all changes, verify with:

```bash
~/.cmdstan/cmdstan-2.38.0/bin/stanc --include-paths=stan --include-paths=stan/tumor \
  stan/tumor/sf-ssm-log-space.stan

~/.cmdstan/cmdstan-2.38.0/bin/stanc --include-paths=stan --include-paths=stan/psa \
  stan/psa/pioneer.stan

~/.cmdstan/cmdstan-2.38.0/bin/stanc --include-paths=stan --include-paths=stan/psa \
  stan/psa/ms-standalone.stan
```

### Testing

1. **Syntax check**: Both Stan models compile without errors
2. **FE mode smoke test**: Set `enable_ms_level_baseline_hazard = c(1L, 0L, 0L)` in
   pioneer, run standalone posterior combined, verify no divergences and that trial-level
   intercepts are not pooled toward zero
3. **RE backward compatibility**: Set `c(2L, 0L, 0L)`, verify identical behavior to the
   old `c(1L, 0L, 0L)` encoding
4. **GP backward compatibility**: Set `c(3L, 0L, 0L)` in sclc, verify identical to
   old `c(2L, 0L, 0L)`
