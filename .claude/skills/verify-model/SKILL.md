---
name: verify-model
description: After Stan model changes, verify implementation matches mathematical specifications in docs. Use when stan/ files are modified or when user requests model verification.
user-invocable: true
allowed-tools: [Read, Grep, Glob, Bash]
---

# Model Verification Skill

Verifies that Stan model implementation matches the mathematical specifications documented in `docs/MODEL_MATHEMATICAL_SPECIFICATION.md` and related documentation.

## When to Use

- After modifications to files in `stan/` or `stan/modules/`
- When user explicitly requests model verification
- Before committing changes to Stan models
- When reviewing pull requests that touch Stan code

## Verification Process

### 1. Identify What Changed
- Determine which Stan files were modified
- Identify which module(s) are affected (tr, frac, init, other_events, measurement)

### 2. Locate Relevant Documentation
Read the corresponding sections from:
- `docs/MODEL_MATHEMATICAL_SPECIFICATION.md` - Primary mathematical specification
- `docs/ARCHITECTURE.md` - System architecture and parameter conventions
- `docs/ar1_process_noise.md` - If process noise is involved
- `docs/OTHER_EVENTS_MODEL.md` - If other events module is involved

### 3. Verify Implementation Against Specification

Check for alignment in:

**Parameter Names:**
- Population parameters match naming convention: `*_pop` (e.g., `tr_intercept_pop`)
- Trial-level standard deviations: `*_sd_trial_*` (e.g., `tr_sd_trial_intercept`)
- Patient-level standard deviations: `*_sd_patient_*` (e.g., `tr_sd_patient_intercept`)
- Raw effects use non-centered parameterization: `*_raw_trial_*`, `*_raw_patient_*`

**Mathematical Equations:**
- Compare Stan code formulations to documented equations
- Check for algebraically equivalent forms (e.g., `log(x * y)` vs `log(x) + log(y)`)
- Verify coefficient signs (positive/negative)
- Confirm functional forms (linear, exponential, etc.)

**Prior Distributions:**
- Match prior families (normal, lognormal, exponential, etc.)
- Verify hyperparameters align with documented priors
- Check that priors are specified in correct parameterization

**Hierarchical Structure:**
- Verify three-level hierarchy: population → trial → patient
- Confirm non-centered parameterization is used
- Check centering and scaling transformations

**Module Prefixes:**
- `tr_*` - Tumor regression parameters
- `frac_*` - Growth fraction parameters
- `init_*` - Initial state parameters
- `oe_*` - Other events parameters

### 4. Flag Potential Issues

Report any findings in this format:

**✅ VERIFIED:**
- List aspects that correctly match specification

**⚠️ POTENTIAL MISMATCHES:**
- Parameter name: Expected `X` from docs, found `Y` in code
- Equation: Docs show `Z`, code implements `W` (may be equivalent?)
- Prior: Docs specify `distribution(params)`, code uses `different_distribution(params)`

**❓ UNCLEAR:**
- List anything ambiguous that needs human judgment

**📍 REFERENCES:**
- Cite specific equation numbers, section headings, or line numbers from docs
- Cite specific file paths and line numbers from Stan code

### 5. Consider Equivalent Formulations

Recognize these as potentially valid:
- Log-space transformations
- Algebraic rearrangements
- Different but equivalent parameterizations
- Vectorized vs loop implementations

Do not flag these as mismatches unless there's a semantic difference.

## Limitations

This skill provides best-effort verification but cannot:
- Guarantee mathematical equivalence (may require manual proof)
- Detect all subtle semantic bugs
- Verify numerical stability or computational efficiency

## Examples

### Example 1: Tumor Regression Module Change
```
User modified: stan/modules/tr/parameters.stan

1. Read docs/MODEL_MATHEMATICAL_SPECIFICATION.md Section 3.1
2. Compare parameter definitions
3. Verify naming: tr_intercept_pop, tr_sd_trial_intercept, etc.
4. Check that regression equation matches documented form
5. Report findings
```

### Example 2: Prior Specification Change
```
User modified: stan/modules/frac/priors.stan

1. Read docs/MODEL_MATHEMATICAL_SPECIFICATION.md Section 3.2
2. Check documented prior distributions
3. Compare to implemented priors in Stan code
4. Verify hyperparameter values match r/priors.R defaults
5. Report findings
```

## Related Files

- `r/priors.R` - Default prior hyperparameters (should match docs)
- `r/initializers.R` - Initialization values (should respect model structure)
- `targets/sclc_targets.R` - Pipeline configuration (flags and settings)
