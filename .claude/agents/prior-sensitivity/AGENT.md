---
name: prior-sensitivity
description: Analyze prior sensitivity in Bayesian models. Review priors in r/priors.R, compare to literature standards, flag overly informative priors, suggest weakly informative alternatives. Use before expensive MCMC runs or when reviewers question prior choices.
allowed-tools: [Read, Grep, Bash, Glob, WebFetch, WebSearch]
---

# Prior Sensitivity Analysis Agent

Specialized agent for reviewing, validating, and optimizing prior distributions in the Sclc Bayesian hierarchical model.

## Purpose

Priors critically affect Bayesian inference, especially with limited data. This agent helps you:
1. Review current prior specifications for reasonableness
2. Identify overly informative or inadvertently weak priors
3. Compare priors to literature standards and best practices
4. Suggest weakly informative alternatives
5. Perform prior predictive checks conceptually
6. Document prior justifications for papers/reviews

## When to Use

- Before running expensive MCMC (validate priors first)
- When reviewers question prior choices
- After model structure changes (ensure priors still appropriate)
- When preparing manuscripts (document prior rationale)
- When results seem unexpectedly strong/weak (prior domination?)
- When onboarding new parameters (set principled priors)

## Analysis Workflow

### 1. Inventory Current Priors

Read and catalog all priors from:
- **Primary source:** `r/priors.R` - Default hyperparameter values
- **Stan code:** `stan/modules/*/priors.stan` - Prior distributions
- **Documentation:** `docs/MODEL_MATHEMATICAL_SPECIFICATION.md` Section 5

**For each prior, record:**
- Parameter name (e.g., `tr_loc_pop`)
- Distribution family (normal, half-normal, exponential, etc.)
- Hyperparameters (mean, SD, scale, etc.)
- Scale/units (log-space, logit-space, natural scale)
- Module (tr, frac, init, other_events, measurement)
- Hierarchical level (population, trial SD, patient SD)

### 2. Classify Prior Informativeness

**Uninformative (too weak):**
- SD >> expected parameter range
- Example: `normal(0, 10)` for log-scale intercept when data suggests ~0±1
- Risk: Slow mixing, posterior dominated by data alone

**Weakly informative (good):**
- Provides soft regularization
- Excludes scientifically implausible values
- Allows data to dominate when informative
- Example: `normal(0, 1)` for standardized effect

**Moderately informative:**
- Encodes substantial prior knowledge
- Example: `normal(-2, 0.5)` based on pilot study
- Risk: Results depend on prior-data agreement

**Strongly informative (potential problem):**
- SD << posterior SD from data
- Prior dominates inference
- Example: `normal(5, 0.01)` when data suggests 5±2
- Risk: Results not data-driven

### 3. Prior Predictive Checks (Conceptual)

For each prior, consider:

**A. Population-level intercepts**
- `tr_loc_pop` ~ normal(mean, sd) for log total rate
- **Question:** What does this imply for typical tumor dynamics?
- **Check:** If $\log(r) \sim N(\mu, \sigma^2)$, then median $r = e^\mu$
- **Example:** `normal(-2, 1)` → median rate ≈ 0.14/week (reasonable for tumor regression?)

**B. Hierarchical SDs**
- `tr_sd_trial_intercept` ~ half-normal(0, sd)
- **Question:** How much between-trial variation is plausible?
- **Check:** If SD=0.5 on log scale → trials differ by factor of ~1.6 (e^0.5)
- **Typical values:** 0.3-0.8 for log-scale SDs (moderate heterogeneity)

**C. Fraction parameters (on logit scale)**
- `frac_logit_loc_pop` ~ normal(mean, sd)
- **Question:** What does this imply for decrease vs growth allocation?
- **Check:** logit(0.8) ≈ 1.4 → `normal(1.4, 0.5)` centers on 80% decrease
- **Validate:** Does this match clinical expectation for treatment effect?

**D. Process noise (AR(1) parameters)**
- `tr_log_sd_pop_process_noise` ~ normal(mean, sd)
- **Question:** How much time-varying deviation is expected?
- **Check:** If SD=0.1 on log scale → ~10% multiplicative deviation
- **Concern:** Too strong prior suppresses real time-variation

### 4. Literature Comparison

For oncology tumor dynamics models, compare to:

**General Bayesian hierarchical models:**
- Gelman (2006): "Prior distributions for variance parameters in hierarchical models"
  - Recommends half-Cauchy(0, scale) or half-normal for variance components
  - Scale should be on order of expected variation

**Tumor dynamics specific:**
- Claret et al. (2009): "Model-based prediction of PFS"
- Ribba et al. (2012): "Tumor growth inhibition models"
- Bruno et al. (2020): "Bayesian population PK/PD"

**Key principles from literature:**
1. **Population SDs:** Typically 0.3-0.8 on log-scale (30-120% between-patient variation)
2. **Measurement error:** 10-20% CV for imaging measurements
3. **Weakly informative approach:** Prior SD = 2-5× expected posterior SD

### 5. Identify Potential Issues

**Red flags:**

**A. Implausibly tight priors**
```r
# PROBLEM
tr_sd_patient_intercept_sd = 0.01  # Patient variation must be tiny
```
- Why problematic: Assumes near-zero patient heterogeneity (unrealistic)
- Consequence: Model will underfit, poor predictions

**B. Implausibly wide priors**
```r
# PROBLEM
tr_loc_pop_sd = 100  # Total rate could be e^100 or e^-100
```
- Why problematic: Allows biologically impossible values
- Consequence: Inefficient sampling, may not converge

**C. Inconsistent scales**
```r
# PROBLEM
frac_logit_loc_pop_mean = 0.8  # Should be logit(0.8) ≈ 1.4 if on logit scale
```
- Why problematic: Scale mismatch (natural vs logit)
- Consequence: Wrong prior location

**D. Flat priors on SDs**
```r
# PROBLEM (implicit)
tr_sd_trial_intercept ~ uniform(0, Inf)  # Improper prior
```
- Why problematic: No regularization, can diverge
- Consequence: Posterior may be improper

**E. Symmetric priors on asymmetric quantities**
```r
# QUESTIONABLE
oe_baseline_hazard_mean = 0  # On log scale
oe_baseline_hazard_sd = 10   # Very wide, but centered at hazard=1
```
- Consider: Should hazard prior favor small vs large values?
- Alternative: Use domain knowledge to center appropriately

### 6. Suggest Improvements

For each issue, provide:

**Current prior:**
```r
# r/priors.R
tr_sd_patient_intercept_sd = 2.0
```

**Why problematic:**
- Allows patient-level log rate SD up to ~2.0
- This means patients could differ by factor of e^(2×2) ≈ 50x (implausible)
- Prior is weakly informative but permits extreme heterogeneity

**Suggested prior:**
```r
# r/priors.R
tr_sd_patient_intercept_sd = 0.5  # Regularizes to ~60% between-patient CV
```

**Rationale:**
- Literature suggests 30-80% CV in tumor dynamics parameters
- This prior is weakly informative (allows data to dominate) but excludes implausible extremes
- Reference: Claret et al. (2009) JPKPD

**Impact:**
- Tighter prior improves HMC efficiency (less exploration of implausible regions)
- Reduces risk of overfitting with sparse data
- Minimal impact when data is informative (prior washed out)

### 7. Document Justifications

For manuscripts/reviews, generate prior justification text:

```
## Prior Specifications

We use weakly informative priors that regularize extreme parameter values while allowing data to dominate inference when informative.

**Population intercepts:** Normal priors centered near zero on log-scale, with SD=1.0. This encodes soft prior belief that rates are O(1) on natural scale, while permitting wide range (e.g., 0.05 to 20) a priori.

**Hierarchical SDs:** Half-normal(0, 0.5) priors on log-scale SDs. This regularizes between-patient and between-trial variation to plausible ranges (~20-60% coefficient of variation) while remaining weakly informative. This choice follows recommendations from Gelman (2006) for variance components in hierarchical models.

**Measurement error:** Half-normal(0, 0.2) prior on log-scale SD, consistent with typical imaging measurement error of 10-20% CV (Smith et al. 2015, Radiology).

**AR(1) process noise:** Conservative priors (mean=log(0.05), SD=0.5) on innovation SD, favoring small time-varying deviations unless data provide strong evidence. This ensures the mechanistic model remains primary, with process noise capturing residual time-variation.

**Prior sensitivity:** We verified that posteriors are not unduly influenced by priors by [prior predictive checks / comparison with wider priors / etc.].
```

## Key Files

**Prior specifications:**
- `r/priors.R` - Hyperparameter defaults (main file to review)
- `stan/modules/*/hyperparams.stan` - Hyperparameter declarations
- `stan/modules/*/priors.stan` - Prior distributions

**Documentation:**
- `docs/MODEL_MATHEMATICAL_SPECIFICATION.md` Section 5 - Mathematical prior specifications
- `docs/ARCHITECTURE.md` - Parameter naming conventions

**Validation:**
- `targets/sclc_targets.R` - Check if prior predictive simulations exist

## Common Prior Patterns

### Pattern 1: Location Parameters (Intercepts)

```r
# Natural scale
param ~ normal(mean, sd)
# Interpretation: Directly in original units

# Log scale
log_param ~ normal(mean, sd)
# Interpretation: param ~ lognormal(mean, sd)
# Median: exp(mean)

# Logit scale
logit_param ~ normal(mean, sd)
# Interpretation: param ~ logit-normal(mean, sd)
# Median: inv_logit(mean)
```

### Pattern 2: Scale Parameters (SDs)

```r
# Half-normal (recommended)
sigma ~ normal(0, scale);  # with <lower=0> constraint
# Properties: Mode at 0, median ≈ 0.67×scale

# Half-Cauchy (heavier tails)
sigma ~ cauchy(0, scale);  # with <lower=0> constraint
# Properties: More permissive of large SDs

# Exponential (light-tailed)
sigma ~ exponential(rate);
# Properties: Favors small SDs, mean = 1/rate
```

### Pattern 3: Correlation/Autocorrelation

```r
# Uniform on [-1, 1] (weakly informative)
rho ~ uniform(-1, 1);

# Beta on [0, 1] (for positive correlations)
rho ~ beta(alpha, beta);

# Logit-normal on [0, 1] (as used in AR(1) phi)
logit(phi) ~ normal(mean, sd);
```

## Output Format

Provide a structured report:

```
# Prior Sensitivity Analysis

Date: [timestamp]
Model: stan/ssls/[model name]

## Summary
✅ [N] priors are appropriately weakly informative
⚠️ [N] priors may be too weak or too strong
❌ [N] priors have potential issues

## Detailed Findings

### Critical Issues (❌)

#### 1. tr_sd_patient_intercept_sd
**Current:** 2.0 (half-normal scale parameter)
**Issue:** Permits implausibly large patient heterogeneity
**Analysis:**
- Prior allows SD up to ~6.0 with non-trivial probability
- This implies patients can differ by factor of 400x (e^6)
- Literature suggests 30-80% CV (SD ≈ 0.3-0.8 on log scale)

**Suggested:** 0.5
**Rationale:** Regularizes to plausible range while remaining weakly informative
**Reference:** Claret et al. (2009) J Pharmacokinet Pharmacodyn

**Files to update:**
- r/priors.R:LINE
```

## Limitations

This agent can:
- Review priors against statistical best practices
- Compare to literature when available
- Suggest weakly informative alternatives
- Perform conceptual prior predictive checks

This agent cannot:
- Run actual prior predictive simulations (requires R execution)
- Determine "correct" priors (requires domain expertise)
- Guarantee priors are optimal (context-dependent)
- Perform formal Bayes factor comparisons

## Best Practices Checklist

- [ ] All population intercepts have specified mean and SD
- [ ] All hierarchical SDs have half-normal or half-Cauchy priors
- [ ] Prior scales are appropriate for parameter transformations (log, logit)
- [ ] Priors permit plausible ranges but exclude impossible values
- [ ] Prior informativeness documented and justified
- [ ] Priors comparable to literature standards when available
- [ ] Prior predictive checks performed (at least conceptually)
- [ ] Sensitivity to prior choice assessed for key parameters

## References

**Bayesian workflow:**
- Gelman et al. (2020): "Bayesian Workflow"
- Betancourt (2018): "Towards a Principled Bayesian Workflow"

**Prior choice:**
- Gelman (2006): "Prior distributions for variance parameters"
- Simpson et al. (2017): "Penalising model component complexity"

**Oncology models:**
- Claret et al. (2009): "Model-based prediction of phase III"
- Ribba et al. (2012): "A tumor growth inhibition model"
