# Propensity-Weighted Bayesian Borrowing from RWD

## Motivation

The pioneer model combines clinical trial data (FPI-2265-202, ~107 patients) with
real-world data (Flatiron, ~2,047 patients). The question is: how much should RWD patients
influence the posterior?

Previously we used a **trial-level random effect** (n=2 groups) to control borrowing, but:
- With n=2 trials, the between-trial SD is essentially unidentified from data
- The funnel geometry of the trial-level hierarchy slows HMC dramatically
- It doesn't respect covariate overlap — an RWD patient identical to a trial patient
  is treated the same as one with completely different characteristics

## Core Idea

Replace the trial-level hierarchy with a **joint propensity submodel** that learns
`P(would be in trial | covariates)` and uses it to weight each RWD patient's
likelihood contribution. Trial patients always contribute fully (weight = 1).
RWD patients contribute proportionally to how "trial-like" they are.

## Design: Joint Bayesian Propensity Model

### Why joint (not two-stage)?

A two-stage approach (fit propensity model in R → pass fixed weights to Stan) ignores
uncertainty in the propensity scores. The joint approach:
- Propagates propensity uncertainty to endpoint estimates
- Allows the PSA/survival data to inform the propensity model (feedback)
- Avoids arbitrary point-estimate cutoffs
- Is fully Bayesian — one model, one posterior

### Stan Architecture

**New data field:**
```stan
// In _base_data.stan or a new modules/propensity/data.stan:
array[n_patients] int<lower=0, upper=1> is_trial;
```

**New parameters:**
```stan
// ~17 new parameters — logistic regression is convex, trivial for HMC
vector[n_covar] beta_propensity;
```

**Transformed parameters — compute weights:**
```stan
vector<lower=0, upper=1>[n_patients] likelihood_weight;
{
  vector[n_patients] logit_propensity = covar_design_matrix * beta_propensity;
  for (i in 1:n_patients) {
    likelihood_weight[i] = is_trial[i] ? 1.0 : inv_logit(logit_propensity[i]);
  }
}
```

**Model block — propensity likelihood + weighted outcome likelihoods:**
```stan
// Propensity submodel
beta_propensity ~ normal(0, 2.5);  // weakly informative
is_trial ~ bernoulli_logit(covar_design_matrix * beta_propensity);

// PSA likelihood — weighted per patient
for (j in 1:n_hmc_patients) {
  int p = hmc_patient_idx[j];
  target += likelihood_weight[p] * psa_loglik(p);
}

// Multistate likelihood — weighted per patient
// (requires adding weight argument to multistate_lpmf)
```

### Changes Required

| Component | Change | Complexity |
|-----------|--------|------------|
| New module `modules/propensity/` | `data.stan`, `parameters.stan`, `transformed_parameters.stan`, `priors.stan` | New, small |
| `pioneer.stan` | Include propensity module, add weighted PSA likelihood | Medium |
| `multistate.stanfunctions` | Add `vector weight` argument to `multistate_lpmf`, multiply `patient_ll` by `weight[i]` | Medium |
| `r/pioneer/prepare_analysis_data.R` | Add `is_trial = as.integer(trial == "FPI-2265-202")` to stan data | Trivial |
| `r/pioneer/initializers.R` | Add `beta_propensity = rep(0, n_covar)` | Trivial |
| `targets/pioneer_targets.R` | Add `combined_no_trial` model variant (already done on `burden-endpoints` branch) | Done |

### PSA Likelihood Weighting

The PSA likelihood loop in `pioneer.stan` is already per-patient:
```stan
for (j in 1:n_hmc_patients) {
  int p = hmc_patient_idx[j];
  normalized_psa[data_start + n_screen : data_end] ~ sf_log_space_obs(
    states[state_start + n_screen : ...], measure_sd_psa, ...);
}
```

The `~` sampling statement adds to `target`. To weight it, convert to:
```stan
target += likelihood_weight[p] * sf_log_space_obs_lpdf(
    normalized_psa[...] | states[...], measure_sd_psa, ...);
```

### Multistate Likelihood Weighting

`multistate_lpmf` currently accumulates a scalar sum across all patients:
```stan
for (i in 1:n_patients) {
  real patient_ll = 0;
  // ... compute patient_ll ...
  total_ll += patient_ll;
}
return total_ll;
```

Add a weight vector parameter:
```stan
real weighted_multistate_lpmf(
  array[] int final_state,
  vector weight,    // NEW: per-patient weight
  // ... existing args ...
) {
  for (i in 1:n_patients) {
    real patient_ll = 0;
    // ... existing logic unchanged ...
    total_ll += weight[i] * patient_ll;  // CHANGED: weighted
  }
  return total_ll;
}
```

### What the Weights Do NOT Affect

- **Priors** — the hierarchical priors on population/arm/patient parameters are unweighted.
  The propensity weight controls how much each patient's DATA influences the posterior,
  not the prior shrinkage structure.
- **Generated quantities** — endpoint computations (PFS, KM, OS) use ALL patients unweighted.
  The weights only affect parameter estimation.

## Design Choices to Make

### 1. Weight function

The simplest is `w_i = P(trial | x_i)` via logistic regression. Alternatives:
- **Overlap weights**: `w_i = P(trial | x_i) * (1 - P(trial | x_i))` — upweights patients
  in the overlap region, downweights extremes on both sides
- **Trimmed weights**: `w_i = max(epsilon, P(trial | x_i))` — floor to prevent zero-weight patients
- **Power parameter**: `w_i = P(trial | x_i)^alpha` where `alpha` is learned — controls
  how aggressively dissimilar patients are downweighted

### 2. Should the propensity model use the QR-transformed covariates?

The PSA and multistate models use QR-decomposed covariates. The propensity model could:
- Use the same QR-transformed matrix (shares the decomposition)
- Use the original covariate matrix (more interpretable coefficients)

Using the original is simpler and more interpretable. The propensity model is just logistic
regression — it doesn't need QR for numerical stability.

### 3. Feature flag

Add `enable_propensity_weighting` flag. When off, `likelihood_weight = ones_vector(n_patients)`.
This allows A/B comparison of weighted vs unweighted fits.

### 4. Standalone model

The standalone multistate model could also use propensity weighting. It already has
`covar_design_matrix` in its data. Adding the propensity submodel there would be
straightforward.

## Relationship to Other Approaches

| Approach | Trial-level RE | Propensity weights | Laplace | Status |
|----------|:-:|:-:|:-:|---|
| `combined` (current) | ✓ | ✗ | ✗ | Running (job #267), slow warmup |
| `combined` + Laplace | ✓ | ✗ | ✓ | Running (job #279), stuck at init |
| `combined_no_trial` | ✗ | ✗ | ✗ | Running (job #281) |
| `combined_no_trial` + Laplace | ✗ | ✗ | ✓ | Running (job #280), stuck at init |
| **Propensity-weighted** (this) | ✗ | ✓ | ✗ | To implement |
| **Propensity + Laplace** | ✗ | ✓ | ✓ | Future |

The propensity approach replaces the trial-level RE with something more principled AND
potentially resolves the Laplace initialization issues (by downweighting dissimilar RWD
patients whose extreme covariates cause QR overflow).

## Implementation Plan

1. Create `modules/propensity/` with the standard 4-file module pattern
2. Add `is_trial` to R data prep and Stan data
3. Add `beta_propensity` parameter + `likelihood_weight` computation
4. Convert PSA likelihood from `~` to `target +=` with weight
5. Add `weight` argument to `multistate_lpmf`
6. Add `enable_propensity_weighting` flag
7. Add initializer for `beta_propensity`
8. Test with `combined_no_trial` model variant
9. Compare posterior predictive KM curves: weighted vs unweighted

## Context from This Session

- The `karim/burden-endpoints` branch has the unified burden endpoint module,
  multistate alignment with SCLC, `combined_no_trial` model variant,
  PCWG3 constants, `calculate_all_patients_pcwg3()`, and the Stan include
  architecture refactor (`_base_hierarchy_data.stan` / `_hmc_routing_transformed_data.stan`).
- The standalone multistate model compiles and runs well (E-BFMI ~1.0).
- The full combined model struggles with HMC geometry (E-BFMI < 0.02 on 3/4 chains).
- Laplace runs are stuck at initialization — likely QR overflow from extreme RWD covariates.
- Propensity weighting could help both problems: downweighting extreme RWD patients
  improves geometry AND prevents QR overflow.
