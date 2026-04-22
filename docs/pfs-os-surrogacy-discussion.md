# PFS-OS Relationship: Surrogacy Discussion

**Date**: 2026-04-09  
**Branch**: `karim/burden-endpoints-504`  
**Context**: Analysis page `quarto/pioneer/website/analysis/pfs-os-relationship.qmd`

---

## What Started This

The original question: "Earlier PFS — does it predict earlier OS?" Led to building a PFS-OS
relationship page. Initially included:
- PFS vs OS scatter (posterior medians per patient)
- Conditional OS KM by PFS quartile (per MCMC draw)
- Spearman rank correlation between PFS and OS

## Why the Original Analysis Was Contaminated

**Key insight**: OS ≥ PFS *by construction* — if you progressed at week X, you must have
survived at least to week X. Any positive PFS-OS correlation is partly a mathematical floor
effect, not a biological signal. The Spearman analysis doesn't disentangle:

1. Mechanical correlation (OS ≥ PFS always)
2. Biological signal (earlier progression → earlier death)

## The Real Question: Is PFS a Surrogate for OS?

The goal reframed as the regulatory question: **does the treatment effect on OS flow entirely
through PFS?**

## The Illness-Death Model as a Causal DAG

The three-state model defines three causal pathways from treatment to death:

```
Treatment → h₀₁(t) → PFS → h₁₂(t) → OS
     │                              ↑
     └──→ h₀₂(t) ─────────────────→ OS (direct death)
```

| Transition | Pathway | Surrogacy condition |
|------------|---------|---------------------|
| h₀₁(t)    | 0→1 (progression) | Must have arm effect — this IS the surrogate signal |
| h₀₂(t)    | 0→2 (direct death) | Must be arm-independent |
| h₁₂(t)    | 1→2 (post-progression) | Must be arm-independent |

PFS is a **perfect surrogate** iff arm effects on h₀₂ ≈ 0 AND arm effects on h₁₂ ≈ 0.

## Critical Finding: Surrogacy Is Baked In By Construction

### How arm effects enter the model (`combined_markov_fe_arm_level_only`)

**PSA dynamics** (tr, frac, init modules):
```r
enable_level_intercept_tr   = c(trial = 0, arm = FE, patient = RE)
enable_level_intercept_frac = c(trial = 0, arm = FE, patient = RE)
enable_level_intercept_init = c(trial = 0, arm = FE, patient = RE)
```
Arms have **fixed-effect intercepts** on tumor regression, growth fraction, and initial state.
This is the only place treatment affects the PSA dynamics.

**Multistate hazards** (h₀₁, h₀₂, h₁₂):
```r
enable_ms_level_baseline_hazard = c(trial = 0, arm = 0, patient = 0)  # NO arm baseline
enable_ms_pop_time_invariant_cov = TRUE   # covariates = age, race, ecog, labs — NOT arm
enable_ms_level_cov = c(trial = 0, arm = 0, patient = 0)              # NO arm slopes
```

The covariate formula (`covar_formula` in `targets/pioneer_targets.R` lines 91–106) is:
```r
~ age + race + ecogbl + baseline_albumin + baseline_alt + baseline_chloride + 
  baseline_creatinine + baseline_neutrophils + baseline_nlr + baseline_AST + 
  baseline_ALP + baseline_leukocytes + bone_mets + lymph_mets + visceral_mets
```
**`arm` is not in it.** The `time_invariant_coef_qr_01`, `_02`, `_12` parameters act on
patient-level baseline covariates only — never on arm indicators.

### Implication

**h₀₁, h₀₂, and h₁₂ are all arm-independent in this model.** Treatment affects OS only
through PSA dynamics → PSA progression timing → PFS timing. Once you've progressed,
your post-progression survival (h₁₂) is the same regardless of which arm you were in.
Once you're at risk for direct death (h₀₂), same story.

This means the model **structurally assumes** PFS surrogacy for OS. The posterior predictive
distribution cannot falsify the surrogacy conditions because they are not parameters in the
model.

## Subtle Nuance: Timing Effect via Markov Clock

Even with arm-independent h₁₂, post-progression survival can still *appear* to differ by arm.
For `ms_time_scale_12 = markov`, h₁₂ depends on calendar time t. Arms that have different
PSA dynamics will progress at different calendar times — and h₁₂(t) evaluated at those
different progression times may give different values. This is not a treatment effect on
post-progression survival; it is the model propagating timing differences through a
non-constant h₁₂(t) baseline. Worth noting in the page to avoid misinterpretation.

## What Was Actually Built

The page `quarto/pioneer/website/analysis/pfs-os-relationship.qmd` was rewritten with:

1. **Causal DAG framing** — explains the three-state model and surrogacy conditions
2. **Arm-specific PFS and OS KM** — using pre-computed `cond_spop_pfs_km_est` /
   `cond_spop_os_km_est` (variable == "arms") from `pioneer_cond_*_km_rvar_posterior_*`
3. **Post-progression survival by arm** — computed from patient-level PPC draws
   (`all_pioneer_ppc_rvar_*`), filtering to trial patients, per-draw KM of (OS - PFS)
   among progressors, wrapped into rvar
4. **ΔPFS vs ΔOS scatter** — arm-level median PFS and OS from quantile targets,
   differences vs Arm 1 (50 kBq/kg), identity line for perfect surrogacy
5. **PFS vs OS scatter** — per-patient posterior medians faceted by arm (exit-route colored)

**Problem**: The narrative claims to "test" the surrogacy conditions, but the model
structurally assumes them. The code is correct; the framing needs revision.

## Open Questions / Next Steps

1. **Narrative reframe**: Change from "testing surrogacy" to "examining the structural
   surrogacy assumption and its implications" — what does the model predict for arm-specific
   OS under this assumption?

2. **Should arm effects be added to h₁₂ or h₀₂?** This would be a model extension:
   - Add `enable_ms_level_baseline_hazard = c(trial = 0, arm = RE, patient = 0)` for
     transitions 1→2 and/or 0→2
   - Would allow the model to learn whether post-progression survival actually differs by arm
   - Would require new priors (`ms_*_level_intercept_sd` in `r/pioneer/priors.R`)
   - This is the right way to properly test surrogacy vs assume it

3. **ΔPFS vs ΔOS is still informative** — even under structural surrogacy, the magnitudes
   of ΔOS relative to ΔPFS show how much of the OS treatment effect the model attributes
   to PFS timing differences (mediation via PSA dynamics) vs the fixed h₀₂/h₁₂ population
   baseline mixing

4. **The scatter and conditional KM** are useful for understanding the patient-level
   PFS-OS joint distribution, even if they don't test surrogacy

## Key Files

| File | Role |
|------|------|
| `targets/pioneer_targets.R` lines 88–106 | `covar_formula` — no arm |
| `targets/pioneer_targets.R` lines 504–535 | `default_settings` — multistate flags |
| `targets/pioneer_targets.R` lines 422–440 | `tar_map` model variants |
| `stan/modules/multistate/flags.stan` | `enable_ms_level_baseline_hazard` |
| `stan/modules/multistate/transformed_parameters.stan` lines 118, 247, 361 | `linpred_pop_*` computed from QR covariates |
| `r/pioneer/prepare_analysis_data.R` lines 247–295 | `prepare_flatiron_covar_matrix()` |
| `quarto/pioneer/website/analysis/pfs-os-relationship.qmd` | The page (code correct, narrative needs reframe) |
