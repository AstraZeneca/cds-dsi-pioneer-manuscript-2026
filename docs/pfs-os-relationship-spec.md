# PFS-OS Relationship Page: Design Spec (Updated)

**Date**: 2026-04-09  
**Branch**: `karim/burden-endpoints-504`  
**File**: `quarto/pioneer/website/analysis/pfs-os-relationship.qmd`  
**Registered in**: `_quarto.yml` under "Full Model > Results", after `burden-endpoints.qmd`

---

## Evolution of Goals

### Original intent
Understand whether earlier PFS predicts earlier OS. Three sections built and committed:
1. PFS vs OS scatter (posterior median per patient, colored by exit route)
2. Conditional OS KM by PFS quartile (per-draw, 5 groups: 0→2 + Q1–Q4 progressors)
3. Spearman rank correlation between PFS and OS per draw

### Why the original was unsatisfying
The Spearman analysis is **mechanically contaminated**: OS ≥ PFS by construction, so any
positive correlation is partly a mathematical floor effect (not biological signal). The
conditional KM by quartile is interesting but doesn't directly address the causal question.

### Reframed goal
**Is PFS a valid surrogate endpoint for OS?** This is the regulatory question — does the
treatment effect on OS flow entirely through the treatment effect on PFS?

---

## Illness-Death Model as a Causal DAG

The model defines three pathways from treatment to death:

```
Treatment → h₀₁(t) → PFS → h₁₂(t) → OS
     │                              ↑
     └──→ h₀₂(t) ─────────────────→ OS (direct death)
```

PFS surrogacy conditions (Prentice 1989 criteria, translated to this model):
- **h₁₂(t) arm-independent**: post-progression survival does not differ by arm
- **h₀₂(t) arm-independent**: direct death rate does not differ by arm

If both hold, treatment can only affect OS by affecting PFS timing → PFS is a valid surrogate.

---

## Critical Structural Finding

**PFS surrogacy is baked into the model by construction for `combined_markov_fe_arm_level_only`.**

### How arm effects enter the model

| Module | Arm effect | Mechanism |
|--------|-----------|-----------|
| tr (tumor regression) | **YES** — FE | `enable_level_intercept_tr[arm] = fe` |
| frac (growth fraction) | **YES** — FE | `enable_level_intercept_frac[arm] = fe` |
| init (initial state) | **YES** — FE | `enable_level_intercept_init[arm] = fe` |
| h₀₁ (0→1 MS hazard) | **NO** | `enable_ms_level_baseline_hazard[arm] = 0` |
| h₀₂ (0→2 direct death) | **NO** | same |
| h₁₂ (1→2 post-prog) | **NO** | same |

The covariate design matrix (`time_invariant_coef_qr_*` applied to all three transitions)
does NOT include arm — the formula only has baseline patient characteristics:
```r
covar_formula <- ~ age + race + ecogbl + baseline_albumin + ... + visceral_mets
```
(`targets/pioneer_targets.R` lines 91–106)

### Consequence
Treatment effects reach OS only via PSA dynamics → PSA progression timing → PFS. The
model cannot learn arm-specific h₀₂ or h₁₂ — these are pooled across arms. The posterior
predictive distribution cannot falsify surrogacy because it never had the parameters to
encode violation.

### Timing nuance (Markov clock)
`ms_time_scale_12 = markov` → h₁₂ depends on calendar time t. Arms that progress at
different calendar times (due to PSA dynamics differences) face different baseline h₁₂(t)
values. This can produce arm differences in post-progression survival in the PPC even though
h₁₂ has no arm parameter — it reflects timing differences, not treatment effects on
post-progression death.

---

## Current Page State (Code Correct, Narrative Wrong)

The page was rewritten with surrogacy framing and contains five sections:

### 1. Causal Structure (`#sec-causal-structure`)
Markdown explaining the DAG and surrogacy conditions. **Problem**: the callout says "we test
these conditions using the PPC" — but we can't, the model assumes them.

### 2. Arm-Specific PFS and OS KM (`#sec-arm-km`)
Side-by-side PFS and OS KM by arm using pre-computed targets:
- `pioneer_cond_pfs_km_rvar_posterior_{full_model}` (variable == "arms")
- `pioneer_cond_os_km_rvar_posterior_{full_model}` (variable == "arms")
Uses `cond_spop_*_km_est`, overlaid per arm with `stat_lineribbon`. **Code correct.**

### 3. Post-Progression Survival by Arm (`#sec-post-prog`)
Computed from `all_pioneer_ppc_rvar_{full_model}`, filtering to trial patients
(excluding "Pluvicto monotherapy"). Per draw: identify progressors (PFS event, OS > PFS),
compute OS − PFS, fit KM per arm, wrap in `rvar`. **Code correct.** But interpretation
callout says "overlapping curves → arm-independent h₁₂" — which is true in principle but
the curves are expected to largely overlap because the model forces this. Needs a note about
the Markov timing nuance.

### 4. Treatment Effect Mediation (`#sec-mediation`)
Arm-specific median PFS and OS from quantile targets, ΔPFS and ΔOS vs Arm 1 (50 kBq/kg),
identity scatter + surrogacy ratio table. **Code correct and still informative** — even under
structural surrogacy, this shows how arm effects propagate from PSA dynamics through to OS.

### 5. PFS vs OS Scatter (`#sec-scatter`)
Per-patient posterior medians, faceted by arm, colored by exit route. **Code correct.**
Still useful as an overview of the PFS-OS joint distribution.

---

## What Needs Changing

### Narrative (priority)
1. Section 1: Reframe from "testing surrogacy" → "examining the structural surrogacy
   assumption built into this model variant"
2. Section 3: Add note about Markov timing nuance — partial arm separation is expected from
   timing differences, not treatment effects on h₁₂
3. Title/subtitle: Update to reflect the model-assumption framing rather than "testing"

### Possible model extension (future work)
To actually test surrogacy, add arm effects to h₁₂ and/or h₀₂:
```r
enable_ms_level_baseline_hazard = c(trial = 0, arm = RE, patient = 0)
```
This would require:
- New priors for `ms_*_level_intercept_sd` for each enabled arm-level transition
- New model variant in `tar_map`
- Refit + compare DIC/ELPD with vs without arm effects on post-progression hazard

---

## Data Dependencies

| Target | Purpose |
|--------|---------|
| `all_pioneer_ppc_rvar_{full_model}` | Patient-level spop PFS/OS rvars |
| `all_analysis_data_{full_model}` | Arm/trial metadata (joined to PPC) |
| `pioneer_cond_pfs_km_rvar_posterior_{full_model}` | Arm-level PFS KM |
| `pioneer_cond_os_km_rvar_posterior_{full_model}` | Arm-level OS KM |
| `pioneer_cond_pfs_quant_posterior_{full_model}` | Arm-level median PFS |
| `pioneer_cond_quant_os_posterior_{full_model}` | Arm-level median OS |

All filtered to `variable == "arms"` where applicable.

---

## Related Docs

- `docs/pfs-os-surrogacy-discussion.md` — Full narrative of the discussion thread
- `quarto/pioneer/website/analysis/subgroup-km.qmd` — Reference for reading cond KM targets
- `quarto/pioneer/website/analysis/subgroup-burden-endpoints.qmd` — Reference for quantile targets
