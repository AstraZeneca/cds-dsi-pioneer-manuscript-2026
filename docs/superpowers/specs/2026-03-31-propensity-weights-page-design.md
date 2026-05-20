# Propensity Weights Analysis Page — Design Spec

**Date:** 2026-03-31  
**Branch:** `karim/burden-endpoints`  
**Author:** Karim Naguib

---

## Overview

Add a Quarto analysis page that visualises the posterior distribution of
propensity likelihood weights estimated by the Pioneer model. Weights are
already computed as a Stan transformed parameter (`likelihood_weight[n_patients]`)
and written to every MCMC CSV — this work exposes them through the targets
pipeline and renders them in the website.

---

## Goals

- Surface `likelihood_weight` from the MCMC output as a first-class pipeline
  target.
- Show the posterior distribution of RWD patient weights, faceted by trial/source,
  using the existing `stat_distogram` infrastructure.
- Posterior only (no prior comparison needed).

---

## Scope

In scope:
- One new `tar_target` inside the inner `type` tar_map.
- One new `all_` combiner target at the model level.
- One new Quarto page `analysis/propensity-weights.qmd`.
- Navigation registration in `_quarto.yml`.

Out of scope:
- Prior vs posterior comparison.
- Per-patient weight tables or high-weight patient identification.
- Models without propensity weighting (`no_hist`) — page loads the full model
  variant (`full_model`) which always has propensity enabled.

---

## Targets Changes

### 1. Raw draws target

Location: inside the **inner `type` tar_map** in `targets/pioneer_targets.R`,
alongside `pioneer_draws_km`, `pioneer_draws_psa`, etc.

```r
tar_target(
  pioneer_draws_propensity_weights,
  select_draws(pioneer_fit, likelihood_weight)
)
```

This selects `likelihood_weight[1]` … `likelihood_weight[n_patients]` from the
MCMC CSV files. No new Stan code required — the variable is already a transformed
parameter.

### 2. Model-level combiner

Location: the "Combined prior + posterior" block in `targets/pioneer_targets.R`,
following the same pattern as `all_pioneer_draws_km`, etc.

```r
tar_target(
  all_pioneer_draws_propensity_weights,
  bind_rows(
    pioneer_draws_propensity_weights_prior,
    pioneer_draws_propensity_weights_posterior
  ) |> mutate(model = model)
)
```

---

## Quarto Page

**File:** `quarto/pioneer/website/analysis/propensity-weights.qmd`

### Frontmatter

```yaml
---
title: "Propensity Weights"
subtitle: "Posterior distribution of per-patient likelihood weights"
---
```

### Setup chunk

```r
source("_setup.R")

# Always use "combined" — it's the only model with propensity enabled.
# full_model defaults to "no_hist" which has all weights = 1.0.
propensity_model <- "combined"

weights_raw <- tar_read_raw(
  str_glue("all_pioneer_draws_propensity_weights_{propensity_model}"),
  store = pioneer_store
)

analysis_data <- tar_read_raw(
  str_glue("all_analysis_data_{propensity_model}"),
  store = pioneer_store
) |>
  select(-any_of("visit_data")) |>
  mutate(source = if_else(trial == "FPI-2265-202", "Trial", "RWD"))

# Posterior only; spread to per-patient rvar; exclude trial patients (always weight = 1)
weights_rvar <- weights_raw |>
  filter(fit_type == "posterior") |>
  spread_rvars(likelihood_weight[i]) |>
  bind_cols(analysis_data) |>
  filter(source == "RWD")
```

### Plot

A distogram of `likelihood_weight` using `stat_distogram`, faceted by `trial`,
x-axis = weight (0–1), y-axis = patient count. Uses the project AZ colour palette.

The distogram shows uncertainty in the weights across MCMC draws as ribbons
(50%, 80%, 95% CIs), with the median as a line — identical to how patient-level
rates are shown in the SCLC-01 parameter inference page.

### Narrative

Short callout note explaining: trial patients are excluded (always weight = 1);
the distribution shows how "trial-like" each Flatiron patient's covariates are;
right-skewed distribution is expected; effective sample size
(`sum(posterior_mean_weights)`) is reported as a scalar.

---

## Navigation (`_quarto.yml`)

### `render` list

Add after `analysis/parameter-inference.qmd`:
```yaml
- "analysis/propensity-weights.qmd"
```

### Sidebar (Full Model → Model Parameters section)

Add after `analysis/parameter-inference.qmd`:
```yaml
- analysis/propensity-weights.qmd
```

---

## File Checklist

| File | Change |
|------|--------|
| `targets/pioneer_targets.R` | Add `pioneer_draws_propensity_weights` inside inner tar_map |
| `targets/pioneer_targets.R` | Add `all_pioneer_draws_propensity_weights` combiner |
| `quarto/pioneer/website/analysis/propensity-weights.qmd` | New file |
| `quarto/pioneer/website/_quarto.yml` | Add to render list + sidebar |

---

## Non-Goals / Deferred

- Per-covariate coefficient plots (`beta_propensity`) — deferred; can be added
  to `parameter-inference.qmd` later.
- `no_hist` / `combined_ungated` model comparisons.
- High-weight patient identification table.
