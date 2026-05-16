# Propensity Weights Analysis Page — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a `likelihood_weight` draws target to the pipeline and a Quarto analysis page that displays a distogram of posterior propensity weights, faceted by Flatiron cohort (trial), for the `combined` model variant.

**Architecture:** (1) A new `tar_target` inside the inner `type` tar_map extracts `likelihood_weight[i]` from the MCMC CSVs using the existing `select_draws()` utility. (2) A model-level combiner merges prior + posterior draws following the established `all_pioneer_*` pattern. (3) A new Quarto page reads the combined draws, spreads them into per-patient rvars via `spread_rvars`, joins patient metadata from the store, and plots a `stat_distogram` faceted by trial — all using existing infrastructure.

**Tech Stack:** R, targets, CmdStanR, tidybayes (`spread_rvars`), ggplot2, ggdist (`stat_distogram`, `geom_lineribbon`), posterior (rvar)

---

## File Map

| File | Action | What changes |
|------|--------|--------------|
| `targets/pioneer_targets.R` | Modify | Add `pioneer_draws_propensity_weights` inside inner tar_map (~line 546 block) |
| `targets/pioneer_targets.R` | Modify | Add `all_pioneer_draws_propensity_weights` combiner in the "Combined prior + posterior" block (~line 989) |
| `quarto/pioneer/website/analysis/propensity-weights.qmd` | Create | New Quarto analysis page |
| `quarto/pioneer/website/_quarto.yml` | Modify | Add page to `render` list and Full Model sidebar |

---

## Task 1: Add raw draws target inside inner tar_map

**Files:**
- Modify: `targets/pioneer_targets.R`

The inner `tar_map` iterates over `type = c("prior", "posterior")` and produces
targets like `pioneer_draws_km_posterior_combined`. Add a new draws target for
`likelihood_weight` in the same block.

- [ ] **Step 1: Locate the insertion point**

In `targets/pioneer_targets.R`, find the `# PSA Draws Extraction` comment
block (around line 544–598). The new target goes after `pioneer_draws_ppc`
(the last draws target before the rvar-processing targets).

- [ ] **Step 2: Add the target**

Insert after the `pioneer_draws_ppc` tar_target block (before the
`pioneer_ppc_rvar` target at ~line 683):

```r
      # Propensity likelihood weights (per-patient, capped density ratios)
      tar_target(
        pioneer_draws_propensity_weights,
        select_draws(pioneer_fit, likelihood_weight)
      ),
```

- [ ] **Step 3: Verify syntax**

```bash
Rscript -e 'source("targets/pioneer_targets.R")' 2>&1 | head -20
```

Expected: no parse errors. (The script will emit warnings about missing env vars
like `TAR_PROJECT` — that's fine.)

- [ ] **Step 4: Commit**

```bash
git add targets/pioneer_targets.R
git commit -m "feat(pioneer): add pioneer_draws_propensity_weights target"
```

---

## Task 2: Add model-level combiner target

**Files:**
- Modify: `targets/pioneer_targets.R`

The "Combined prior + posterior" block (~line 989) assembles `all_pioneer_*`
targets by binding the prior and posterior variants. Add the propensity weights
combiner to this block.

- [ ] **Step 1: Locate the insertion point**

Find the line:
```r
    # Combined prior + posterior (full model) ================================
```

The block contains ~30 `tar_target(all_pioneer_*)` entries. Add the new
combiner after `all_pioneer_ppc_rvar` (the last "draws" combiner before the
quantile combiners, ~line 1024).

- [ ] **Step 2: Add the combiner**

```r
    tar_target(
      all_pioneer_draws_propensity_weights,
      bind_rows(
        pioneer_draws_propensity_weights_prior,
        pioneer_draws_propensity_weights_posterior
      ) |> mutate(model = model)
    ),
```

- [ ] **Step 3: Verify syntax**

```bash
Rscript -e 'source("targets/pioneer_targets.R")' 2>&1 | head -20
```

Expected: no parse errors.

- [ ] **Step 4: Commit**

```bash
git add targets/pioneer_targets.R
git commit -m "feat(pioneer): add all_pioneer_draws_propensity_weights combiner"
```

---

## Task 3: Create the Quarto analysis page

**Files:**
- Create: `quarto/pioneer/website/analysis/propensity-weights.qmd`

The page loads posterior-only `likelihood_weight` draws for the `combined` model,
joins patient metadata to identify RWD patients, then plots a `stat_distogram`
faceted by `trial`.

Note on `stat_distogram`: it extends `ggdist:::StatLineribbon`. It expects a
column mapped to `ydist` that is an rvar (from `spread_rvars`). The stat
internally calls `rvar_sample_hist(ydist, breaks, freq = TRUE)` to bin the
distributions across patients into a per-bin count rvar, then `geom_lineribbon`
draws median + CI ribbons for each bin as a step histogram. x-axis = weight
value (0–1), y-axis = patient count with credible interval ribbons.

- [ ] **Step 1: Create the file**

```r
---
title: "Propensity Weights"
subtitle: "Posterior distribution of per-patient likelihood weights (combined model)"
date-modified: "2026-03-31"
number-sections: true
---

```{r setup}
#| include: false

source("_setup.R")

# Always use "combined" — it's the only variant with propensity enabled.
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

# Posterior only; spread to per-patient rvar; exclude trial patients (always 1.0)
weights_rvar <- weights_raw |>
  filter(fit_type == "posterior") |>
  spread_rvars(likelihood_weight[i]) |>
  bind_cols(analysis_data) |>
  filter(source == "RWD")
```

## Weight Distribution {#sec-weights}

Each RWD (Flatiron) patient receives a posterior likelihood weight
$w_i = \min(1,\, P(\mathbf{x}_i \mid \text{trial}) / P(\mathbf{x}_i \mid \text{RWD}))$
that reflects how "trial-like" their baseline covariates are. Trial patients are
always clamped to 1.0 and excluded from these plots.

::: {.callout-note}
## Interpreting the distogram

- **x-axis**: likelihood weight (0 = no contribution; 1 = full contribution)
- **y-axis**: number of RWD patients in each weight bin
- **Ribbons**: 50%, 80%, 95% posterior credible intervals on the count per bin
- A right-skewed distribution is expected — most RWD patients have low overlap
  with the trial population; a small subset is highly trial-like
:::

```{r}
#| label: fig-weight-distogram
#| fig-cap: !expr str_glue("Posterior distribution of propensity likelihood weights for RWD patients — {propensity_model} model. Faceted by Flatiron cohort.")
#| fig-width: 9
#| fig-height: 5

ggplot(weights_rvar, aes(ydist = likelihood_weight)) +
  stat_distogram(
    breaks = seq(0, 1, by = 0.05),
    .width = c(0.5, 0.8, 0.95),
    fill = AZ_navy,
    color = AZ_navy,
    alpha = 0.3
  ) +
  facet_wrap(vars(trial), scales = "free_y") +
  scale_x_continuous("Likelihood weight", limits = c(0, 1),
                     breaks = seq(0, 1, by = 0.2)) +
  scale_y_continuous("Patient count") +
  labs(
    caption = "Ribbons show 50%, 80%, 95% posterior credible intervals on patient count per bin."
  ) +
  theme(plot.caption = element_text(hjust = 0, size = rel(0.9), margin = margin(t = 6)))
```

## Effective Sample Size {#sec-ess}

The effective number of RWD patients contributing to the posterior is
$\text{ESS} = \sum_{i \in \text{RWD}} w_i$. An ESS much smaller than the
nominal RWD sample size indicates that only a few patients have meaningful
covariate overlap with the trial.

```{r}
#| label: tbl-ess
#| tbl-cap: !expr str_glue("Effective sample size (ESS = sum of weights) by Flatiron cohort — {propensity_model} model, posterior.")

weights_rvar |>
  summarise(
    n_rwd = n(),
    ess = median(E(likelihood_weight)),
    ess_pct = round(100 * ess / n_rwd, 1),
    .by = trial
  ) |>
  rename(
    Trial = trial,
    `N (RWD)` = n_rwd,
    `ESS (posterior mean)` = ess,
    `ESS %` = ess_pct
  ) |>
  knitr::kable(digits = 1)
```
```

(Wrap the R code chunks in proper backtick fences in the file — the plan uses
plain indentation here to avoid confusing the markdown renderer.)

- [ ] **Step 2: Verify the file parses**

```bash
Rscript -e 'knitr::parse_block(readLines("quarto/pioneer/website/analysis/propensity-weights.qmd"))' 2>&1 | head -10
```

Expected: no errors (only potential warnings about undefined vars at parse time).

- [ ] **Step 3: Commit**

```bash
git add quarto/pioneer/website/analysis/propensity-weights.qmd
git commit -m "feat(pioneer): add propensity-weights Quarto analysis page"
```

---

## Task 4: Register the page in `_quarto.yml`

**Files:**
- Modify: `quarto/pioneer/website/_quarto.yml`

Two insertions: the `render` list and the Full Model sidebar.

- [ ] **Step 1: Add to the render list**

In `_quarto.yml`, find:
```yaml
    - "analysis/parameter-inference.qmd"
```
Add immediately after:
```yaml
    - "analysis/propensity-weights.qmd"
```

- [ ] **Step 2: Add to the sidebar**

Find the Full Model → Model Parameters sidebar section:
```yaml
        - section: "Model Parameters"
          contents:
            - analysis/parameter-inference.qmd
```
Add immediately after `parameter-inference.qmd`:
```yaml
            - analysis/propensity-weights.qmd
```

- [ ] **Step 3: Verify YAML syntax**

```bash
Rscript -e 'yaml::read_yaml("quarto/pioneer/website/_quarto.yml")' 2>&1 | head -5
```

Expected: no errors.

- [ ] **Step 4: Commit**

```bash
git add quarto/pioneer/website/_quarto.yml
git commit -m "feat(pioneer): register propensity-weights page in navigation"
```

---

## Self-Review Checklist

After completing all tasks, verify:

- [ ] `Rscript -e 'source("targets/pioneer_targets.R")'` parses without errors
- [ ] `yaml::read_yaml("quarto/pioneer/website/_quarto.yml")` parses without errors
- [ ] `propensity-weights.qmd` is listed in both `render:` and the Full Model sidebar in `_quarto.yml`
- [ ] `pioneer_draws_propensity_weights` appears in the inner `type` tar_map (before the rvar-processing targets)
- [ ] `all_pioneer_draws_propensity_weights` appears in the "Combined prior + posterior" block
- [ ] Page uses `propensity_model <- "combined"` (not `full_model`)
- [ ] Page filters to `source == "RWD"` before plotting (trial patients excluded)
