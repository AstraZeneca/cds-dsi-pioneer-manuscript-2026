# Publication Website Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Create a Quarto analysis website at `quarto/publication/website/` for the publication pipeline, mirroring the sclc site with minimal duplication.

**Architecture:** Follows the established pioneer pattern — own `_quarto.yml` + `_setup.R` sourcing `quarto/_shared/_setup.R`, own `az-theme.scss` importing `_shared/az-colors.scss`. No documentation section (sclc is canonical). Publication pages call low-level plot primitives directly (no PDL1/DCO filtering needed since publication has one dataset with no stratification). Missing pipeline targets are added to `targets/publication_targets.R` before building the site.

**Tech Stack:** R, Quarto, targets, tidybayes (`spread_rvars`, `gather_rvars`), ggdist, `r/plot_functions.R`, `r/sclc/plot_functions.R`

---

## File Map

**Create:**
- `r/publication/target_stores.R` — defines `pub_store` path
- `quarto/publication/website/_quarto.yml` — site config, navbar
- `quarto/publication/website/_metadata.yml` — date-modified metadata
- `quarto/publication/website/_setup.R` — sources shared setup + pub store + sclc plot fns
- `quarto/publication/website/az-theme.scss` — imports `_shared/az-colors`, no changes
- `quarto/publication/website/styles.css` — identical to sclc (hero/card/figure CSS)
- `quarto/publication/website/index.qmd` — landing page
- `quarto/publication/website/analysis/_metadata.yml`
- `quarto/publication/website/analysis/data-exploration.qmd`
- `quarto/publication/website/analysis/survival-analysis.qmd`
- `quarto/publication/website/analysis/overall-survival.qmd`
- `quarto/publication/website/analysis/clinical-endpoints.qmd`
- `quarto/publication/website/analysis/parameter-inference.qmd`
- `quarto/publication/website/analysis/patient-dynamics.qmd`
- `quarto/publication/website/analysis/diagnostics.qmd`

**Modify:**
- `targets/publication_targets.R` — add ORR, pfs-n, patient-level, noise, bpi targets + combined aggregates

---

## Task 1: Add missing pipeline targets — ORR and PFS-n

**Files:**
- Modify: `targets/publication_targets.R`

These two targets live inside the existing `tar_map(tibble(type = ...), ...)` block, alongside the existing `tumor_ssls_km_rvar` and `tumor_ssls_trial_pfs_quant` targets. Add them just before the `# Population parameters` comment.

- [ ] **Step 1: Add ORR and pfs-n targets inside `tar_map`**

In `targets/publication_targets.R`, find the block ending with:
```r
    ),

    # Population parameters -----------------------------------------------------
```

Insert before that comment:
```r
    tar_target(
      tumor_ssls_orr_rvar,
      tumor_ssls_draws_endpoints |>
        recover_types(select(all_analysis_data, trial)) |>
        spread_rvars(
          sample_target_orr[trial],
          spop_target_orr[trial]
        ) |>
        mutate(fit_type = type)
    ),

    tar_target(
      tumor_ssls_forecast_target_pfs_n_rvar,
      tumor_ssls_draws_endpoints |>
        recover_types(select(all_analysis_data, trial)) |>
        spread_rvars(
          sample_target_pfs_n[trial, n],
          spop_target_pfs_n[trial, n],
          sample_ms_pfs_n[trial, n],
          spop_ms_pfs_n[trial, n],
          sample_pfs_n[trial, n],
          spop_pfs_n[trial, n]
        ) |>
        left_join(pfs_timepoints_pub, by = "n") |>
        mutate(fit_type = type)
    ),
```

- [ ] **Step 2: Add combined aggregates**

In the `# Combined prior + posterior` section at the bottom of the file, add after the existing `all_tumor_ssls_trial_pfs_quant` target:
```r
  tar_target(
    all_tumor_ssls_orr_rvar,
    bind_rows(tumor_ssls_orr_rvar_prior, tumor_ssls_orr_rvar_posterior)
  ),
  tar_target(
    all_tumor_ssls_forecast_target_pfs_n_rvar,
    bind_rows(
      tumor_ssls_forecast_target_pfs_n_rvar_prior,
      tumor_ssls_forecast_target_pfs_n_rvar_posterior
    )
  ),
```

- [ ] **Step 3: Validate R syntax**

```bash
Rscript --no-init-file -e 'source("targets/publication_targets.R")' 2>&1 | head -20
```
Expected: no errors (will warn about missing objects like `tumor_ssls_draws_endpoints` — that's fine, they only exist at runtime).

- [ ] **Step 4: Commit**

```bash
git add targets/publication_targets.R
git commit -m "feat(publication): add ORR and pfs-n endpoint targets"
```

---

## Task 2: Add missing pipeline targets — noise SD and patient binned params

**Files:**
- Modify: `targets/publication_targets.R`

These are needed for `parameter-inference.qmd`. They live inside the `tar_map` alongside the existing `tumor_ssls_rates_rvar`.

- [ ] **Step 1: Add noise_sd, rates_bpi, decrease_prop_bpi inside `tar_map`**

Add after the existing `tumor_ssls_rates_rvar` target (just before `tumor_ssls_coef`):

```r
    tar_target(
      tumor_ssls_noise_sd_rvar,
      gather_rvars(tumor_ssls_draws_pop, measure_sd_sld) |>
        mutate(fit_type = type)
    ),

    tar_map(
      tibble(level = c("patient")),
      names = "level",

      tar_target(
        tumor_ssls_rates_bpi,
        get_tumor_ssls_level_param_binned(
          tumor_ssls_draws_patient_params,
          level,
          param = str_c(
            "{level}_",
            c("log_decrease_rate", "log_growth_rate",
              "log_growth_rate_residual", "log_decrease_rate_residual")
          ),
          type,
          breaks = seq(-3, 3, 0.05),
          inv_link_breaks = seq(0, 25, 0.5)
        )
      ),

      tar_target(
        tumor_ssls_decrease_prop_bpi,
        get_tumor_ssls_level_param_binned(
          tumor_ssls_draws_patient_params,
          level,
          param = str_c("frac_logit_loc_{level}"),
          type,
          breaks = seq(-5, 5, 0.05),
          inv_link = rvar_plogis,
          inv_link_breaks = seq(0, 1, 0.01)
        )
      )
    ),
```

- [ ] **Step 2: Add combined aggregates**

In the combined section at the bottom, add after `all_tumor_ssls_coef`:

```r
  tar_target(
    all_tumor_ssls_noise_sd_rvar,
    bind_rows(tumor_ssls_noise_sd_rvar_prior, tumor_ssls_noise_sd_rvar_posterior)
  ),
  tar_target(
    all_tumor_ssls_patient_rates_bpi,
    bind_rows(
      tumor_ssls_rates_bpi_patient_prior,
      tumor_ssls_rates_bpi_patient_posterior
    )
  ),
  tar_target(
    all_tumor_ssls_patient_decrease_prop_bpi,
    bind_rows(
      tumor_ssls_decrease_prop_bpi_patient_prior,
      tumor_ssls_decrease_prop_bpi_patient_posterior
    )
  ),
```

- [ ] **Step 3: Validate R syntax**

```bash
Rscript --no-init-file -e 'source("targets/publication_targets.R")' 2>&1 | head -20
```

- [ ] **Step 4: Commit**

```bash
git add targets/publication_targets.R
git commit -m "feat(publication): add noise SD and patient binned parameter targets"
```

---

## Task 3: Add missing pipeline targets — patient dynamics (SLD + RECIST)

**Files:**
- Modify: `targets/publication_targets.R`

These require a nested `tar_map` over `event_type` (right_censored / uncensored) *inside* the outer `tar_map` over `type`. This is how sclc does it. The outer `tar_map` already ends with `tumor_ssls_coef`. Add the nested map after `tumor_ssls_decrease_prop_bpi` (before closing the outer `tar_map`).

- [ ] **Step 1: Add nested tar_map for patient dynamics inside outer tar_map**

In `targets/publication_targets.R`, find the inner `tar_map(tibble(level = c("patient")), ...)` block you just added and after its closing `)`, still inside the outer `tar_map`, add:

```r
    tar_map(
      tibble(
        event_type   = c("right_censored", "uncensored"),
        event_cond   = c(expr(right_censored), expr(!right_censored)),
        event_slicer = c(\(d, n) d, \(d, n) slice_sample(d, n = n))
      ),
      names = "event_type",

      tar_target(
        state_patient_subsample,
        get_state_patients(
          all_analysis_data,
          by = trial,
          cond = event_cond,
          slicer = event_slicer,
          sample_size = 40
        )
      ),

      tar_target(
        tumor_ssls_rep_sld_rvar,
        get_sld(tumor_ssls_draws_sld_recist, state_patient_subsample) |>
          mutate(fit_type = type),
        resources = tar_resources(crew = tar_resources_crew(controller = "many samples"))
      ),

      tar_target(
        tumor_ssls_forecast_sld_rvar,
        get_forecast_sld(
          tumor_ssls_draws_sld_recist,
          all_analysis_data,
          state_patient_subsample,
          forecast_extent = extend_max_all_t
        ) |>
          mutate(fit_type = type),
        resources = tar_resources(crew = tar_resources_crew(controller = "many samples"))
      ),

      tar_target(
        tumor_ssls_staged_sld_rvar,
        bind_rows(
          obs = tumor_ssls_rep_sld_rvar |> rename(patient_sld = rep_patient_sld),
          forecast = tumor_ssls_forecast_sld_rvar |> rename(patient_sld = forecast_patient_sld),
          .id = "stage"
        )
      ),

      tar_target(
        tumor_ssls_recist_rvar,
        get_recist(tumor_ssls_draws_sld_recist, state_patient_subsample) |>
          mutate(fit_type = type),
        resources = tar_resources(crew = tar_resources_crew(controller = "many samples"))
      ),

      tar_target(
        tumor_ssls_forecast_recist_rvar,
        get_forecast_recist(
          tumor_ssls_draws_sld_recist,
          all_analysis_data,
          state_patient_subsample,
          forecast_extent = extend_max_all_t
        ) |>
          mutate(fit_type = type),
        resources = tar_resources(crew = tar_resources_crew(controller = "many samples"))
      ),

      tar_target(
        tumor_ssls_staged_recist_rvar,
        bind_rows(
          obs = tumor_ssls_recist_rvar |>
            filter(!is.na(response)) |>
            rename(recist = rep_recist),
          forecast = tumor_ssls_forecast_recist_rvar |>
            rename(recist = forecast_obs_recist),
          .id = "stage"
        )
      )
    ),
```

No combined aggregates for the SLD/RECIST targets — the patient-dynamics page reads the individual per-event-type targets directly (see Task 12), matching how sclc does it.

- [ ] **Step 3: Validate R syntax**

```bash
Rscript --no-init-file -e 'source("targets/publication_targets.R")' 2>&1 | head -20
```

- [ ] **Step 4: Validate targets graph structure (requires a running R session with TAR_PROJECT=publication)**

```bash
TAR_PROJECT=publication Rscript --no-init-file -e '
  source(".Rprofile"); init_project()
  source("targets/publication_targets.R")
  names(publication_targets) |> head(20)
'
```
Expected: list of target names including `tumor_ssls_staged_sld_rvar_right_censored_posterior`, etc.

- [ ] **Step 5: Commit**

```bash
git add targets/publication_targets.R
git commit -m "feat(publication): add patient dynamics (SLD, RECIST) pipeline targets"
```

---

## Task 4: Create R infrastructure — target stores

**Files:**
- Create: `r/publication/target_stores.R`

- [ ] **Step 1: Create the file**

```r
# Target Store Path Definitions for Publication Project
#
# Prerequisites:
# - output_path must be defined (set by init_project())
# - TAR_RUN environment variable should be set

pub_store <- file.path(output_path, "publication", Sys.getenv("TAR_RUN", "main"), "_targets")
```

- [ ] **Step 2: Commit**

```bash
git add r/publication/target_stores.R
git commit -m "feat(publication): add target_stores.R"
```

---

## Task 5: Scaffold the website — config and theme files

**Files:**
- Create: `quarto/publication/website/_quarto.yml`
- Create: `quarto/publication/website/_metadata.yml`
- Create: `quarto/publication/website/az-theme.scss`
- Create: `quarto/publication/website/styles.css`
- Create: `quarto/publication/website/analysis/_metadata.yml`

- [ ] **Step 1: Create `_quarto.yml`**

```yaml
project:
  type: website
  output-dir: _site
  execute-dir: project
  render:
    - "*.qmd"
    - "**/*.qmd"

website:
  title: "PIONEER: Publication"
  description: "Bayesian hierarchical modeling for tumor dynamics and survival — publication analysis"
  favicon: ../../_shared/images/pioneer-helmet.png
  navbar:
    background: primary
    left:
      - icon: house-fill
        href: index.qmd
        aria-label: "Home"
      - text: Data
        href: analysis/data-exploration.qmd
      - text: Analysis
        menu:
          - text: Results
          - text: Clinical Endpoints
            href: analysis/clinical-endpoints.qmd
          - text: Survival Analysis
            href: analysis/survival-analysis.qmd
          - text: Overall Survival
            href: analysis/overall-survival.qmd
          - text: "---"
          - text: Model Details
          - text: Parameter Inference
            href: analysis/parameter-inference.qmd
          - text: Patient Dynamics
            href: analysis/patient-dynamics.qmd
          - text: "---"
          - text: Validation
          - text: Diagnostics
            href: analysis/diagnostics.qmd
    right:
      - icon: book
        href: https://rstudio-connect.seml.scp.astrazeneca.net/content/084c017e-8ec4-4efd-ae02-e4a86c9dc51f/
        aria-label: "SCLC-01 Documentation"
      - icon: github
        href: https://github.com/azu-oncology-rd/cds-dsi-pioneer-sclc-01-2025
  sidebar:
    - title: "Data"
      style: docked
      contents:
        - analysis/data-exploration.qmd
    - title: "Analysis"
      style: docked
      contents:
        - section: "Results"
          contents:
            - analysis/clinical-endpoints.qmd
            - analysis/survival-analysis.qmd
            - analysis/overall-survival.qmd
        - section: "Model Details"
          contents:
            - analysis/parameter-inference.qmd
            - analysis/patient-dynamics.qmd
        - section: "Validation"
          contents:
            - analysis/diagnostics.qmd
  page-navigation: true
  search:
    location: navbar
    type: overlay
  page-footer:
    left: "Last updated: {{< meta date-modified >}}"
    center: "PIONEER: Publication"
    right: "Analysis v1.0.0"

metadata-files:
  - _metadata.yml

format:
  html:
    theme:
      light: [cosmo, az-theme.scss]
      dark: [darkly, az-theme.scss]
    css: styles.css
    toc: true
    toc-depth: 3
    toc-expand: 2
    code-fold: true
    code-tools: true
    code-copy: true
    lightbox: true
    number-sections: true
    fig-cap-location: bottom
    fig-align: center
    fig-responsive: true
    fig-width: 6
    fig-height: 4.5
    out-width: "100%"
    reference-location: margin
    date-format: "MMMM D, YYYY"

bibliography: ../../../bib/pioneer.bib

execute:
  freeze: auto
  warning: false
```

- [ ] **Step 2: Create `_metadata.yml`**

```yaml
date-modified: last-modified
```

- [ ] **Step 3: Create `az-theme.scss`**

```scss
// Import shared AZ color palette
@import "../../_shared/az-colors";

/*-- scss:defaults --*/

$font-family-sans-serif: Arial, "Helvetica Neue", Helvetica, sans-serif;

$primary: $az-navy;
$secondary: $az-platinum;
$success: $az-green;
$info: $az-turquoise;
$warning: $az-gold;
$danger: $az-pink;
$light: #f8f9fa;
$dark: $az-darkgrey;

$link-color: $az-navy;
$link-hover-color: $az-plum;

$navbar-bg: $az-navy;
$navbar-fg: white;

/*-- scss:rules --*/

:root {
  --az-plum: #{$az-plum};
  --az-gold: #{$az-gold};
  --az-turquoise: #{$az-turquoise};
  --az-darkpurple: #{$az-darkpurple};
  --az-green: #{$az-green};
  --az-navy: #{$az-navy};
  --az-platinum: #{$az-platinum};
  --az-darkgrey: #{$az-darkgrey};
  --az-pink: #{$az-pink};
  --az-lightpurple: #{$az-lightpurple};
}

.navbar {
  background-color: $az-navy !important;
  border-bottom: 2px solid $az-gold;
}

.navbar-brand, .navbar-nav .nav-link {
  color: white !important;
}

.navbar-nav .nav-link:hover {
  color: $az-turquoise !important;
}

.sidebar-navigation .sidebar-item .sidebar-item-text:hover {
  color: $az-navy;
}

.callout-note      { border-left-color: $az-navy; }
.callout-tip       { border-left-color: $az-turquoise; }
.callout-warning   { border-left-color: $az-gold; }
.callout-important { border-left-color: $az-pink; }

.feature-card:hover { border-color: $az-navy; }

.hero-section { border-bottom: 3px solid $az-gold; }
.hero-image img { border: 3px solid $az-gold; }

hr { border-color: $az-gold; opacity: 0.5; }

.table thead th { border-bottom: 2px solid $az-gold; }
```

- [ ] **Step 4: Copy `styles.css` from sclc**

```bash
cp quarto/sclc/website/styles.css quarto/publication/website/styles.css
```

- [ ] **Step 5: Create `analysis/_metadata.yml`**

```yaml
execute:
  echo: true
  warning: false
  message: false
```

- [ ] **Step 6: Commit**

```bash
git add quarto/publication/website/
git commit -m "feat(publication): scaffold website config, theme, and CSS"
```

---

## Task 6: Create `_setup.R`

**Files:**
- Create: `quarto/publication/website/_setup.R`

- [ ] **Step 1: Create the file**

```r
here::i_am("quarto/publication/website/_setup.R")

source(here::here("quarto", "_shared", "_setup.R"))

source(here::here("r", "sclc", "plot_functions.R"))
source(here::here("r", "sclc", "table_functions.R"))
source(here::here("r", "sclc", "accuracy.R"))

source(here::here("r", "publication", "target_stores.R"))

targets::tar_config_set(store = pub_store)

ggplot2::theme_set(
  theme_minimal(base_family = "Arial") +
  ggplot2::theme(
    plot.margin = ggplot2::margin(5, 10, 20, 20, "pt"),
    plot.caption = ggplot2::element_text(hjust = 0, margin = ggplot2::margin(15, 0, 0, 0, "pt"))
  )
)
```

- [ ] **Step 2: Commit**

```bash
git add quarto/publication/website/_setup.R
git commit -m "feat(publication): add website _setup.R"
```

---

## Task 7: Create `index.qmd`

**Files:**
- Create: `quarto/publication/website/index.qmd`

- [ ] **Step 1: Create the file**

```qmd
---
title: "PIONEER: Publication"
subtitle: "Bayesian Hierarchical Modeling for Tumor Dynamics and Survival"
page-layout: full
toc: false
number-sections: false
---

```{r}
#| label: setup
#| include: false

source("_setup.R")
```

<style>
.quarto-title-block { display: none; }
</style>

::: {.column-page}
::: {.hero-section}
::: {.hero-content}
::: {.hero-text}
# PIONEER Publication Analysis

Bayesian hierarchical model for tumor dynamics and progression-free survival — publication dataset (497 patients, 6 covariates).
:::
::: {.hero-image}
![](../../_shared/images/pioneer-helmet.png)
:::
:::
:::
:::

::: {.column-page}
## Quick Navigation

::: {.card-grid}

::: {.feature-card}
### [Data Exploration](analysis/data-exploration.qmd)
Overview of the publication dataset: patient characteristics, tumor measurements, and survival follow-up.
:::

::: {.feature-card}
### [Clinical Endpoints](analysis/clinical-endpoints.qmd)
Overall Response Rate (ORR), Median PFS, and PFS at fixed timepoints from the publication model.
:::

::: {.feature-card}
### [Survival Analysis](analysis/survival-analysis.qmd)
Kaplan-Meier survival curves comparing model forecasts to observed PFS.
:::

::: {.feature-card}
### [Overall Survival](analysis/overall-survival.qmd)
OS model forecasts and observed Kaplan-Meier curves.
:::

::: {.feature-card}
### [Parameter Inference](analysis/parameter-inference.qmd)
Population-level and patient-level posterior parameter distributions.
:::

::: {.feature-card}
### [Patient Dynamics](analysis/patient-dynamics.qmd)
Individual patient tumor trajectories and RECIST predictions.
:::

:::
:::

::: {.column-page}
## About This Analysis

This website presents results from the publication pipeline running the SF-SSM Stan model on 497 patients with 6 covariates (age, sex, ECOG, haemoglobin, log-LDH, albumin). No PDL1 stratification or external historical trial borrowing. For model documentation, see the [SCLC-01 site](https://rstudio-connect.seml.scp.astrazeneca.net/content/084c017e-8ec4-4efd-ae02-e4a86c9dc51f/).
:::
```

- [ ] **Step 2: Commit**

```bash
git add quarto/publication/website/index.qmd
git commit -m "feat(publication): add landing page"
```

---

## Task 8: Create analysis pages — Data Exploration and Diagnostics

**Files:**
- Create: `quarto/publication/website/analysis/data-exploration.qmd`
- Create: `quarto/publication/website/analysis/diagnostics.qmd`

These two pages are the simplest — data-exploration reads raw `all_analysis_data`, diagnostics reads the fit object directly.

- [ ] **Step 1: Create `data-exploration.qmd`**

```qmd
---
title: "Data Exploration"
subtitle: "Overview of the publication analysis dataset"
---

```{r}
#| label: setup
#| include: false

source("_setup.R")
all_analysis_data <- targets::tar_read(all_analysis_data, store = pub_store)
```

This page provides an overview of the publication dataset.

## Patient Characteristics {#sec-patients}

```{r}
#| label: tbl-patient-characteristics
#| tbl-cap: Summary of patient baseline characteristics in the publication dataset.

all_analysis_data |>
  distinct(usubjid, age, male, ecog, hgb, ldh_log, albumin) |>
  summary()
```

## Enrollment Timeline {#sec-enrollment}

```{r}
#| label: fig-enrollment-timeline
#| fig-cap: Patient treatment start dates in the publication dataset.
#| fig-width: 10
#| fig-height: 4

all_analysis_data |>
  distinct(usubjid, trtsdt) |>
  ggplot(aes(x = trtsdt)) +
  geom_histogram(binwidth = 14, fill = "#003865") +
  scale_x_date("Treatment Start Date") +
  scale_y_continuous("Patients enrolled") +
  labs(title = "Publication cohort enrollment") +
  NULL
```

## SLD Measurements {#sec-sld}

```{r}
#| label: fig-sld-observations
#| fig-cap: Observed SLD trajectories for all patients.
#| fig-width: 10
#| fig-height: 6

all_analysis_data |>
  filter(!is.na(sld)) |>
  ggplot(aes(x = t, y = sld, group = usubjid)) +
  geom_line(alpha = 0.2, linewidth = 0.3, color = "#003865") +
  scale_x_continuous("Weeks") +
  scale_y_continuous("SLD [cm]") +
  NULL
```

## Survival Follow-up {#sec-survival}

```{r}
#| label: fig-survival-followup
#| fig-cap: PFS and follow-up time distribution.
#| fig-width: 8
#| fig-height: 4

all_analysis_data |>
  distinct(usubjid, pfs, right_censored) |>
  ggplot(aes(x = pfs, fill = !right_censored)) +
  geom_histogram(binwidth = 4, position = "stack") +
  scale_x_continuous("PFS [Weeks]", breaks = seq(0, 200, 26)) +
  scale_y_continuous("Patients") +
  scale_fill_manual("Event", values = c("TRUE" = "#F0AB00", "FALSE" = "#003865"),
                    labels = c("TRUE" = "Progression/Death", "FALSE" = "Censored")) +
  NULL
```
```

- [ ] **Step 2: Create `diagnostics.qmd`**

```qmd
---
title: "Sampling Diagnostics"
subtitle: "MCMC convergence and sampling quality assessment"
---

```{r}
#| label: setup
#| include: false

source("_setup.R")
```

This page presents MCMC sampling diagnostics for the publication model fit.

::: {.callout-note}
## Interpretation
- **Divergences**: Should be 0 or very low
- **E-BFMI**: Should be > 0.3
- **R-hat**: Should be < 1.01 for all parameters
- **ESS**: Effective sample size — should be sufficiently large
:::

## Diagnostics Summary {#sec-summary}

```{r}
#| label: diagnostics-summary

targets::tar_read(tumor_ssls_nuts_summary_posterior, store = pub_store)
```

## Convergence Check {#sec-convergence}

```{r}
#| label: convergence-check

targets::tar_read(tumor_ssls_convergence_posterior, store = pub_store)
```

## Detailed Inspection {#sec-detailed}

For detailed parameter-level diagnostics, trace plots, and pair plots, load the fit object:

```r
fit <- targets::tar_read(tumor_ssls_res_posterior, store = pub_store)
fit$summary()
fit$diagnostic_summary()
```
```

- [ ] **Step 3: Commit**

```bash
git add quarto/publication/website/analysis/data-exploration.qmd \
        quarto/publication/website/analysis/diagnostics.qmd
git commit -m "feat(publication): add data-exploration and diagnostics pages"
```

---

## Task 9: Create analysis pages — Survival and Overall Survival

**Files:**
- Create: `quarto/publication/website/analysis/survival-analysis.qmd`
- Create: `quarto/publication/website/analysis/overall-survival.qmd`

These call `plot_km` directly — the low-level primitive that takes generic data.

- [ ] **Step 1: Create `survival-analysis.qmd`**

```qmd
---
title: "Survival Analysis"
subtitle: "Kaplan-Meier forecasts and observed PFS"
---

```{r}
#| label: setup
#| include: false

source("_setup.R")
km_data     <- targets::tar_read(all_tumor_ssls_km_rvar, store = pub_store)
obs_km_data <- targets::tar_read(km_trial_pfs,           store = pub_store)
```

## PFS Kaplan-Meier Curves {#sec-pfs-km}

::: {.panel-tabset}

### Sample-Based

```{r}
#| label: fig-pfs-km-sample
#| fig-cap: Sample-based PFS Kaplan-Meier curves (posterior). Shaded bands show 50% and 90% credible intervals.

km_data |>
  filter(fit_type == "posterior") |>
  plot_km(
    km_est_var = sample_pfs_km_est,
    obs_km_data = obs_km_data |>
      filter(fct_match(btype, "ub")) |>
      select(!quantiles) |>
      unnest(km_data)
  )
```

### Superpopulation

```{r}
#| label: fig-pfs-km-spop
#| fig-cap: Superpopulation PFS Kaplan-Meier curves (posterior).

km_data |>
  filter(fit_type == "posterior") |>
  plot_km(
    km_est_var = spop_pfs_km_est,
    obs_km_data = obs_km_data |>
      filter(fct_match(btype, "ub")) |>
      select(!quantiles) |>
      unnest(km_data)
  )
```

### Prior vs Posterior

```{r}
#| label: fig-pfs-km-prior-post
#| fig-cap: Prior vs posterior PFS KM comparison.

km_data |>
  plot_km(
    km_est_var = sample_pfs_km_est,
    group = fit_type,
    alpha_group = fit_type
  )
```

:::
```

- [ ] **Step 2: Create `overall-survival.qmd`**

```qmd
---
title: "Overall Survival"
subtitle: "OS Kaplan-Meier forecasts and observed data"
---

```{r}
#| label: setup
#| include: false

source("_setup.R")
os_km_data  <- targets::tar_read(all_tumor_ssls_km_os_rvar, store = pub_store)
obs_os_data <- targets::tar_read(km_trial_os,               store = pub_store)
```

## OS Kaplan-Meier Curves {#sec-os-km}

::: {.panel-tabset}

### Sample-Based

```{r}
#| label: fig-os-km-sample
#| fig-cap: Sample-based OS Kaplan-Meier curves (posterior).

os_km_data |>
  filter(fit_type == "posterior") |>
  plot_km(
    km_est_var = sample_os_km_est,
    obs_km_data = obs_os_data |>
      filter(fct_match(btype, "ub")) |>
      select(!quantiles) |>
      unnest(km_data)
  )
```

### Superpopulation

```{r}
#| label: fig-os-km-spop
#| fig-cap: Superpopulation OS Kaplan-Meier curves (posterior).

os_km_data |>
  filter(fit_type == "posterior") |>
  plot_km(
    km_est_var = spop_os_km_est,
    obs_km_data = obs_os_data |>
      filter(fct_match(btype, "ub")) |>
      select(!quantiles) |>
      unnest(km_data)
  )
```

### Prior vs Posterior

```{r}
#| label: fig-os-km-prior-post
#| fig-cap: Prior vs posterior OS KM comparison.

os_km_data |>
  plot_km(km_est_var = sample_os_km_est, group = fit_type, alpha_group = fit_type)
```

:::
```

- [ ] **Step 3: Commit**

```bash
git add quarto/publication/website/analysis/survival-analysis.qmd \
        quarto/publication/website/analysis/overall-survival.qmd
git commit -m "feat(publication): add survival and OS analysis pages"
```

---

## Task 10: Create analysis pages — Clinical Endpoints

**Files:**
- Create: `quarto/publication/website/analysis/clinical-endpoints.qmd`

ORR and PFS-n call `plot_outcome_by_pdl1_and_trial` and `plot_survival_n_by_estimate_type` — but those hardcode `fct_match(dco, "jan26")` and PDL1 group filtering. Publication calls `stat_pointinterval` directly on the simpler data structure.

- [ ] **Step 1: Create `clinical-endpoints.qmd`**

```qmd
---
title: "Clinical Endpoints"
subtitle: "Overall Response Rate and Progression-Free Survival results"
---

```{r}
#| label: setup
#| include: false

source("_setup.R")
orr_data      <- targets::tar_read(all_tumor_ssls_orr_rvar,                    store = pub_store)
pfs_quant     <- targets::tar_read(all_tumor_ssls_trial_pfs_quant,              store = pub_store)
pfs_n_data    <- targets::tar_read(all_tumor_ssls_forecast_target_pfs_n_rvar,   store = pub_store)
pfs_timepoints <- targets::tar_read(pfs_timepoints,                             store = pub_store)
obs_km_data   <- targets::tar_read(km_trial_pfs,                                store = pub_store)
```

## Overall Response Rate (ORR) {#sec-orr}

::: {.panel-tabset}

### Sample-Based

```{r}
#| label: fig-orr-sample
#| fig-cap: Sample-based ORR (posterior).

orr_data |>
  filter(fit_type == "posterior") |>
  ggplot() +
  stat_pointinterval(aes(xdist = sample_target_orr, y = trial, color = fit_type),
                     .width = c(0.5, 0.9)) +
  scale_color_discrete("", type = AZ_palette, label = str_to_title) +
  scale_x_continuous("Overall Response Rate (ORR)", limits = c(0, 1)) +
  labs(y = "") +
  NULL
```

### Prior vs Posterior

```{r}
#| label: fig-orr-prior-post
#| fig-cap: Prior vs posterior ORR comparison.

orr_data |>
  ggplot() +
  stat_pointinterval(aes(xdist = sample_target_orr, y = fit_type, color = fit_type),
                     .width = c(0.5, 0.9)) +
  scale_color_discrete("", type = AZ_palette, label = str_to_title) +
  scale_x_continuous("Overall Response Rate (ORR)", limits = c(0, 1)) +
  labs(y = "") +
  NULL
```

:::

## Median PFS {#sec-median-pfs}

::: {.panel-tabset}

### Sample-Based

```{r}
#| label: fig-median-pfs-sample
#| fig-cap: Sample-based median PFS (posterior). Weeks converted to months for display.

pfs_quant |>
  filter(fit_type == "posterior", quantile == 0.5) |>
  ggplot() +
  stat_pointinterval(aes(xdist = sample_pfs_quant, y = trial, color = fit_type),
                     .width = c(0.5, 0.9)) +
  scale_x_continuous("Median PFS [Months]",
                     breaks = months_to_weeks(seq(0, 30, 3)),
                     label  = label_weeks_to_months) +
  scale_color_discrete("", type = AZ_palette, label = str_to_title) +
  labs(y = "") +
  NULL
```

### Prior vs Posterior

```{r}
#| label: fig-median-pfs-prior-post
#| fig-cap: Prior vs posterior median PFS comparison.

pfs_quant |>
  filter(quantile == 0.5) |>
  ggplot() +
  stat_pointinterval(aes(xdist = sample_pfs_quant, y = fit_type, color = fit_type),
                     .width = c(0.5, 0.9)) +
  scale_x_continuous("Median PFS [Months]",
                     breaks = months_to_weeks(seq(0, 30, 3)),
                     label  = label_weeks_to_months) +
  scale_color_discrete("", type = AZ_palette, label = str_to_title) +
  labs(y = "") +
  NULL
```

:::

## PFS at Fixed Timepoints {#sec-pfs-n}

```{r}
#| label: fig-pfs-n
#| fig-cap: PFS probability at fixed timepoints (posterior, sample-based).

pfs_n_data |>
  filter(fit_type == "posterior") |>
  ggplot() +
  stat_pointinterval(
    aes(ydist = sample_pfs_n, x = timepoint, color = fit_type,
        group = interaction(fit_type, timepoint)),
    position = "dodge", .width = c(0.5, 0.8)
  ) +
  scale_x_continuous("Timepoint [Months]", breaks = pfs_timepoints$timepoint) +
  scale_y_continuous("PFS Probability", limits = c(0, 1)) +
  scale_color_discrete("", type = AZ_palette, label = str_to_title) +
  NULL
```
```

- [ ] **Step 2: Commit**

```bash
git add quarto/publication/website/analysis/clinical-endpoints.qmd
git commit -m "feat(publication): add clinical endpoints page"
```

---

## Task 11: Create analysis pages — Parameter Inference

**Files:**
- Create: `quarto/publication/website/analysis/parameter-inference.qmd`

Uses `plot_prior_post_dens`, `plot_level_rates`, `plot_level_decrease_prop` — all generic primitives.

- [ ] **Step 1: Create `parameter-inference.qmd`**

```qmd
---
title: "Parameter Inference"
subtitle: "Population and patient-level posterior distributions"
---

```{r}
#| label: setup
#| include: false

source("_setup.R")
rates_rvar   <- targets::tar_read(all_tumor_ssls_rates_rvar,            store = pub_store)
noise_sd     <- targets::tar_read(all_tumor_ssls_noise_sd_rvar,         store = pub_store)
patient_rates <- targets::tar_read(all_tumor_ssls_patient_rates_bpi,     store = pub_store)
decrease_prop <- targets::tar_read(all_tumor_ssls_patient_decrease_prop_bpi, store = pub_store)
coef_data    <- targets::tar_read(all_tumor_ssls_coef,                  store = pub_store)
```

## Population-Level Rates {#sec-pop-rates}

```{r}
#| label: fig-pop-rates
#| fig-cap: Population-level log growth and regression rates — prior vs posterior.

rates_rvar |>
  plot_prior_post_dens(.value, normalize = "groups") +
  labs(x = "") +
  facet_wrap(vars(.variable), scales = "free",
             labeller = labeller(.variable = \(l) str_remove(l, "log_"))) +
  NULL
```

## Measurement Noise {#sec-noise}

```{r}
#| label: fig-noise-sd
#| fig-cap: Measurement noise SD (SLD observation error) — prior vs posterior.

noise_sd |>
  plot_prior_post_hist(normalize = "groups") +
  labs(x = "Measurement SD") +
  NULL
```

## Patient-Level Rate Distributions {#sec-patient-rates}

```{r}
#| label: fig-patient-rates
#| fig-cap: Distribution of patient-level posterior growth and regression rates.

patient_rates |>
  filter(fit_type == "posterior") |>
  plot_level_rates()
```

## Growth Fraction {#sec-growth-fraction}

```{r}
#| label: fig-decrease-prop
#| fig-cap: Distribution of patient-level initial growth fraction (posterior).

decrease_prop |>
  filter(fit_type == "posterior") |>
  plot_level_decrease_prop()
```

## Covariate Effects {#sec-covariate-effects}

```{r}
#| label: fig-covariate-effects
#| fig-cap: Covariate coefficient posterior distributions (QR-space).

coef_data |>
  filter(fit_type == "posterior") |>
  ggplot() +
  stat_pointinterval(aes(xdist = .value, y = factor(n), color = .name),
                     .width = c(0.5, 0.9)) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  scale_color_discrete("Parameter", type = AZ_palette) +
  labs(x = "Coefficient (QR space)", y = "Covariate index") +
  NULL
```
```

- [ ] **Step 2: Commit**

```bash
git add quarto/publication/website/analysis/parameter-inference.qmd
git commit -m "feat(publication): add parameter inference page"
```

---

## Task 12: Create analysis pages — Patient Dynamics

**Files:**
- Create: `quarto/publication/website/analysis/patient-dynamics.qmd`

Uses `plot_sld_forecast_censored` and `plot_recist_censored` — both generic primitives.

- [ ] **Step 1: Create `patient-dynamics.qmd`**

The page reads individual per-event-type targets directly (same pattern as sclc's patient-dynamics page — `all_tumor_ssls_staged_sld_rvar` doesn't carry an event_type column, so individual targets are used).

```qmd
---
title: "Patient Dynamics"
subtitle: "Individual tumor trajectories and RECIST predictions"
---

```{r}
#| label: setup
#| include: false

source("_setup.R")
sld_rc    <- targets::tar_read(tumor_ssls_staged_sld_rvar_right_censored_posterior, store = pub_store)
sld_unc   <- targets::tar_read(tumor_ssls_staged_sld_rvar_uncensored_posterior,     store = pub_store)
recist_rc <- targets::tar_read(tumor_ssls_staged_recist_rvar_right_censored_posterior, store = pub_store)
```

## SLD Forecasts — Right-Censored Patients {#sec-sld-censored}

Individual SLD trajectories for patients who are right-censored (progression-free at data cutoff), with model forecasts extended beyond the last observed measurement.

```{r}
#| label: fig-sld-censored
#| fig-cap: SLD trajectories and forecasts for right-censored patients (posterior).
#| fig-width: 10
#| fig-height: 12

sld_rc |> plot_sld_forecast_censored()
```

## SLD Posterior Predictive Check — Uncensored Patients {#sec-sld-ppc}

```{r}
#| label: fig-sld-ppc
#| fig-cap: Observed vs posterior-predicted SLD for uncensored patients.
#| fig-width: 10
#| fig-height: 12

sld_unc |> plot_sld_forecast_censored()
```

## RECIST Predictions — Right-Censored Patients {#sec-recist-censored}

```{r}
#| label: fig-recist-censored
#| fig-cap: RECIST response predictions for right-censored patients (posterior).
#| fig-width: 10
#| fig-height: 12

recist_rc |> plot_recist_censored()
```
```

- [ ] **Step 2: Commit**

```bash
git add quarto/publication/website/analysis/patient-dynamics.qmd
git commit -m "feat(publication): add patient dynamics page"
```

---

## Task 13: Update `.claude/rules/r-guidelines.md` with publication site commands

**Files:**
- Modify: `.claude/rules/r-guidelines.md`

- [ ] **Step 1: Add publication site entries to the r-guidelines.md**

In `.claude/rules/r-guidelines.md`, find the section under `### Building the Website Locally` and add after the pioneer lines:

```bash
# From project root — Publication
quarto render quarto/publication/website
```

Find the section `### Publishing to RStudio Connect` and add after the pioneer deploySite line:

```r
Rscript -e 'rsconnect::deploySite(siteDir = "quarto/publication/website", server = "az-connect", account = "kmjq089")'
```

Find the `### Key Configuration Files` section and add:
```
- `quarto/publication/website/az-theme.scss` — Publication site theme (same as pioneer)
```

- [ ] **Step 2: Commit**

```bash
git add .claude/rules/r-guidelines.md
git commit -m "docs: add publication website build and publish commands to r-guidelines"
```

---

## Task 14: Smoke-test the site renders

- [ ] **Step 1: Check TAR_PROJECT is set correctly in `.Renviron`**

```bash
grep TAR_PROJECT .Renviron
```
Expected: `TAR_PROJECT=publication`

- [ ] **Step 2: Verify the pub store exists with completed targets**

```bash
ls /mnt/data/analysis-results/${DOMINO_STARTING_USERNAME}/publication/*/  2>/dev/null | head -20
```

- [ ] **Step 3: Render the site (from project root)**

```bash
quarto render quarto/publication/website 2>&1 | tail -30
```
Expected: `Output created: _site/index.html` with no fatal errors. Freeze errors for missing targets are acceptable at this stage — note which pages fail.

- [ ] **Step 4: Commit any freeze files generated**

```bash
git add quarto/publication/website/_freeze/
git commit -m "chore(publication): add initial freeze cache from smoke-test render"
```
