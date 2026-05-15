# Simple Survival Website Section — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add an "OS Survival Only Model" sidebar section to the Pioneer Quarto website with a single page showing Week-13 landmark OS KM curves and quantile estimates for the `re_arm` variant.

**Architecture:** Three coordinated changes — add `simple_survival_store` to `target_stores.R` (sourced by `_setup.R` automatically), add the new sidebar section and render entry to `_quarto.yml`, and create `analysis/simple-survival.qmd` that loads from that store and plots with existing `plot_pioneer_km()` and `plot_pioneer_quant_by_estimate_type()`.

**Tech Stack:** Quarto, R, targets, ggplot2, ggdist, patchwork. All plot functions already exist in `r/pioneer/plot_functions.R`.

---

## File Map

| Action | File |
|--------|------|
| Modify | `r/pioneer/target_stores.R` |
| Modify | `quarto/pioneer/website/_quarto.yml` |
| Create | `quarto/pioneer/website/analysis/simple-survival.qmd` |

No new R functions needed — all plotting is handled by existing helpers.

---

## Task 1: Add `simple_survival_store` to target stores

**Files:**
- Modify: `r/pioneer/target_stores.R`

- [ ] **Step 1: Add the store variable**

Append to the end of `r/pioneer/target_stores.R`:

```r
# Simple survival (landmark OS, 0->2 only) store.
# Override with TAR_SIMPLE_SURVIVAL_RUN env var (default: "simple-survival").
simple_survival_store <- file.path(
  output_path, "pioneer",
  Sys.getenv("TAR_SIMPLE_SURVIVAL_RUN", "simple-survival"),
  "_targets"
)
```

- [ ] **Step 2: Verify the file parses without error**

```bash
cd /mnt/code
Rscript -e '
  output_path <- "/mnt/data/analysis-results/karim_naguib"
  source("r/pioneer/target_stores.R")
  cat("simple_survival_store:", simple_survival_store, "\n")
'
```

Expected output:
```
simple_survival_store: /mnt/data/analysis-results/karim_naguib/pioneer/simple-survival/_targets
```

- [ ] **Step 3: Commit**

```bash
git add r/pioneer/target_stores.R
git commit -m "feat(website): add simple_survival_store to target_stores.R"
```

---

## Task 2: Register the page in `_quarto.yml`

**Files:**
- Modify: `quarto/pioneer/website/_quarto.yml`

- [ ] **Step 1: Add the render entry**

In the `render:` list (after the existing `analysis/` entries and before `documentation/`), add:

```yaml
    - "analysis/simple-survival.qmd"
```

The render list after the change (relevant portion):

```yaml
    - "analysis/subgroup-burden-endpoints.qmd"
    - "analysis/simple-survival.qmd"
    # - "analysis/model-validation.qmd"  # needs all_pioneer_ppc_rvar (not yet built)
    - "documentation/build-info.qmd"
```

- [ ] **Step 2: Add the sidebar section and navbar entry**

After the closing of the `full` sidebar (after line `# - section: "Model Validation" ...`), add a new sidebar entry. Also add a navbar entry after `"Full Model"`. The two changes:

**Navbar** — add after the `Full Model` entry:

```yaml
      - text: "OS Only Model"
        href: analysis/simple-survival.qmd
```

**Sidebar** — add after the `full` sidebar block (after its last `#` comment line):

```yaml
    - id: os-only
      title: "OS Survival Only Model"
      style: docked
      contents:
        - analysis/simple-survival.qmd
```

- [ ] **Step 3: Verify YAML parses**

```bash
cd /mnt/code
python3 -c "
import yaml
with open('quarto/pioneer/website/_quarto.yml') as f:
    d = yaml.safe_load(f)
sidebars = d['website']['sidebar']
names = [s.get('title') for s in sidebars]
print('Sidebars:', names)
nav = [e.get('text') for e in d['website']['navbar']['left'] if 'text' in e]
print('Navbar:', nav)
"
```

Expected output:
```
Sidebars: ['Data', 'Full Model', 'OS Survival Only Model', 'Documentation']
Navbar: ['Data', 'Full Model', 'OS Only Model', 'Documentation']
```

- [ ] **Step 4: Commit**

```bash
git add quarto/pioneer/website/_quarto.yml
git commit -m "feat(website): add OS Survival Only Model section to navigation"
```

---

## Task 3: Create the simple-survival page

**Files:**
- Create: `quarto/pioneer/website/analysis/simple-survival.qmd`

- [ ] **Step 1: Create the page**

```bash
cat > /mnt/code/.worktrees/karim/simple-survival/quarto/pioneer/website/analysis/simple-survival.qmd << 'QMDEOF'
---
title: "Overall Survival"
subtitle: "Week-13 landmark OS model (0→2 only, re_arm variant)"
---

```{r setup}
#| include: false
source("_setup.R")

km_os_data <- bind_rows(
  tar_read_raw("simple_survival_km_os_rvar_prior_re_arm",     store = simple_survival_store),
  tar_read_raw("simple_survival_km_os_rvar_posterior_re_arm", store = simple_survival_store)
)

quant_os_data <- bind_rows(
  tar_read_raw("simple_survival_quant_os_prior_re_arm",     store = simple_survival_store),
  tar_read_raw("simple_survival_quant_os_posterior_re_arm", store = simple_survival_store)
)

# Compute observed OS KM from the landmark cohort (times already rebased to Week 13).
# NOTE: we cannot use pioneer_obs_os_km here — that KM is from treatment start (t=0),
# whereas model predictions are from the Week-13 landmark. Instead we compute the
# conditional KM inline from simple_survival_analysis_data.
landmark_data <- tar_read(simple_survival_analysis_data, store = simple_survival_store)

obs_os_km <- landmark_data |>
  transmute(
    trial,
    t     = coalesce(death_week, patient_max_t),
    event = as.integer(!is.na(death_week))
  ) |>
  group_by(trial) |>
  group_modify(\(d, k) {
    fit <- survival::survfit(survival::Surv(t, event) ~ 1, data = d)
    tibble(t = fit$time, s = fit$surv)
  }) |>
  ungroup()
```

This page presents model-predicted overall survival (OS) Kaplan-Meier curves from
the Week-13 landmark analysis. Only the 0→2 (direct death) transition is modelled.

::: {.callout-note}
## Landmark Design

Patients are conditioned on surviving to Week 13 (Day 85, Cycle 4 Day 1). All event
times are rebased so that $t = 0$ corresponds to the landmark. First post-treatment PSA
(earliest measurement from Day 64 onwards) is included as an additional covariate.
The 0→1 progression transition is disabled; OS is driven solely by the 0→2 hazard.
The `re_arm` variant places arm-level random effects on the baseline 0→2 hazard.
:::

## Sample-Based OS {#sec-os-sample}

::: {.panel-tabset}

### Posterior

```{r}
#| label: fig-ss-os-sample-posterior
#| fig-cap: Posterior sample-based OS Kaplan-Meier curves by trial. Times are measured from the Week-13 landmark. Ribbons show 80% credible intervals.
#| fig-width: 8

km_os_data |>
  filter(fit_type == "posterior") |>
  plot_pioneer_km(km_est_var = sample_os_km_est, obs_km_data = obs_os_km, endpoint = "os")
```

### Prior vs Posterior

```{r}
#| label: fig-ss-os-sample-both
#| fig-cap: Prior and posterior sample-based OS Kaplan-Meier curves by trial. Prior is shown with lighter shading.
#| fig-width: 8

plot_pioneer_km(km_os_data, km_est_var = sample_os_km_est, obs_km_data = obs_os_km, endpoint = "os")
```

:::

## Superpopulation OS {#sec-os-spop}

Superpopulation estimates represent the expected survival for the underlying patient
population, averaging over individual patient characteristics.

::: {.panel-tabset}

### Posterior

```{r}
#| label: fig-ss-os-spop-posterior
#| fig-cap: Posterior superpopulation OS Kaplan-Meier curves by trial.
#| fig-width: 8

km_os_data |>
  filter(fit_type == "posterior") |>
  plot_pioneer_km(km_est_var = spop_os_km_est, obs_km_data = obs_os_km, endpoint = "os")
```

### Prior vs Posterior

```{r}
#| label: fig-ss-os-spop-both
#| fig-cap: Prior and posterior superpopulation OS Kaplan-Meier curves by trial.
#| fig-width: 8

plot_pioneer_km(km_os_data, km_est_var = spop_os_km_est, obs_km_data = obs_os_km, endpoint = "os")
```

:::

## OS Quantiles {#sec-os-quantiles}

Posterior estimates of the 25th, 50th, and 75th survival time quantiles (Q25, Q50, Q75)
by trial. Times are measured from the Week-13 landmark.

### Sample-Based

```{r}
#| label: fig-ss-os-quant-sample
#| fig-cap: Posterior sample-based OS time quantiles by trial. Points show posterior medians; bars show 50% (thick) and 80% (thin) credible intervals.
#| fig-width: 8

plot_pioneer_quant_by_estimate_type("sample", quant_os_data, endpoint = "os")
```

### Superpopulation

```{r}
#| label: fig-ss-os-quant-spop
#| fig-cap: Posterior superpopulation OS time quantiles by trial.
#| fig-width: 8

plot_pioneer_quant_by_estimate_type("spop", quant_os_data, endpoint = "os")
```

::: {.related-pages}
## Related Pages

- [Overall Survival](overall-survival.qmd) — Full model OS (all transitions)
- [Survival Analysis](survival-analysis.qmd) — Full model PFS
:::
QMDEOF
```

**Note:** The heredoc above will be interpreted by the shell. Alternatively, use the Write tool to create the file with the exact content shown.

- [ ] **Step 2: Verify target names and observed KM computation**

```bash
cd /mnt/code
Rscript -e '
  output_path <- "/mnt/data/analysis-results/karim_naguib"
  source("r/pioneer/target_stores.R")
  library(targets)
  library(dplyr)

  # Check all four data targets exist in the store
  meta <- tar_meta(store = simple_survival_store)
  targets_needed <- c(
    "simple_survival_km_os_rvar_prior_re_arm",
    "simple_survival_km_os_rvar_posterior_re_arm",
    "simple_survival_quant_os_prior_re_arm",
    "simple_survival_quant_os_posterior_re_arm",
    "simple_survival_analysis_data"
  )
  found <- targets_needed %in% meta$name
  for (i in seq_along(targets_needed)) cat(targets_needed[i], "->", found[i], "\n")

  # Check observed KM computes correctly
  landmark_data <- tar_read(simple_survival_analysis_data, store = simple_survival_store)
  obs_os_km <- landmark_data |>
    transmute(trial, t = coalesce(death_week, patient_max_t), event = as.integer(!is.na(death_week))) |>
    group_by(trial) |>
    group_modify(\(d, k) {
      fit <- survival::survfit(survival::Surv(t, event) ~ 1, data = d)
      tibble(t = fit$time, s = fit$surv)
    }) |>
    ungroup()
  cat("obs_os_km rows:", nrow(obs_os_km), "\n")
  cat("trials:", paste(unique(obs_os_km$trial), collapse = ", "), "\n")
  cat("t range:", min(obs_os_km$t), "to", max(obs_os_km$t), "weeks\n")
'
```

Expected: all five targets show `-> TRUE`, `obs_os_km` has >0 rows, t range starts near 0.

- [ ] **Step 3: Commit**

```bash
git add quarto/pioneer/website/analysis/simple-survival.qmd
git commit -m "feat(website): add simple-survival OS page (Week-13 landmark, re_arm)"
```

---

## Task 4: Render and verify

**Files:** No changes — render only.

- [ ] **Step 1: Render the new page in isolation**

```bash
cd /mnt/code
quarto render quarto/pioneer/website/analysis/simple-survival.qmd
```

Expected: renders without error, produces `_site/analysis/simple-survival.html`.

Check for common failure modes in the output:
- `Error in tar_read_raw` → target name wrong or store path wrong
- `object 'simple_survival_store' not found` → `target_stores.R` edit missing or not sourced
- `could not find function "plot_pioneer_quant_by_estimate_type"` → `_setup.R` not sourcing pioneer plot functions

- [ ] **Step 2: Spot-check rendered HTML**

```bash
grep -c "fig-ss-os" quarto/pioneer/website/_site/analysis/simple-survival.html
```

Expected: `4` (one match per figure label).

- [ ] **Step 3: Render the full site**

```bash
cd /mnt/code
quarto render quarto/pioneer/website
```

Expected: completes without error. The `_site/` directory will contain the updated site with the new "OS Survival Only Model" section visible in the sidebar.

- [ ] **Step 4: Commit freeze cache**

After a successful render, Quarto writes computed outputs to `_freeze/`. Commit them so subsequent renders skip re-execution:

```bash
git add quarto/pioneer/website/_freeze/analysis/simple-survival/
git add quarto/pioneer/website/_site/
git commit -m "chore(website): render simple-survival OS page"
```
