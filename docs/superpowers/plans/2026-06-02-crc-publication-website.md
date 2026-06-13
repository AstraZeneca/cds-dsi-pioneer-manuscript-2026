# CRC Publication Website Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Split `quarto/publication/website/` into two independent disease subsites (`sclc/` and `crc/`) with a shared `pub_tar_read()` helper that handles the `tar_map` disease suffix transparently.

**Architecture:** Rename the existing `website/` directory to `sclc/`. Create `quarto/publication/_setup_common.R` with a `pub_tar_read()` helper that appends `pub_disease_suffix` to target names. Create a parallel `crc/` subsite (copied from `sclc/`) with stub pages and `pub_disease_suffix <- "crc"`. SCLC uses `pub_disease_suffix <- ""` to keep the existing store working without a rebuild.

**Tech Stack:** Quarto, R, `targets`, `here`, `rsconnect`

---

## File Map

| Action | Path | Purpose |
|--------|------|---------|
| Rename dir | `quarto/publication/website/` → `quarto/publication/sclc/` | SCLC subsite |
| Modify | `quarto/publication/sclc/_setup.R` | Set `pub_disease_suffix <- ""`; source `_setup_common.R` |
| Modify | `quarto/publication/sclc/_quarto.yml` | Update title to include "SCLC" |
| Modify | `quarto/publication/sclc/index.qmd` | Update hero text |
| Modify | `quarto/publication/sclc/analysis/*.qmd` (all 7) | Replace `tar_read()` calls with `pub_tar_read()` |
| Create | `quarto/publication/_setup_common.R` | `pub_tar_read()` helper |
| Create | `quarto/publication/crc/` (full subsite) | CRC subsite — all files |
| Modify | `.claude/rules/r-guidelines.md` | Update render/deploy commands |

---

## Task 1: Create the shared `_setup_common.R` helper

**Files:**
- Create: `quarto/publication/_setup_common.R`

- [ ] **Step 1: Create the file**

```r
# quarto/publication/_setup_common.R
# pub_tar_read() reads a target from pub_store, appending pub_disease_suffix
# when non-empty. Set pub_disease_suffix before sourcing this file.
#
# pub_disease_suffix <- ""     → reads unsuffixed name (legacy SCLC store)
# pub_disease_suffix <- "crc"  → reads <name>_crc

pub_tar_read <- function(name, store = pub_store) {
  full_name <- if (nchar(pub_disease_suffix) > 0) {
    paste0(name, "_", pub_disease_suffix)
  } else {
    name
  }
  targets::tar_read(!!rlang::sym(full_name), store = store)
}
```

- [ ] **Step 2: Verify the file exists**

```bash
ls quarto/publication/_setup_common.R
```

Expected: file listed.

- [ ] **Step 3: Commit**

```bash
git add quarto/publication/_setup_common.R
git commit -m "feat(pub-website): add shared pub_tar_read() helper"
```

---

## Task 2: Rename `website/` to `sclc/`

**Files:**
- Rename: `quarto/publication/website/` → `quarto/publication/sclc/`

- [ ] **Step 1: Rename the directory with git mv**

```bash
git mv quarto/publication/website quarto/publication/sclc
```

- [ ] **Step 2: Verify**

```bash
ls quarto/publication/
```

Expected: `sclc/` present, `website/` gone.

- [ ] **Step 3: Commit**

```bash
git commit -m "feat(pub-website): rename website/ to sclc/"
```

---

## Task 3: Update SCLC `_setup.R` to use the helper

**Files:**
- Modify: `quarto/publication/sclc/_setup.R`

- [ ] **Step 1: Replace the file content**

Current content:
```r
here::i_am("quarto/publication/website/_setup.R")

source(here::here("quarto", "_shared", "_setup.R"))

source(here::here("r", "sclc", "plot_functions.R"))
source(here::here("r", "sclc", "table_functions.R"))
source(here::here("r", "sclc", "accuracy.R"))

source(here::here("r", "publication", "target_stores.R"))
```

Replace with:
```r
pub_disease_suffix <- ""  # blank: SCLC store predates tar_map disease suffix

here::i_am("quarto/publication/sclc/_setup.R")

source(here::here("quarto", "_shared", "_setup.R"))

source(here::here("r", "sclc", "plot_functions.R"))
source(here::here("r", "sclc", "table_functions.R"))
source(here::here("r", "sclc", "accuracy.R"))

source(here::here("r", "publication", "target_stores.R"))
source(here::here("quarto", "publication", "_setup_common.R"))
```

- [ ] **Step 2: Commit**

```bash
git add quarto/publication/sclc/_setup.R
git commit -m "feat(pub-website): update SCLC _setup.R — add pub_disease_suffix and helper"
```

---

## Task 4: Update SCLC `_quarto.yml` title

**Files:**
- Modify: `quarto/publication/sclc/_quarto.yml`

- [ ] **Step 1: Change the website title**

Find this line in `_quarto.yml`:
```yaml
  title: "PIONEER: Publication"
```

Replace with:
```yaml
  title: "PIONEER: Publication — SCLC"
```

Also find and update the page footer center text:
```yaml
    center: "PIONEER: Publication"
```
→
```yaml
    center: "PIONEER: Publication — SCLC"
```

- [ ] **Step 2: Commit**

```bash
git add quarto/publication/sclc/_quarto.yml
git commit -m "feat(pub-website): set SCLC subsite title"
```

---

## Task 5: Update SCLC `index.qmd` hero text

**Files:**
- Modify: `quarto/publication/sclc/index.qmd`

- [ ] **Step 1: Update the hero heading and description**

Find:
```markdown
# PIONEER Publication Analysis

Bayesian hierarchical model for tumor dynamics and progression-free survival — publication dataset (497 patients, 6 covariates).
```

Replace with:
```markdown
# PIONEER Publication — SCLC

Bayesian hierarchical model for tumor dynamics and progression-free survival — SCLC publication dataset (497 patients, 6 covariates).
```

Also find the "About This Analysis" paragraph:
```markdown
This website presents results from the publication pipeline running the SF-SSM Stan model on 497 patients with 6 covariates (age, sex, ECOG, haemoglobin, log-LDH, albumin).
```

Replace with:
```markdown
This website presents SCLC results from the publication pipeline running the SF-SSM Stan model on 497 patients with 6 covariates (age, sex, ECOG, haemoglobin, log-LDH, albumin).
```

- [ ] **Step 2: Commit**

```bash
git add quarto/publication/sclc/index.qmd
git commit -m "feat(pub-website): update SCLC index hero text"
```

---

## Task 6: Migrate SCLC analysis pages to `pub_tar_read()`

**Files:**
- Modify: `quarto/publication/sclc/analysis/clinical-endpoints.qmd`
- Modify: `quarto/publication/sclc/analysis/data-exploration.qmd`
- Modify: `quarto/publication/sclc/analysis/diagnostics.qmd`
- Modify: `quarto/publication/sclc/analysis/overall-survival.qmd`
- Modify: `quarto/publication/sclc/analysis/parameter-inference.qmd`
- Modify: `quarto/publication/sclc/analysis/patient-dynamics.qmd`
- Modify: `quarto/publication/sclc/analysis/survival-analysis.qmd`

Replace every `targets::tar_read(<name>, store = pub_store)` with
`pub_tar_read("<name>")` in the setup chunk of each page.
The `store = pub_store` argument is dropped because it is the default in
`pub_tar_read()`.

- [ ] **Step 1: Update `clinical-endpoints.qmd` setup chunk**

```r
source("_setup.R")
orr_data       <- pub_tar_read("all_tumor_ssls_orr_rvar")
pfs_quant      <- pub_tar_read("all_tumor_ssls_trial_pfs_quant")
pfs_n_data     <- pub_tar_read("all_tumor_ssls_forecast_target_pfs_n_rvar")
pfs_timepoints <- tibble::tibble(n = 1:5, timepoint = c(6, 9, 12, 15, 18))
obs_km_data    <- pub_tar_read("km_trial_pfs")
```

- [ ] **Step 2: Update `data-exploration.qmd` setup chunk**

```r
source("_setup.R")
all_analysis_data <- pub_tar_read("all_analysis_data")
```

- [ ] **Step 3: Update `diagnostics.qmd` setup chunk and inline calls**

Setup chunk (no reads at setup time — reads are inline):
```r
source("_setup.R")
```

Replace the two inline `tar_read` calls in code chunks:
```r
# chunk: diagnostics-summary
pub_tar_read("tumor_ssls_nuts_summary_posterior")
```
```r
# chunk: convergence-check
pub_tar_read("tumor_ssls_convergence_posterior")
```

The inline code example block (inside a fenced ` ```r ` block, not executed)
should also be updated for consistency:
```r
fit <- pub_tar_read("tumor_ssls_res_posterior")
fit$summary()
fit$diagnostic_summary()
```

- [ ] **Step 4: Update `overall-survival.qmd` setup chunk**

```r
source("_setup.R")
os_km_data  <- pub_tar_read("all_tumor_ssls_km_os_rvar")
obs_os_data <- pub_tar_read("km_trial_os") |>
  filter(fct_match(btype, "ub")) |>
  select(!quantiles) |>
  unnest(km_data)
```

- [ ] **Step 5: Update `parameter-inference.qmd` setup chunk**

```r
source("_setup.R")
rates_rvar    <- pub_tar_read("all_tumor_ssls_rates_rvar")
noise_sd      <- pub_tar_read("all_tumor_ssls_noise_sd_rvar")
patient_rates <- pub_tar_read("all_tumor_ssls_patient_rates_bpi")
decrease_prop <- pub_tar_read("all_tumor_ssls_patient_decrease_prop_bpi")
coef_data     <- pub_tar_read("all_tumor_ssls_coef")
```

- [ ] **Step 6: Update `patient-dynamics.qmd` setup chunk**

```r
source("_setup.R")
sld_rc    <- pub_tar_read("tumor_ssls_staged_sld_rvar_right_censored_posterior")
sld_unc   <- pub_tar_read("tumor_ssls_staged_sld_rvar_uncensored_posterior")
recist_rc <- pub_tar_read("tumor_ssls_staged_recist_rvar_right_censored_posterior")
```

- [ ] **Step 7: Update `survival-analysis.qmd` setup chunk**

```r
source("_setup.R")
km_data     <- pub_tar_read("all_tumor_ssls_km_rvar")
obs_km_data <- pub_tar_read("km_trial_pfs") |>
  filter(fct_match(btype, "ub")) |>
  select(!quantiles) |>
  unnest(km_data)
```

- [ ] **Step 8: Commit all page changes**

```bash
git add quarto/publication/sclc/analysis/
git commit -m "feat(pub-website): migrate SCLC pages to pub_tar_read() helper"
```

---

## Task 7: Verify SCLC subsite renders

- [ ] **Step 1: Render the SCLC subsite from the project root**

```bash
quarto render quarto/publication/sclc
```

Expected: renders without error; all pages appear in `quarto/publication/sclc/_site/`.

- [ ] **Step 2: Commit freeze cache updates if any pages re-executed**

```bash
git add quarto/publication/sclc/_freeze/
git commit -m "chore(pub-website): update SCLC freeze cache after helper migration"
```

(Skip if no freeze files changed.)

---

## Task 8: Create CRC subsite structure

**Files:**
- Create: `quarto/publication/crc/` (full directory tree)

- [ ] **Step 1: Copy the SCLC subsite as the CRC template**

```bash
cp -r quarto/publication/sclc quarto/publication/crc
```

- [ ] **Step 2: Remove SCLC freeze cache from CRC (CRC has no data yet)**

```bash
rm -rf quarto/publication/crc/_freeze
```

- [ ] **Step 3: Remove the existing SCLC Connect deployment manifest from CRC**

```bash
rm -f quarto/publication/crc/rsconnect/
```

- [ ] **Step 4: Verify the CRC directory**

```bash
ls quarto/publication/crc/
```

Expected: `_quarto.yml  _setup.R  analysis/  az-theme.scss  images/  index.qmd  styles.css`

- [ ] **Step 5: Commit the skeleton**

```bash
git add quarto/publication/crc/
git commit -m "feat(pub-website): add CRC subsite skeleton (copied from sclc)"
```

---

## Task 9: Configure CRC `_setup.R`

**Files:**
- Modify: `quarto/publication/crc/_setup.R`

- [ ] **Step 1: Replace the file content**

```r
pub_disease_suffix <- "crc"

here::i_am("quarto/publication/crc/_setup.R")

source(here::here("quarto", "_shared", "_setup.R"))

source(here::here("r", "sclc", "plot_functions.R"))
source(here::here("r", "sclc", "table_functions.R"))
source(here::here("r", "sclc", "accuracy.R"))

source(here::here("r", "publication", "target_stores.R"))
source(here::here("quarto", "publication", "_setup_common.R"))
```

- [ ] **Step 2: Commit**

```bash
git add quarto/publication/crc/_setup.R
git commit -m "feat(pub-website): configure CRC _setup.R with pub_disease_suffix=crc"
```

---

## Task 10: Update CRC `_quarto.yml` title

**Files:**
- Modify: `quarto/publication/crc/_quarto.yml`

- [ ] **Step 1: Change the website title**

Find:
```yaml
  title: "PIONEER: Publication — SCLC"
```

Replace with:
```yaml
  title: "PIONEER: Publication — CRC"
```

Also update the page footer:
```yaml
    center: "PIONEER: Publication — SCLC"
```
→
```yaml
    center: "PIONEER: Publication — CRC"
```

- [ ] **Step 2: Commit**

```bash
git add quarto/publication/crc/_quarto.yml
git commit -m "feat(pub-website): set CRC subsite title"
```

---

## Task 11: Write CRC stub pages

**Files:**
- Modify: `quarto/publication/crc/index.qmd`
- Modify: `quarto/publication/crc/analysis/clinical-endpoints.qmd`
- Modify: `quarto/publication/crc/analysis/data-exploration.qmd`
- Modify: `quarto/publication/crc/analysis/diagnostics.qmd`
- Modify: `quarto/publication/crc/analysis/overall-survival.qmd`
- Modify: `quarto/publication/crc/analysis/parameter-inference.qmd`
- Modify: `quarto/publication/crc/analysis/patient-dynamics.qmd`
- Modify: `quarto/publication/crc/analysis/survival-analysis.qmd`

- [ ] **Step 1: Update `crc/index.qmd` hero text**

Find and replace the hero heading and description:
```markdown
# PIONEER Publication — SCLC

Bayesian hierarchical model for tumor dynamics and progression-free survival — SCLC publication dataset (497 patients, 6 covariates).
```
→
```markdown
# PIONEER Publication — CRC

Bayesian hierarchical model for tumor dynamics and progression-free survival — CRC publication dataset.
```

Find and replace the "About This Analysis" paragraph:
```markdown
This website presents SCLC results from the publication pipeline running the SF-SSM Stan model on 497 patients with 6 covariates (age, sex, ECOG, haemoglobin, log-LDH, albumin).
```
→
```markdown
This website presents CRC results from the publication pipeline running the SF-SSM Stan model. Results will appear here once the CRC fit (`tumor_ssls_res_posterior_crc`) completes.
```

- [ ] **Step 2: Replace all 7 analysis pages with stubs**

For each of the 7 pages (`clinical-endpoints.qmd`, `data-exploration.qmd`,
`diagnostics.qmd`, `overall-survival.qmd`, `parameter-inference.qmd`,
`patient-dynamics.qmd`, `survival-analysis.qmd`), replace the entire setup
chunk (lines starting from ` ```{r}` down to the closing ` ``` ` of the setup
chunk) with just `source("_setup.R")`, and replace the page body content with
the pending callout. Keep the YAML front matter unchanged.

Template for each page (substitute the existing title/subtitle front matter):

`clinical-endpoints.qmd`:
```qmd
---
title: "Clinical Endpoints"
subtitle: "Overall Response Rate and Progression-Free Survival results"
---

```{r}
#| label: setup
#| include: false

source("_setup.R")
```

> **CRC results pending** — the `tumor_ssls_res_posterior_crc` fit is currently
> running. This page will be updated once results are available.
```

`data-exploration.qmd`:
```qmd
---
title: "Data Exploration"
subtitle: "Overview of the publication analysis dataset"
---

```{r}
#| label: setup
#| include: false

source("_setup.R")
```

> **CRC results pending** — the `tumor_ssls_res_posterior_crc` fit is currently
> running. This page will be updated once results are available.
```

`diagnostics.qmd`:
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

> **CRC results pending** — the `tumor_ssls_res_posterior_crc` fit is currently
> running. This page will be updated once results are available.
```

`overall-survival.qmd`:
```qmd
---
title: "Overall Survival"
subtitle: "OS Kaplan-Meier forecasts and observed data"
---

```{r}
#| label: setup
#| include: false

source("_setup.R")
```

> **CRC results pending** — the `tumor_ssls_res_posterior_crc` fit is currently
> running. This page will be updated once results are available.
```

`parameter-inference.qmd`:
```qmd
---
title: "Parameter Inference"
subtitle: "Population and patient-level posterior distributions"
---

```{r}
#| label: setup
#| include: false

source("_setup.R")
```

> **CRC results pending** — the `tumor_ssls_res_posterior_crc` fit is currently
> running. This page will be updated once results are available.
```

`patient-dynamics.qmd`:
```qmd
---
title: "Patient Dynamics"
subtitle: "Individual tumor trajectories and RECIST predictions"
---

```{r}
#| label: setup
#| include: false

source("_setup.R")
```

> **CRC results pending** — the `tumor_ssls_res_posterior_crc` fit is currently
> running. This page will be updated once results are available.
```

`survival-analysis.qmd`:
```qmd
---
title: "Survival Analysis"
subtitle: "Kaplan-Meier forecasts and observed PFS"
---

```{r}
#| label: setup
#| include: false

source("_setup.R")
```

> **CRC results pending** — the `tumor_ssls_res_posterior_crc` fit is currently
> running. This page will be updated once results are available.
```

- [ ] **Step 3: Commit**

```bash
git add quarto/publication/crc/
git commit -m "feat(pub-website): write CRC stub pages with pending callout"
```

---

## Task 12: Verify CRC subsite renders

- [ ] **Step 1: Render the CRC subsite from the project root**

```bash
quarto render quarto/publication/crc
```

Expected: renders without error; all pages show the pending callout in
`quarto/publication/crc/_site/`.

- [ ] **Step 2: Commit freeze cache**

```bash
git add quarto/publication/crc/_freeze/
git commit -m "chore(pub-website): add CRC freeze cache (stub pages)"
```

---

## Task 13: Update documentation

**Files:**
- Modify: `.claude/rules/r-guidelines.md`

- [ ] **Step 1: Update the render commands in r-guidelines.md**

Find:
```bash
# From project root — Publication
quarto render quarto/publication/website
```

Replace with:
```bash
# From project root — Publication (SCLC)
quarto render quarto/publication/sclc

# From project root — Publication (CRC)
quarto render quarto/publication/crc
```

Find the deploy commands section:
```r
Rscript -e 'rsconnect::deploySite(siteDir = "quarto/publication/website", server = "az-connect", account = "kmjq089")'
```

Replace with:
```r
Rscript -e 'rsconnect::deploySite(siteDir = "quarto/publication/sclc", server = "az-connect", account = "kmjq089")'
Rscript -e 'rsconnect::deploySite(siteDir = "quarto/publication/crc",  server = "az-connect", account = "kmjq089")'
```

Update the Published URLs note:
```
- Publication: https://rstudio-connect.seml.scp.astrazeneca.net/content/fdee3ad4-616a-42b4-b258-47dfb14ba50f/
```
→
```
- Publication (SCLC): https://rstudio-connect.seml.scp.astrazeneca.net/content/fdee3ad4-616a-42b4-b258-47dfb14ba50f/
- Publication (CRC): (new Connect item — created on first deploy)
```

- [ ] **Step 2: Commit**

```bash
git add .claude/rules/r-guidelines.md
git commit -m "docs: update render/deploy commands for sclc/ and crc/ subsites"
```
