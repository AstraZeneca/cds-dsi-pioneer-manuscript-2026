---
paths:
  - "**/*.R"
  - "**/*.qmd"
---

## R Coding Guidelines

- Use modern pipe operator `|>` (not `%>%`)
- Follow tidyverse style guide
- Prefer `purrr` and `dplyr` over base R loops
- Use `testthat` for unit tests
- **NEVER hardcode subject IDs** (usubjid, patient_id, etc.) - always use dynamic selection or filtering
- **NO hardcoded subject IDs in GitHub issues** — describe the problem (e.g., "5 patients have NA pfs") without listing specific IDs. Add diagnostic R code as a comment on the issue instead.

## Quarto Website

### Location and Purpose

Each trial has its own Quarto website under `quarto/<trial>/website/`:
- **SCLC-01**: `quarto/sclc/website/` — titled "PIONEER: SCLC-01 Forecasting"
- **Pioneer**: `quarto/pioneer/website/` — titled "PIONEER: Pioneer Forecasting"

Shared resources (theme, common setup) live in `quarto/_shared/`.

### Building the Website Locally

**IMPORTANT**: Always render from the project root directory (`/mnt/code`), not from inside the trial website directory. This is because `_quarto.yml` sets `execute-dir: project`, which means R code needs access to the project-level renv and data.

```bash
# From project root — SCLC-01
quarto render quarto/sclc/website

# From project root — Pioneer
quarto render quarto/pioneer/website

# From project root — Publication (SCLC)
quarto render quarto/publication/sclc

# From project root — Publication (CRC)
quarto render quarto/publication/crc
```

Rendered sites are output to `quarto/<trial>/website/_site/`. To preview:
```bash
quarto preview quarto/sclc/website
quarto preview quarto/pioneer/website
```

### Publishing to RStudio Connect

See `docs/PUBLISHING.md` for full publishing instructions (API keys, rsconnect setup, troubleshooting).

**Quick reference:**
```bash
# Render first (always from project root)
quarto render quarto/sclc/website
quarto render quarto/pioneer/website

# Deploy via R
Rscript -e 'rsconnect::deploySite(siteDir = "quarto/sclc/website", server = "az-connect", account = "kmjq089")'
Rscript -e 'rsconnect::deploySite(siteDir = "quarto/pioneer/website", server = "az-connect", account = "kmjq089")'
Rscript -e 'rsconnect::deploySite(siteDir = "quarto/publication/sclc", server = "az-connect", account = "kmjq089")'
Rscript -e 'rsconnect::deploySite(siteDir = "quarto/publication/crc",  server = "az-connect", account = "kmjq089")'
```

**Published URLs:**
- Pioneer: https://rstudio-connect.seml.scp.astrazeneca.net/content/4ec50122-5d33-43fa-831c-3df0243084f7/
- Publication (SCLC): https://rstudio-connect.seml.scp.astrazeneca.net/content/fdee3ad4-616a-42b4-b258-47dfb14ba50f/
- Publication (CRC): (new Connect item — created on first deploy)

### Key Configuration Files
- `_quarto.yml` - Site configuration, navigation, theme settings (per-trial)
- `quarto/_shared/az-colors.scss` - AstraZeneca color scheme (navy, gold, turquoise, etc.)
- `quarto/sclc/website/az-theme.scss` - SCLC-01 theme (extends az-colors.scss)
- `quarto/publication/sclc/az-theme.scss` — Publication SCLC site theme (same as pioneer)
- `quarto/publication/crc/az-theme.scss` — Publication CRC site theme
- `styles.css` - Custom CSS for hero section and layout
- `_freeze/` - Cache directory for executed R code (speeds up rebuilds)
- `quarto/_shared/images/` - PIONEER helmet logo

### Theme and Styling

The site uses **AZ corporate colors** defined in `quarto/_shared/az-colors.scss` (navy `#003865`, gold `#F0AB00`, turquoise `#68D2DF`, plum `#830051`, pink `#D0006F`, platinum `#9DB0AC`).

### Rebuilding After Changes

If you update plot functions in `r/plot_functions.R`, you may need to clear the freeze cache:
```bash
rm -rf quarto/sclc/website/_freeze/analysis/
```

Then re-render the affected pages.

### TikZ Diagrams
- Compile diagrams: `/home/ubuntu/.TinyTeX/bin/x86_64-linux/pdflatex -output-directory=quarto/sclc/website/images quarto/sclc/website/images/<name>.tex`
- Convert to SVG: `pdf2svg quarto/sclc/website/images/<name>.pdf quarto/sclc/website/images/<name>-tikz.svg`
- Clean artifacts: `rm -f quarto/sclc/website/images/<name>.{aux,log,pdf}`
- Use `/tikz-diagram` skill for guided workflow with prerequisite checks

### Documentation Sync
- Run `/sync-docs` after code changes to identify outdated documentation
- Updates doc-mappings.yaml when new source files are added to major features
- Architecture diagrams in `quarto/sclc/website/images/` may need regeneration when model structure changes

## Quarto and Documentation Conventions

- **Always use "SCLC-01"** when referring to the trial in user-facing text (documentation, plots, presentations)
- Use lowercase "sclc" only for code identifiers (variable names, trial codes, file paths)
- Example: Write "SCLC-01 trial" in figure captions, but `filter(trial == "sclc")` in R code
- **Use automatic section numbering**: Set `number-sections: true` in frontmatter, don't use manual numbers (1.1, 2.3) in headings
- **Cross-references**: Use section IDs `{#sec-name}` and reference with `@sec-name`, never hardcode "Section X.Y.Z"
- **Model specification is the blueprint**: The model specification is split across three pages, each authoritative for its domain:
  - `quarto/sclc/website/documentation/tumor-dynamics-specification.qmd` — tumor state-space model, RECIST, tumor priors
  - `quarto/sclc/website/documentation/multistate-specification.qmd` — illness-death model, transitions, multistate priors
  - `quarto/sclc/website/documentation/clinical-endpoints-specification.qmd` — PFS, OS, posterior inference
  - `quarto/sclc/website/documentation/model-architecture.qmd` — overview, hierarchy, notation, feature flags
  - When changing Stan code, update the relevant specification page to reflect the change
  - When the specification defines behavior (e.g., index conventions, endpoint formulas, routing logic), the code must not violate those definitions without updating the spec first
  - If a proposed code change contradicts the specification, flag the discrepancy before implementing
  - Treat the specification as a contract: it documents what the model *should* do, not just what it *happens* to do
