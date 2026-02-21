---
name: data-pipeline-confirmation
enabled: true
event: bash
pattern: Rscript.*/(data_preparation_pipeline|main_pipeline|wrangle_data|prepare.*data\.R|pipeline.*\.R)
action: block
---

⚠️ **Data Preparation Pipeline Script Detected**

You're about to run a data preparation or pipeline script.

**Detected command pattern:**
- Data preparation pipeline scripts
- Main pipeline execution scripts
- Data wrangling/cooking scripts

**Why this requires confirmation:**
- These scripts can take a long time to run
- They may overwrite existing processed data
- They should only be run when necessary (e.g., after raw data updates)

**Before proceeding, verify:**
- [ ] Do you really need to regenerate this data?
- [ ] Is there already a processed CSV file available?
- [ ] Have you checked if targets can rebuild just what's needed?

**Better alternatives:**
- Use `targets::tar_make()` to rebuild only what changed
- Check if existing CSV files already have the data you need
- Ask user before running long data processing scripts

If you really need to run this, please confirm with the user first.
