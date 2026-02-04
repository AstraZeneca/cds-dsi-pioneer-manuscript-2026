## Description

<!-- Briefly describe what this PR does and why -->

## Changes

<!-- List the main changes in this PR -->
-
-

## Testing

<!-- Describe how you tested these changes -->
- [ ] Code runs without errors
- [ ] Tests pass (if applicable)
- [ ] Manually verified changes work as expected

---

## Pre-Merge Checklist

**IMPORTANT:** All items below are mandatory. Do not submit for review until all boxes are checked.

### Critical Rules (Will Block Merge)

- [ ] **NO hardcoded subject IDs** - Verified no patient IDs (usubjid, E0302003, etc.) in:
  - [ ] Production code
  - [ ] Debug/example code
  - [ ] Comments
  - [ ] Test fixtures
  - Used dynamic selection instead: `slice_sample()`, `filter()`, or parameterized IDs

- [ ] **NO install.packages() calls** - Removed all runtime package installation
  - Dependencies added to `renv.lock` via `renv::snapshot()` instead

- [ ] **Use native pipe `|>`** - Replaced all `%>%` with `|>` throughout

- [ ] **NO syntax errors** - Verified code parses:
  ```bash
  # R code
  Rscript -e 'source("r/your_file.R")'

  # Stan code (if applicable)
  stanc --include-paths=stan,stan/ssls your_model.stan
  ```

### Code Quality

- [ ] Follows tidyverse style guide
- [ ] No global assignments (`<<-`) in functions
- [ ] Function parameter defaults don't reference undefined variables
- [ ] Documentation matches implementation (roxygen, function returns, etc.)
- [ ] No `tar_config_set()` - used explicit `store=` arguments

### Pre-Submission Validation

- [ ] **Ran automated code review** (highly recommended):
  ```bash
  # In Claude Code CLI
  /pr-review-toolkit:review-pr code
  ```

- [ ] **Searched for common violations**:
  ```bash
  # Check for hardcoded IDs
  grep -rE "(E[0-9]{7,}|[A-Z][0-9]{4}[0-9]{3}[0-9]{3})" r/ --include="*.R"

  # Check for install.packages
  grep -r "install.packages" r/

  # Check for old pipe operator
  grep -r "%>%" r/
  ```

---

## Reviewer Notes

<!-- Any specific areas you'd like reviewers to focus on? -->

---

<details>
<summary>See full guidelines in CLAUDE.md</summary>

For complete coding standards and PR requirements, see [CLAUDE.md](../CLAUDE.md#pull-request-checklist).

</details>
