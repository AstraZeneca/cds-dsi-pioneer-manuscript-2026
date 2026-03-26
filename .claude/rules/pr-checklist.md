## Pull Request Checklist

Before creating or submitting a PR, verify all items below. These are **mandatory requirements**:

### Critical Rules (Will Block Merge)

- [ ] **NO hardcoded subject IDs anywhere** - No patient IDs, usubjid, subject identifiers, or any specific subject codes in:
  - Production code
  - Debug statements
  - Example code
  - Comments or documentation
  - Test fixtures (use synthetic IDs instead)
  - Use dynamic selection: `slice_sample(n=1)`, `filter()`, or parameterize the ID

- [ ] **NO install.packages() calls** - Dependencies must be managed via `renv`:
  - Remove all `install.packages()` and `if (!require(pkg)) install.packages(pkg)`
  - Add new packages to `renv.lock` using `renv::snapshot()`

- [ ] **Use native pipe `|>` not `%>%`** - Modern R pipe operator required throughout

- [ ] **NO syntax errors** - Code must parse without errors:
  - Run basic syntax check: `Rscript -e 'source("your_file.R")'`
  - For Stan: `stanc --include-paths=stan,stan/ssls your_model.stan`

### Code Quality Requirements

- [ ] **Follow tidyverse style guide** for R code

- [ ] **NO global assignments (`<<-`)** in functions - Return values instead

- [ ] **Proper function defaults** - Don't reference undefined variables in parameter defaults

- [ ] **Documentation matches implementation** - Roxygen docs must match actual return values and parameters

- [ ] **NO tar_config_set()** - Use explicit `store=` argument in all targets functions

### Recommended Pre-PR Workflow

```bash
# 1. Before committing - Quick validation
Rscript -e 'styler::style_file("r/your_modified_file.R")'  # Auto-format
grep -r "install.packages" r/  # Check for install calls
grep -r "%>%" r/  # Check for old pipe operator

# 2. Search for hardcoded subject IDs (common patterns)
grep -rE "(E[0-9]{7,}|[A-Z][0-9]{4}[0-9]{3}[0-9]{3})" r/ --include="*.R"

# 3. Test your changes
Rscript -e 'source("r/your_file.R")'  # Basic syntax check
Rscript -e 'testthat::test_dir("tests/testthat")'  # Run tests

# 4. Before creating PR - Automated review
# In Claude Code CLI:
/pr-review-toolkit:review-pr code
```

### Using Automated PR Review

```bash
# In Claude Code
/pr-review-toolkit:review-pr code  # Check code quality and compliance
/pr-review-toolkit:review-pr tests  # Verify test coverage
/pr-review-toolkit:review-pr all    # Comprehensive review (recommended)
```
