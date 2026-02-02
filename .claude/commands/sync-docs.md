# Documentation Sync Skill

## Purpose

This skill helps keep documentation synchronized with source code changes by:
1. Identifying which source files have changed
2. Detecting new files that aren't in any mapping (and suggesting where to add them)
3. Finding documentation that relates to those files (via `.claude/doc-mappings.yaml`)
4. Reviewing both source and doc for discrepancies
5. Suggesting updates to keep docs current

## Usage

```bash
# Check recent git changes against main branch
/sync-docs

# Check only staged changes
/sync-docs --staged

# Check specific file(s)
/sync-docs --file path/to/file.R
/sync-docs --file stan/ssls/sf-ssls-lfo.stan

# Full audit of all mappings (comprehensive review)
/sync-docs --all

# Check uncommitted changes only
/sync-docs --uncommitted
```

## How It Works

### 1. Load Mappings

Read `.claude/doc-mappings.yaml` to understand which source files map to which documentation files.

### 2. Identify Changed Files

Depending on the invocation mode:

- **Default** (`/sync-docs`): Run `git diff main...HEAD --name-only` to find files changed on current branch
- **--staged**: Run `git diff --staged --name-only` to find staged files
- **--uncommitted**: Run `git diff HEAD --name-only` to find uncommitted changes
- **--file**: Use the explicitly provided file path
- **--all**: Skip this step and review ALL mappings

### 3. Match Against Mappings

For each changed file, check if it matches any `sources` pattern in the mappings using glob pattern matching.

### 4. Check for Unmapped New Files

Before reviewing source-doc pairs, identify NEW files that don't match any existing mapping patterns:

a. **Detect new files**: Use `git diff --diff-filter=A` to find newly added files
b. **Test against all mapping patterns**: Check if each new file matches any `sources` pattern
c. **Flag unmapped files**: Identify files that should probably be tracked but aren't in any mapping

For each unmapped new file, determine:
- **File type and location**: Is it a Stan module? R utility? Documentation?
- **Related mappings**: Which existing mapping(s) might this file relate to?
- **Mapping suggestions**:
  - If it clearly fits an existing mapping → suggest adding the pattern
  - If it's a new area → suggest creating a new mapping entry

### 5. Review Source-Doc Pairs

For each matched mapping:

a. **Read the relevant source files** (all files matching the patterns for that mapping)
b. **Read the corresponding documentation files**
c. **Compare and identify gaps**, such as:
   - New parameters in Stan not mentioned in docs
   - Changed function signatures not reflected in docs
   - New files/modules not described
   - Outdated examples or code snippets
   - Mathematical formulations that changed
   - Prior specifications that changed

### 6. Generate Report

Output a structured report:

```markdown
## Documentation Sync Report

### Summary
- X mappings affected by changes
- Y documentation files need review
- Z new files not in any mapping

---

## ⚠️  Unmapped New Files

The following new files don't match any existing mapping patterns:

### 1. `path/to/new_file.R`

**File type:** R utility function
**Likely related to:** "Parameter Inference Analysis" or "Model Validation Analysis"

**Recommendation:** Add to `.claude/doc-mappings.yaml`

```yaml
# Option 1: Add to existing mapping
- name: "Parameter Inference Analysis"
  sources:
    - "r/posterior.R"
    - "r/new_file.R"  # <-- ADD THIS
```

Or:

```yaml
# Option 2: Create new mapping if this is a new functional area
- name: "New Feature Area"
  description: "What this new file/area does"
  sources:
    - "r/new_file.R"
  docs:
    - "quarto/website/analysis/new-page.qmd"
```

---

## 📋 Documentation Updates Needed

### 1. [Mapping Name]

**Changed source files:**
- path/to/file1.stan
- path/to/file2.R

**Related documentation:**
- quarto/website/documentation/example.qmd

**Findings:**
- [ ] New parameter `tr_foo_bar` added in priors.R but not documented
- [ ] Function signature changed: `compute_accuracy()` now takes `cutoff` argument
- [ ] Module structure changed: new file `ar1_helpers.stan` not mentioned

**Suggested actions:**
1. Add documentation for `tr_foo_bar` in Section 2.3
2. Update function reference in validation page
3. Add ar1_helpers.stan to architecture diagram

---

### 2. [Next Mapping]
...
```

### 7. Interactive Follow-up

After generating the report, ask the user:

```
Would you like me to:
1. Update specific documentation sections?
2. Create a detailed update plan?
3. Just keep this as a reference?
```

## Mapping File Format

The `.claude/doc-mappings.yaml` file structure:

```yaml
mappings:
  - name: "Mapping Name"
    description: "What this mapping covers"
    sources:
      - "path/pattern/**/*.ext"  # Glob patterns
      - "specific/file.R"
    docs:
      - "quarto/website/documentation/page.qmd"
```

### Pattern Matching Rules

- `**` matches any directory recursively
- `*` matches any filename
- Exact paths match exactly
- All patterns are evaluated from project root

## Implementation Notes

### Performance Optimization

- **Don't read everything**: Only load source/doc files for affected mappings
- **Use focused context**: Read specific sections if docs are large
- **Batch related changes**: Group related findings together

### Smart Diffing

Look for these common sync issues:

1. **New parameters**: Stan parameters not in docs
2. **Changed priors**: Prior specifications updated but docs still show old defaults
3. **Function signatures**: R function arguments changed
4. **Module structure**: New Stan module files not documented
5. **Mathematical formulas**: Model equations changed
6. **Code examples**: Snippets in docs don't match current code

### Unmapped File Detection

When detecting new files that don't match any mapping:

**Consider these factors:**
- File type (`.R`, `.stan`, `.qmd`, etc.)
- Location in directory structure (`r/`, `stan/ssls/modules/`, `targets/`, etc.)
- File purpose (utility, analysis, plotting, data prep, model specification)
- Related existing files (is there a similar file in a mapping?)

**Provide smart suggestions:**
- If it clearly extends an existing area → suggest adding to that mapping
- If it's a new functional area → suggest creating a new mapping
- If it's a one-off script or temp file → note it can be ignored

**Files to skip:**
- Test files (`tests/**`)
- Build artifacts (`_targets/`, `_freeze/`, `_site/`)
- Temporary files (`temp.R`, `scratch.stan`)
- Configuration files (`.gitignore`, `.Rprofile`, etc.)
- Data files (`.csv`, `.parquet`, etc.)

### False Positive Reduction

Some changes don't require doc updates:
- Refactoring without behavior change
- Internal helper functions
- Comment changes
- Formatting/style changes

Use judgment to distinguish meaningful changes from trivial ones.

## Examples

### Example 1: Parameter Added

**Change detected:**
```r
# r/priors.R
tr_sd_process_noise = 0.1  # NEW LINE
```

**Finding:**
```
- [ ] New prior hyperparameter `tr_sd_process_noise` added but not documented in model-specification.qmd Section 2.4
```

**Suggested action:**
```
Add to model-specification.qmd after line 250:
- $\sigma_{\text{process}}$: Process noise standard deviation (default: 0.1)
```

### Example 2: Function Signature Changed

**Change detected:**
```r
# r/accuracy.R
-compute_accuracy(pred, obs)
+compute_accuracy(pred, obs, cutoff = NULL)
```

**Finding:**
```
- [ ] Function signature changed in r/accuracy.R: new `cutoff` parameter not reflected in leave-future-out-cross-validation.qmd
```

**Suggested action:**
```
Update documentation to explain the optional cutoff parameter for early data analysis.
```

### Example 3: New Stan Module

**Change detected:**
```
New file: stan/ssls/modules/tr/ar1_helpers.stan
```

**Finding:**
```
- [ ] New helper file ar1_helpers.stan not mentioned in architecture.qmd module structure
```

**Suggested action:**
```
Add to architecture.qmd Section on module organization:
- `ar1_helpers.stan` - Helper functions for AR(1) process noise computation
```

### Example 4: Unmapped New File

**Change detected:**
```
New file: r/sensitivity_analysis.R
```

**Finding:**
```
⚠️  This new file doesn't match any existing mapping patterns!

File type: R utility function
Likely related to: Model validation or analysis workflow
```

**Suggested action:**
```yaml
# Option 1: Add to existing "Model Validation Analysis" mapping
- name: "Model Validation Analysis"
  sources:
    - "r/accuracy.R"
    - "r/sclc/accuracy.R"
    - "r/plot_lfo_diagram.R"
    - "r/sensitivity_analysis.R"  # <-- ADD THIS

# Option 2: Create new mapping if this starts a new analysis area
- name: "Sensitivity Analysis"
  description: "Prior and parameter sensitivity analysis tools"
  sources:
    - "r/sensitivity_analysis.R"
  docs:
    - "quarto/website/analysis/sensitivity-analysis.qmd"  # (create this)
```

## Configuration

The skill reads configuration from `.claude/doc-mappings.yaml`. To add new mappings:

1. Identify source files that should trigger doc updates
2. Identify the corresponding documentation files
3. Add a new mapping entry with appropriate glob patterns
4. Test with `/sync-docs --file <your-source-file>`

## Tips for Effective Use

1. **Run regularly**: Check `/sync-docs` before committing changes to catch missing doc updates
2. **Use --staged**: Quick check before commit: `/sync-docs --staged`
3. **File-specific checks**: When editing a file, run `/sync-docs --file path/to/file` to see immediate impact
4. **Full audits**: Periodically run `/sync-docs --all` for comprehensive review
5. **Keep mappings updated**: The skill will suggest adding new files to mappings - review and update `.claude/doc-mappings.yaml` as suggested
6. **New file alerts**: When you add new source files, run `/sync-docs` to see if they need to be added to a mapping

## Integration with Workflow

### Pre-commit Hook (optional)

Add to `.claude/settings.json`:

```json
{
  "hooks": {
    "pre-commit": "echo '⚠️  Don't forget to check /sync-docs before committing!'"
  }
}
```

### Pull Request Workflow

Before creating a PR:
1. Run `/sync-docs` to identify doc gaps
2. Update documentation as needed
3. Re-render the Quarto website
4. Include doc updates in the same PR

---

## Task Instructions for Claude

When this skill is invoked:

1. **Parse arguments** to determine mode (default, --staged, --file, --all, --uncommitted)

2. **Load mappings** from `.claude/doc-mappings.yaml`

3. **Find changed files** using appropriate git command based on mode

4. **Check for unmapped new files:**
   - Use `git diff --diff-filter=A --name-only` to find newly added files
   - For each new file, test if it matches any existing mapping's source patterns
   - If it doesn't match any pattern, flag it as unmapped
   - Analyze the file path/name to suggest which mapping it belongs to or if it needs a new mapping
   - Provide concrete YAML snippets to add it to `.claude/doc-mappings.yaml`
   - Skip files that clearly don't need tracking (test files, temp files, build artifacts, etc.)

5. **Match files to mappings** using glob pattern matching

6. **For each affected mapping:**
   - Read the changed source files (focus on the changes if possible)
   - Read the related documentation files
   - Identify discrepancies, missing content, outdated information
   - Generate specific, actionable findings

7. **Generate a structured report** as described above, with two sections:
   - First: "⚠️  Unmapped New Files" (if any)
   - Second: "📋 Documentation Updates Needed" (for matched mappings)

8. **Ask user for next steps**: update mappings, update docs, create plan, or just reference

9. **Be efficient**: Don't read unnecessary files; focus on affected mappings only

10. **Be specific**: Point to exact sections, line numbers, parameter names

11. **Be helpful**: Suggest concrete changes, not just "this needs updating"

12. **For unmapped files**: Provide ready-to-use YAML that can be directly added to the mappings file
