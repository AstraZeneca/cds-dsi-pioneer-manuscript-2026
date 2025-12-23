# GitHub Copilot Instructions for This Project

## Code Style Preferences

### Prefer Tidyverse Functions
- Use `readr::read_lines()` instead of `readLines()`
- Use `stringr::str_starts()` instead of `startsWith()`
- Use `stringr::str_glue()` instead of `sprintf()` or `paste0()`
- Use `stringr::str_c()` instead of `paste()`
- Use `purrr::map()` family instead of `lapply()`, `sapply()`, etc.
- Use pipes (`|>`) for sequential operations
- Use modern lambda syntax `\(x)` for anonymous functions

## Critical File Editing Safety Rules

### ⚠️ NEVER Truncate Files
When editing files, especially large R or Python files with multiple functions:

1. **Always read the ENTIRE file first** before making edits
2. **Verify the end of the file** is included in your context
3. **Check line counts** - if a file appears to end abruptly, read more
4. **Use `read_file` without limits** to see the complete file structure
5. **Never assume** a file ends where your initial read stops

### When Making Bulk Formatting Changes

If performing operations like:
- Reformatting function calls (multi-line → single-line)
- Changing argument patterns across multiple functions
- Updating code style/aesthetics
- Any "search and replace" style operations

**REQUIRED STEPS:**
1. Read the file from start to END using `read_file` (check the last line number)
2. Note the total number of functions/sections in the file
3. After editing, verify ALL functions are still present
4. Check the line count hasn't decreased unexpectedly
5. Use `grep_search` to verify critical functions still exist

### Example Safety Check

```bash
# Before editing: Count functions
grep -c "^[a-zA-Z_].*<- function" r/plot_functions.R

# After editing: Verify count matches
grep -c "^[a-zA-Z_].*<- function" r/plot_functions.R
```

### For R Function Files Specifically

Files like `r/plot_functions.R`, `r/table_functions.R`, `r/util.R` often contain:
- Many sequential function definitions
- Functions at the END of the file that are easy to lose
- Important utilities that aren't imported elsewhere

**Before committing changes to these files:**
- Scan for function names at the end of the file
- Verify the file ends with `# nolint end` or similar expected marker
- Check git diff line counts: large deletions (-700 lines) are suspicious

## Recovery Procedure

If functions are accidentally deleted:
1. Use `git log -p -S "function_name"` to find when it was deleted
2. Use `git show <commit>^:path/to/file.R` to see the file before deletion
3. Extract and restore the missing functions
4. Verify all related files still work

## File-Specific Warnings

### `r/plot_functions.R`
- Contains many plotting utilities
- Functions like `plot_lfo_elpd_diff` are at the end
- Always verify the last function after bulk edits

### `r/sclc/plot_functions.R`
- Project-specific plotting functions
- Check both this AND parent `plot_functions.R`

## Testing After Major Edits

After reformatting or bulk changes:
1. Source the file in R to check for syntax errors
2. Run `grep "^[a-zA-Z_].*function" file.R | wc -l` before/after
3. Check for any functions called in `.qmd` files that might now be missing
