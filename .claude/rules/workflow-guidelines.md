## Bash and Command Execution

- **NEVER pipe long-running commands to `head`, `tail`, or similar** when running in background - it prevents real-time output monitoring
- If you need to capture output while preserving streaming, use `tee` instead: `command | tee output.log`
- Background tasks automatically capture output to a file - don't truncate it with pipes
- Example:
  ```bash
  # ❌ DON'T: User can't see real-time progress
  Rscript -e 'targets::tar_make()' | head -100

  # ✓ DO: Full streaming output visible
  Rscript -e 'targets::tar_make()'

  # ✓ ALTERNATIVE: If you need to save output too
  Rscript -e 'targets::tar_make()' | tee build.log
  ```

## Data Pipeline Safety

- **NEVER run data preparation scripts without explicit confirmation from the user**
- A hookify rule (`.claude/hookify.data-pipeline-confirmation.local.md`) blocks accidental execution of:
  - `r/data_preparation_pipeline/*/main_pipeline_data_*.R`
  - Scripts containing `wrangle_data`, `prepare.*data`, or `pipeline.*` patterns
- **Why this matters**: Data pipeline scripts regenerate CSV files and can take hours to run
- **Before running**: Check if processed CSV files already exist at the expected location
- **Better alternative**: Use `targets::tar_make()` with specific target names to rebuild only what changed
- **Pattern**: Data pipelines should only run when raw data is updated, not for routine analysis

## Git Worktree Workflow

When working with multiple git worktrees, follow this pattern to avoid duplicate commits:

1. **Make changes in your current worktree** - Don't `cd` to other worktrees to make the same changes
2. **Commit locally** - Commit your work in the worktree where you made the changes
3. **Merge from the target worktree** - Switch to the other worktree and merge or rebase

**Example:**
```bash
# In worktree A (/mnt/code/worktrees/code-feature-x): make changes and commit
git add .claude/settings.json
git commit -m "Update plugin configuration"

# In worktree B (/mnt/code on main branch): merge the changes
cd /mnt/code
git merge feature-x
git push origin main
```

**Common mistake:**
```bash
# ❌ DON'T: Making the same change in multiple worktrees separately
cd /mnt/code && edit file && git commit
cd /mnt/code/worktrees/code-feature-x && edit file && git commit  # Duplicate!

# ✓ DO: Make change once, then merge
cd /mnt/code/worktrees/code-feature-x && edit file && git commit
cd /mnt/code && git merge feature-x
```

## General Coding Rules

- **No backward-compatibility aliases**: Do not create variable or function aliases for backward compatibility unless explicitly requested. When renaming, update all references directly instead of adding shims or aliases.
