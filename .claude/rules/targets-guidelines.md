---
paths:
  - "**/targets/*.R"
  - "_targets.yaml"
---

## Targets Coding Guidelines

- **NEVER use `tar_config_set(store = ...)`** - it changes global state and causes conflicts
- Always use explicit `store` argument: `tar_read(name, store = "path/_targets")`
- Same applies to all targets functions: `tar_meta()`, `tar_load()`, etc.
- **IMPORTANT**: `TAR_RUN` environment variable does NOT work with `tar_make()` - always use explicit `store="/path/_targets"` argument
- `_targets.yaml` sclc store uses `!expr` with `DOMINO_STARTING_USERNAME` and `TAR_RUN` — never hardcode username or branch in this file
- `tumor_ssls_draws_pop` selection: `time_invariant_coef_qr_*` and `time_varying_coef_*` params don't follow the `_pop` suffix — they need `matches("^(time_invariant|time_varying)_coef")` added to the `select_draws` call
- **NEVER inline complex code in targets** - extract to helper functions in `r/` directory
  - Target commands should be simple function calls, not multi-line code blocks
  - Example: Use `tar_target(name, my_function(arg))` not `tar_target(name, { ... complex code ... })`
  - Helper functions belong in appropriate `r/` subdirectories (e.g., `r/sclc/plot_functions.R`)
- **`pattern = map()` dependencies**: When adding analysis-data-dependent post-processing to a mapped target (e.g., `cutoff_tumor_ssls_stan_data`), add the analysis data target to the `map()` pattern as well.
- **Track `source()` files with `format = "file"`**: If a target calls `source("path/to/file.R")` inside its expression, targets does NOT detect changes to that file. Add a separate `tar_target(my_script, "path/to/file.R", format = "file")` and reference `my_script` in the `source()` call. See `initializers_fixed_file` and `prepare_ms_standalone_data_script` for the established pattern.
- **`tar_invalidate` and `tar_map` naming**: Invalidation target names must include the full `tar_map` suffix. `tar_invalidate(base_foo)` is a no-op if the actual target is `base_foo_jan26` (inside `tar_map(dco_name)`). Always use the `tar-map-names` skill to get the correct full name before passing it to `-i`.
