# targets MCP Plugin Expansion — Design Spec

**Date:** 2026-03-20
**Author:** Karim Naguib
**Plugin:** `targets-toolkit` (Pioneer Claude Marketplace)
**Current version:** 1.4.4

---

## Overview

Expand the `targets-toolkit` MCP plugin from 10 tools to 19 tools (after Phase 4), adding coverage of ~55 user-facing functions from the `{targets}` R package that are appropriate for interactive MCP use. Currently the plugin covers only basic inspection, metadata, and loading. This expansion adds full pipeline observability, branching analysis, configuration queries, debugging support, destructive store management, and pipeline execution.

---

## Background

### Current State (10 tools)

| Tool | Purpose |
|---|---|
| `tar_read` | Read a target's value |
| `tar_load` | Load targets into r_session |
| `tar_meta` | Read target metadata |
| `tar_progress` | Read target execution progress |
| `tar_outdated` | Check which targets will rerun |
| `tar_manifest` | List all targets with commands |
| `tar_invalidate` | Force targets to rerun (delete metadata) |
| `list_objects` | List objects in r_session |
| `eval_expr` | Evaluate R expression in r_session |
| `reset_session` | Clear r_session environment |

### Architecture

R-based MCP server using `callr::r_session` for isolation. All tool calls execute inside a pre-warmed, persistent R session. Auto-restarts on crash. Size limit of 1 MB for returned objects (larger objects return `str()` summary).

Server: `/home/ubuntu/.claude/plugins/marketplaces/pioneer-claude-marketplace/plugins/targets-toolkit/mcp/targets-mcp/server.R`

### Scope of Expansion

The `{targets}` package exports ~120 symbols. Of these:
- ~60 are appropriate for MCP (interactive inspection/management)
- ~40 are developer assertions (`tar_assert_*`) — excluded
- ~10 are runtime-only (work inside running targets) — excluded
- ~10 are RStudio addins, knitr integration, or Shiny apps — excluded
- ~10 are deprecated/superseded — excluded
- `tar_config_set()` — explicitly forbidden by CLAUDE.md

---

## Design

### Approach: Grouped Composite Tools

Related functions are grouped into composite tools with an `action` parameter. This keeps the total tool count manageable (~17) while providing full coverage. All tools use `action` as the discriminating parameter (not `mode`, `operation`, etc.) for consistency.

---

## Tool Specifications

### Unchanged Tools (6)

These tools are well-designed as-is and require no changes:

- `tar_read` — already uses `tar_read_raw` internally
- `tar_outdated` — single-purpose, frequently used
- `tar_manifest` — single-purpose
- `list_objects` — utility, not targets-specific
- `eval_expr` — utility escape hatch
- `reset_session` — utility

---

### Enhanced Existing Tools (3)

#### `tar_load` (enhanced)

Add loading modes to support `tar_load_everything()` and `tar_load_globals()`.

**Parameters:**

| Param | Type | Default | Description |
|---|---|---|---|
| `store` | string | required | Path to `_targets/` store (ignored for `action: "globals"`) |
| `action` | enum | `"targets"` | `"targets"` \| `"everything"` \| `"globals"` |
| `names` | string | optional | Comma-separated target names (only for `action: "targets"`) |

**Actions:**
- `targets` — current behavior, calls `tar_load(names, store)`
- `everything` — calls `tar_load_everything(store)`, loads all available targets
- `globals` — calls `tar_load_globals()`, loads global objects for debugging. **Note:** `tar_load_globals()` does not accept a `store` argument — the `store` parameter is ignored for this action.

---

#### `tar_meta` (enhanced)

Add `fields` parameter to allow column selection.

**Parameters:**

| Param | Type | Default | Description |
|---|---|---|---|
| `store` | string | required | Path to `_targets/` store |
| `names` | string | optional | Comma-separated target names |
| `fields` | string | optional | Comma-separated column names to return (e.g. `"name,error,time"`) |

No actions added — `tar_meta_delete` moves to `tar_destroy`. Keep this tool read-only.

---

#### `tar_progress` (enhanced)

Add summary and branch-level views.

**Parameters:**

| Param | Type | Default | Description |
|---|---|---|---|
| `store` | string | required | Path to `_targets/` store |
| `action` | enum | `"full"` | `"full"` \| `"summary"` \| `"branches"` |
| `names` | string | optional | Target name filter (only for `action: "branches"`) |

**Actions:**
- `full` — current behavior, calls `tar_progress(store)`
- `summary` — calls `tar_progress_summary(store)`, returns aggregate counts per status
- `branches` — calls `tar_progress_branches(store, names)`, branch-level progress for dynamic targets

---

### New Tools (8)

#### `tar_status` (new)

Query which targets are in a given state. All actions return a character vector of target names, optionally enriched with metadata.

**Parameters:**

| Param | Type | Default | Description |
|---|---|---|---|
| `store` | string | required | Path to `_targets/` store |
| `action` | enum | required | See actions below |
| `time` | string | optional | Datetime string for `newer`/`older` (e.g. `"2026-01-01"`) |
| `pattern` | string | optional | Regex/glob pattern for `described_as` action |
| `include_metadata` | boolean | `false` | If true, join result with `tar_meta(fields = "name,error,time")` and return a data frame instead of a character vector |

**Actions:**

| Action | Function | Description |
|---|---|---|
| `completed` | `tar_completed(store)` | Targets that finished successfully |
| `errored` | `tar_errored(store)` | Targets that errored |
| `canceled` | `tar_canceled(store)` | Targets that were canceled |
| `skipped` | `tar_skipped(store)` | Targets skipped as up-to-date |
| `dispatched` | `tar_dispatched(store)` | Targets currently running |
| `objects` | `tar_objects(store)` | Targets with saved output in store |
| `newer` | `tar_newer(time, store)` | Targets completed after `time` |
| `older` | `tar_older(time, store)` | Targets completed before `time` |
| `described_as` | `tar_objects(names = described_as(pattern), store)` | Targets whose description matches `pattern`. Uses the `described_as()` tidyselect helper — there is no standalone `tar_described_as()` function. Implementation: `tar_objects(names = described_as(pattern), store = store)` |

**`include_metadata` join schema:** When `include_metadata: true`, the return value is a data frame with columns `name` (character), `error` (character or NA), `time` (POSIXct). Join key is `name`.

---

#### `tar_inspect` (new)

Detailed pipeline and target inspection. All actions are read-only.

**Parameters:**

| Param | Type | Default | Description |
|---|---|---|---|
| `store` | string | required | Path to `_targets/` store |
| `action` | enum | required | See actions below |
| `names` | string | optional | Single target name for per-target actions (`deps`, `timestamp`, `traceback`). Multiple comma-separated names supported only for `timestamp`. |
| `format` | enum | `"mermaid"` | Visualization format for `visualize` action: `"mermaid"` \| `"network"` |

**Actions:**

| Action | Function | Description |
|---|---|---|
| `sitrep` | `tar_sitrep(store)` | Cue-by-cue status table for all targets |
| `deps` | `tar_deps(name)` | Code dependencies of a target. **Important:** reads the pipeline script (`_targets.R`), not the store — requires the pipeline script to be present and sourceable. Takes a single `names` value. |
| `timestamp` | `tar_timestamp(names, store)` | Build timestamp(s) for one or more targets. `names` may be comma-separated for multiple targets; results are returned as a named vector. |
| `traceback` | `tar_traceback(names, store)` | Error traceback for a failed target. Takes a single `names` value. |
| `pid` | `tar_pid()` | Main process ID for the most recent pipeline run |
| `process` | `tar_process(store)` | Process metadata (PID, start time, etc.) |
| `visualize` | `tar_mermaid()` or `tar_network()` | Dependency graph. `format: "mermaid"` calls `tar_mermaid(store)` and returns a mermaid.js text string. `format: "network"` calls `tar_network(store)` and returns a named list with two data frames: `$vertices` (columns: `name`, `type`, `label`, `status`) and `$edges` (columns: `from`, `to`). |

---

#### `tar_branches` (new)

Dynamic branching information for targets created with `pattern = map(...)` or similar.

**Parameters:**

| Param | Type | Default | Description |
|---|---|---|---|
| `store` | string | required | Path to `_targets/` store |
| `action` | enum | required | See actions below |
| `name` | string | required | Parent (stem) target name |

**Actions:**

| Action | Function | Description |
|---|---|---|
| `index` | `tar_branch_index(name, store)` | Integer indexes of each branch for a dynamic target |
| `names` | `tar_branch_names(name, store)` | Character names of each branch (e.g. `"my_target_abc123"`) |
| `list` | `tar_branches(name, store)` | Data frame with branch names and dependency branch names for full reconstruction |

---

#### `tar_config` (new)

Configuration and path queries. All read-only.

**Parameters:**

| Param | Type | Default | Description |
|---|---|---|---|
| `store` | string | optional | Path to `_targets/` store (where applicable) |
| `action` | enum | required | See actions below |
| `name` | string | optional | Setting name for `get` and `option_get`; target name for `paths` with `path_type: "target"` |
| `path_type` | enum | optional | `"store"` \| `"script"` \| `"target"` — required when `action: "paths"` |

**Actions:**

| Action | Function | Description |
|---|---|---|
| `get` | `tar_config_get(name)` | Read a named setting from `_targets.yaml` (e.g. `"store"`, `"script"`). Returns all settings if `name` is omitted. |
| `projects` | `tar_config_projects()` | List all project names configured in `_targets.yaml` |
| `yaml` | `tar_config_yaml()` | Returns the parsed contents of `_targets.yaml` as a structured object. Implementation calls `tar_config_yaml()` (which returns the file path) then reads and parses that file — the tool returns contents, not the path. |
| `envvars` | `tar_envvars()` | Show all targets-related environment variables and their current values |
| `option_get` | `tar_option_get(name)` | Get the current value of a global target option (e.g. `"format"`, `"memory"`, `"garbage_collection"`) |
| `paths` | varies | Path resolution: `path_type: "store"` → `tar_path_store()`; `"script"` → `tar_path_script()`; `"target"` → `tar_path_target(name)` |

---

#### `tar_validate` (new)

Validate the pipeline before running. Standalone because it is commonly used as a pre-`tar_make()` check.

**Parameters:**

| Param | Type | Default | Description |
|---|---|---|---|
| `store` | string | required | Path to `_targets/` store |
| `script` | string | optional | Path to the pipeline script (defaults to `_targets.R` in working directory). Useful for multi-project setups where `_targets.yaml` defines non-default script paths. |

Returns: validation result — either a success message or a formatted list of errors found.

---

#### `tar_debug` (new)

Debugging and workspace tools.

**Parameters:**

| Param | Type | Default | Description |
|---|---|---|---|
| `store` | string | optional | Path to `_targets/` store (where applicable) |
| `action` | enum | required | See actions below |
| `name` | string | optional | Target name for `workspace_load` and `instructions` |
| `script` | string | optional | Script path for `source` action |
| `seed` | integer | optional | Seed override for `seed_create` |

**Actions:**

| Action | Function | Description |
|---|---|---|
| `workspace_load` | `tar_workspace(name, store)` | Load a saved workspace + seed for a failed target into the r_session for interactive debugging. **WARNING:** This loads many objects into the shared r_session environment, which will contaminate subsequent tool calls. Always call `reset_session` after finishing debug work. |
| `workspace_list` | `tar_workspaces(store)` | List all saved target workspaces (targets must have been defined with `workspace_on_error = TRUE`) |
| `instructions` | `tar_debug_instructions(name)` | Print step-by-step instructions for debugging a specific failed target |
| `source` | `tar_source(script)` | Source one or more R helper scripts (e.g. `tar_helper` support files) into the r_session |
| `seed_create` | `tar_seed_create(name, seed)` | Create a deterministic reproducible seed for a named target |

---

#### `tar_destroy` (new)

**DESTRUCTIVE operations.** Replaces and extends the existing `tar_invalidate` tool. All actions that modify the store require `confirm: true`.

**Parameters:**

| Param | Type | Default | Description |
|---|---|---|---|
| `store` | string | required | Path to `_targets/` store |
| `action` | enum | required | See actions below |
| `names` | string | optional | Comma-separated target names (supports tidyselect patterns). Not applicable to `destroy` or `unblock`. |
| `scope` | string | optional | Scope for `destroy` action only. Passed as `destroy` arg to `tar_destroy()`. Valid values: `"all"` (default), `"cloud"`, `"local"`, `"meta"`, `"process"`, `"progress"`, `"objects"`, `"scratch"`, `"workspaces"`. |
| `confirm` | boolean | `false` | Must be `true` for all destructive actions. Non-destructive `prune_list` does not require it. |

**Actions:**

| Action | Function | Destructive | Description |
|---|---|---|---|
| `delete` | `tar_delete(names, store)` | Yes | Delete output values (cached objects) for named targets |
| `invalidate` | `tar_invalidate(names, store)` | Yes | Delete metadata records for targets, forcing rerun. Cached values in `_targets/objects/` are preserved. |
| `meta_delete` | `tar_meta_delete(names, store)` | Yes | Delete metadata records entirely |
| `prune` | `tar_prune(store)` | Yes | Remove targets and their output that are no longer part of the pipeline |
| `prune_list` | `tar_prune_list(store)` | No | Dry run — list which targets `prune` would remove. No `confirm` required. |
| `destroy` | `tar_destroy(destroy, store)` | Yes | Delete store contents. Scope controlled by `scope` parameter (default `"all"` wipes entire store). |
| `unblock` | `tar_unblock_process(store)` | Yes | Delete the process lock file (`_targets/meta/process`). Use only when a pipeline has crashed and the store is stuck in a locked state. **Danger:** unblocking a genuinely running pipeline will corrupt its state. |

**Safety implementation:**
- Server validates `confirm == TRUE` before executing any destructive action; returns an informative error otherwise
- Store path validated before any write: directory must exist and contain `_targets/meta/`
- `destroy` action: if `names` is non-empty, return an error (names are not meaningful for full-store destroy)
- `unblock` action note: deletes `_targets/meta/process`, not a data file — document that this is appropriate only for crash recovery, not for stopping a running pipeline

**Migration from `tar_invalidate`:**
- The existing `tar_invalidate` tool is kept in Phase 3 as a deprecated alias
- Its functionality is available as `tar_destroy(action: "invalidate", ...)`
- Removed in Phase 4

---

#### `tar_make` (new)

Execute the pipeline. Long-running — returns immediately with monitoring instructions.

**Parameters:**

| Param | Type | Default | Description |
|---|---|---|---|
| `store` | string | required | Path to `_targets/` store |
| `names` | string | optional | Comma-separated target names to build (subset pipeline; builds only named targets and their dependencies) |
| `confirm` | boolean | `false` | Must be `true` to start execution |

**Behavior:**
1. Validates `confirm: true`
2. Detects if a pipeline is already running for this store using `tar_exist_process(store)` combined with checking whether the PID from `tar_process(store)$pid` is an active OS process. Refuses to start if already running.
3. Starts `tar_make()` in a `callr::r_bg()` background process (non-blocking). The background process's stdout/stderr are routed to a log file at `<store>/../tar_make.log` for post-hoc inspection.
4. Returns immediately: `"Pipeline started. Background process PID: <pid>. Log: <path>. Use tar_progress(store=...) to monitor progress."`

**Note:** `tar_make_clustermq()` and `tar_make_future()` are superseded by `crew`-based parallelism in modern targets — not exposed.

---

## Tool Count by Phase

| Phase | Version | Tools Added | Tools Removed | Running Total |
|---|---|---|---|---|
| Baseline | 1.4.4 | — | — | 10 |
| Phase 1 | 1.5.0 | tar_status, tar_inspect, tar_progress (enhanced) | — | 12 |
| Phase 2 | 1.6.0 | tar_branches, tar_config, tar_debug, tar_validate, tar_load (enhanced), tar_meta (enhanced) | — | 18 |
| Phase 3 | 2.0.0 | tar_destroy | tar_invalidate (deprecated, kept as alias) | 19 |
| Phase 4 | 2.1.0 | tar_make | tar_invalidate (alias removed) | **19** |

*Note: The net change is +9 tools (19 - 10). The existing `tar_invalidate` tool is superseded by `tar_destroy(action: "invalidate")` and removed after one deprecation cycle.*

---

## Phasing Plan

### Phase 1 — High-Value Inspection (v1.5.0)

**Scope:** New read-only tools covering pipeline observability.

| Change | Risk |
|---|---|
| Add `tar_status` | Low |
| Add `tar_inspect` | Low |
| Enhance `tar_progress` (summary + branches actions) | Low |

**Functions covered:** ~20 new functions
**Rationale:** All read-only, zero risk. Addresses the most common gaps: "what failed?", "why is this outdated?", "show me the dependency graph."

---

### Phase 2 — Branching, Config & Debug (v1.6.0)

**Scope:** Complete the read-only toolset.

| Change | Risk |
|---|---|
| Add `tar_branches` | Low |
| Add `tar_config` | Low |
| Add `tar_debug` | Low |
| Add `tar_validate` | Low |
| Enhance `tar_load` (everything + globals actions) | Low |
| Enhance `tar_meta` (fields parameter) | Low |

**Functions covered:** ~15 new functions
**Rationale:** Read-only (except tar_load which only loads into r_session). Covers the "long tail" of inspection and debugging needs. Note: `tar_debug` workspace_load action pollutes the shared r_session — documented in the tool description, mitigated by reset_session.

---

### Phase 3 — Destructive Operations (v2.0.0)

**Scope:** Store cleanup and invalidation.

| Change | Risk |
|---|---|
| Add `tar_destroy` | Medium (breaking) |
| Deprecate `tar_invalidate` (keep as alias for one cycle) | Breaking |

**Functions covered:** ~7 new functions
**Rationale:** Destructive operations need careful safety gating. `tar_invalidate` alias deprecation is a breaking change — major version bump. One version cycle of deprecation before removal.

---

### Phase 4 — Pipeline Execution (v2.1.0)

**Scope:** Run `tar_make()` via MCP.

| Change | Risk |
|---|---|
| Add `tar_make` | High (async, side effects) |
| Remove deprecated `tar_invalidate` alias | Breaking |

**Rationale:** Highest complexity (async `r_bg()` process management, PID tracking, log routing). Users already launch pipelines via `sclc_targets.sh` and Domino jobs — this is a convenience addition, not a blocker. Implement last to ensure solid foundation.

---

## Implementation Notes

### Parameter Naming Convention

- Use `action` as the discriminating parameter everywhere (not `mode`, `operation`, `type`)
- Use `names` for target name arguments (matches R function convention); use `name` (singular) only when exactly one target is accepted
- Comma-separated strings for multi-name inputs (matches existing plugin pattern)
- `store` is always the first parameter in tool definitions

### Tool Description Guidelines

- First sentence: what the tool does
- Second sentence: when to use it vs. related tools
- List all `action` values with one-line descriptions
- For destructive tools: include explicit "DESTRUCTIVE — requires confirm: true" warning at the top of the description

### Testing Per Phase

Each phase follows:
1. Implement tool handlers (`tool_*` functions) in `server.R`
2. Register tool definitions (name, description, parameters) in tool list
3. Manual smoke test: call each action via MCP, verify output format and schema
4. Edge case tests: invalid store path, invalid action, missing required params, destructive actions without `confirm: true`
5. Update `SKILL.md` / `AGENT.md` with new capabilities
6. Bump version in both `plugin.json` AND `marketplace.json`

### Version Bump Requirement

Per CLAUDE.md: both `plugin.json` and `marketplace.json` must be bumped together. The plugin manager won't reinstall without a version bump in both files.

---

## Out of Scope

The following are explicitly excluded:

| Category | Examples | Reason |
|---|---|---|
| Pipeline definition | `tar_target`, `tar_cue`, `tar_format` | Used in `_targets.R`, not interactively |
| Developer assertions | `tar_assert_*` (40+ functions) | Internal validation helpers |
| Runtime-only | `tar_name`, `tar_active`, `tar_envir` | Only work inside a running target |
| Condition signaling | `tar_throw_*`, `tar_warn_*` | Package internals |
| RStudio addins | `rstudio_addin_*` | IDE-specific |
| Knitr integration | `tar_engine_knitr`, `tar_interactive` | Rmd/Qmd authoring |
| Cloud operations | `tar_meta_upload`, `tar_meta_sync` | Not relevant to local workflow |
| Shiny apps | `tar_watch`, `tar_watch_ui` | Not suitable for MCP request/response |
| Blocking loops | `tar_poll` | Blocking console loop incompatible with MCP |
| Config write | `tar_config_set` | CLAUDE.md explicitly forbids usage |
| Config write | `tar_config_unset` | Removes settings from `_targets.yaml` — too risky without a full config management story; excluded until a config-write tool is designed separately |
| Deprecated | `tar_make_future`, `tar_bind`, `tar_pipeline` | Obsolete |
| Visualization widgets | `tar_visnetwork`, `tar_glimpse` | Return HTML widgets; `tar_mermaid`/`tar_network` cover this |
