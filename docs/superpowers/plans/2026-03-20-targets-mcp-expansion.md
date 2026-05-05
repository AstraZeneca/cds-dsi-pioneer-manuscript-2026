# targets MCP Plugin Expansion Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Expand the targets-toolkit MCP plugin from 10 to 19 tools, covering ~55 {targets} R package functions organized as grouped composite tools with `action`-based dispatch.

**Architecture:** All tools follow a 3-part pattern in `server.R`: a `tool_*` function that runs inside a `callr::r_session`, a JSON Schema entry in the `TOOLS` list, and a dispatch case in `dispatch_tool()`. No test harness exists for the MCP server — testing is (1) verify R function in isolation via `Rscript -e`, then (2) smoke-test via MCP tool call.

**Tech Stack:** R, callr, jsonlite, targets package. No new dependencies. (`yaml` package is used in `tar_config` but is already available as a transitive dependency of targets itself.)

**Spec:** `docs/superpowers/specs/2026-03-20-targets-mcp-expansion-design.md`

---

## Files

| File | Change |
|---|---|
| `server.R` | All tool additions/modifications — add to TOOLS list, tool_* functions, dispatch switch |
| `plugin.json` | Version bump after each phase |
| `marketplace.json` | Version bump (must match plugin.json) |
| `SKILL.md` | Document new tool capabilities |
| `AGENT.md` | Update targets-explorer agent capabilities |

All file paths:
- **Server:** `/home/ubuntu/.claude/plugins/marketplaces/pioneer-claude-marketplace/plugins/targets-toolkit/mcp/targets-mcp/server.R`
- **plugin.json:** `/home/ubuntu/.claude/plugins/marketplaces/pioneer-claude-marketplace/plugins/targets-toolkit/.claude-plugin/plugin.json`
- **marketplace.json:** `/home/ubuntu/.claude/plugins/marketplaces/pioneer-claude-marketplace/.claude-plugin/marketplace.json`
- **SKILL.md:** `/home/ubuntu/.claude/plugins/marketplaces/pioneer-claude-marketplace/plugins/targets-toolkit/skills/tar-map-names/SKILL.md` (or nearest plugin-level SKILL.md)
- **AGENT.md:** `/home/ubuntu/.claude/plugins/marketplaces/pioneer-claude-marketplace/plugins/targets-toolkit/agents/targets-explorer/AGENT.md`

**Test store for smoke tests:** `/mnt/data/analysis-results/karim_naguib/sclc/dco3/_targets`

---

## server.R Structure Reference

The 3 places to edit for every new tool (line numbers shift as you add — find by function name):
1. **TOOLS list** — `list(name = "...", description = "...", inputSchema = list(...))`
2. **tool_* function** — `tool_tar_xxx <- function(...) { ... }` — runs inside r_session
3. **dispatch switch** — `tar_xxx = .rs$run(tool_tar_xxx, args = list(...)),`

**Exception — `tar_make`:** Uses a `server_tool_tar_make()` function that runs directly in the **server process** (not via `.rs$run()`). This is because `callr::r_bg()` must be launched from the server process, not from inside the `callr::r_session` worker — nested callr session spawning can deadlock due to inherited file descriptors and pipe handles. In `dispatch_tool()`, `tar_make` gets a special case that calls `server_tool_tar_make(args)` directly without `.rs$run()`.

---

## Phase 1 — High-Value Inspection (v1.5.0)

### Task 1: Enhance `tar_progress` — add summary and branches actions

**Files:**
- Modify: `server.R` (TOOLS entry + tool_tar_progress + dispatch)

- [ ] **Step 1: Verify R functions in isolation**

```bash
Rscript -e '
  library(targets)
  store <- "/mnt/data/analysis-results/karim_naguib/sclc/dco3/_targets"
  print(tar_progress_summary(store = store))
  print(head(tar_progress_branches(store = store)))
'
```
Expected: data frames printed without error. Note column names.

- [ ] **Step 2: Replace the `tool_tar_progress` function in server.R**

Find the existing `tool_tar_progress` function and replace its body to handle the new `action` and `names` parameters:

```r
tool_tar_progress <- function(store, action = "full", names = NULL) {
  action <- match.arg(action, c("full", "summary", "branches"))
  switch(action,
    full = {
      progress <- targets::tar_progress(store = store)
      list(success = TRUE, data = as.list(progress), row_count = nrow(progress))
    },
    summary = {
      s <- targets::tar_progress_summary(store = store)
      list(success = TRUE, data = as.list(s))
    },
    branches = {
      name_vec <- if (!is.null(names)) trimws(strsplit(names, ",")[[1]]) else NULL
      pb <- targets::tar_progress_branches(names = name_vec, store = store)
      list(success = TRUE, data = as.list(pb), row_count = nrow(pb))
    }
  )
}
```

- [ ] **Step 3: Update the TOOLS entry for `tar_progress`**

Find the `tar_progress` entry in the TOOLS list and replace the `inputSchema` to add `action` and `names`:

```r
list(
  name = "tar_progress",
  description = paste(
    "Read pipeline execution progress. Shows which targets are running, completed, errored, or skipped.",
    "Actions: 'full' (default) - full per-target progress table;",
    "'summary' - aggregate counts per status (how many completed/errored/etc);",
    "'branches' - branch-level progress for dynamic targets (use names to filter by parent target)."
  ),
  inputSchema = list(
    type = "object",
    properties = list(
      store = list(type = "string", description = "Path to _targets/ store directory"),
      action = list(
        type = "string",
        enum = list("full", "summary", "branches"),
        description = "Progress view: 'full' (default), 'summary', or 'branches'"
      ),
      names = list(
        type = "string",
        description = "Comma-separated parent target names to filter branches (only for action='branches')"
      )
    ),
    required = list("store")
  )
)
```

- [ ] **Step 4: Update the dispatch case for `tar_progress`**

Find and replace the `tar_progress` dispatch line:

```r
tar_progress = .rs$run(tool_tar_progress,
  args = list(store = args$store, action = args$action %||% "full", names = args$names)),
```

- [ ] **Step 5: Smoke test**

```bash
# In a Claude Code session with the plugin loaded, run these MCP tool calls:
# (or restart the MCP server and call via eval_expr in the MCP session)

# Test 1: summary action
# Call tar_progress with store=<store>, action="summary"
# Expected: JSON with aggregate counts (skipped, completed, errored, etc.)

# Test 2: branches action
# Call tar_progress with store=<store>, action="branches"
# Expected: branch-level data frame or empty if no dynamic targets
```

- [ ] **Step 6: Commit**

```bash
git -C /home/ubuntu/.claude/plugins/marketplaces/pioneer-claude-marketplace \
  add plugins/targets-toolkit/mcp/targets-mcp/server.R
git -C /home/ubuntu/.claude/plugins/marketplaces/pioneer-claude-marketplace \
  commit -m "feat(targets-mcp): enhance tar_progress with summary and branches actions"
```

---

### Task 2: Add `tar_status` tool

**Files:**
- Modify: `server.R`

- [ ] **Step 1: Verify R functions in isolation**

```bash
Rscript -e '
  library(targets)
  store <- "/mnt/data/analysis-results/karim_naguib/sclc/dco3/_targets"
  print(tar_completed(store = store))
  print(tar_errored(store = store))
  print(tar_objects(store = store))
  # Test newer/older with a date
  print(tar_newer(Sys.time() - 86400 * 30, store = store))  # last 30 days
'
```
Expected: character vectors printed without error.

- [ ] **Step 2: Add `tool_tar_status` function to server.R**

Add after the last `tool_*` function and before `dispatch_tool`:

```r
tool_tar_status <- function(store, action, time = NULL, pattern = NULL,
                             include_metadata = FALSE) {
  action <- match.arg(action, c("completed", "errored", "canceled", "skipped",
                                 "dispatched", "objects", "newer", "older",
                                 "described_as"))

  names_vec <- switch(action,
    completed   = targets::tar_completed(store = store),
    errored     = targets::tar_errored(store = store),
    canceled    = targets::tar_canceled(store = store),
    skipped     = targets::tar_skipped(store = store),
    dispatched  = targets::tar_dispatched(store = store),
    objects     = targets::tar_objects(store = store),
    newer = {
      if (is.null(time)) stop("time parameter required for 'newer' action")
      targets::tar_newer(time = as.POSIXct(time), store = store)
    },
    older = {
      if (is.null(time)) stop("time parameter required for 'older' action")
      targets::tar_older(time = as.POSIXct(time), store = store)
    },
    described_as = {
      if (is.null(pattern)) stop("pattern parameter required for 'described_as' action")
      targets::tar_objects(names = targets::described_as(pattern), store = store)
    }
  )

  if (isTRUE(include_metadata) && length(names_vec) > 0) {
    meta <- targets::tar_meta(
      names = names_vec,
      fields = c("name", "error", "time"),
      store = store
    )
    list(success = TRUE, data = as.list(meta), row_count = nrow(meta))
  } else {
    list(success = TRUE, names = names_vec, count = length(names_vec))
  }
}
```

- [ ] **Step 3: Add `tar_status` entry to TOOLS list**

Add after the `tar_progress` TOOLS entry:

```r
list(
  name = "tar_status",
  description = paste(
    "Query which targets are in a given state. Returns target names, optionally enriched with metadata.",
    "Actions: 'completed' - finished successfully; 'errored' - failed with an error;",
    "'canceled' - were canceled; 'skipped' - up-to-date and skipped; 'dispatched' - currently running;",
    "'objects' - have saved output in store; 'newer' - completed after a given time (requires time param);",
    "'older' - completed before a given time (requires time param);",
    "'described_as' - match a description pattern (requires pattern param).",
    "Set include_metadata=true to also return name/error/time columns in one call."
  ),
  inputSchema = list(
    type = "object",
    properties = list(
      store = list(type = "string", description = "Path to _targets/ store directory"),
      action = list(
        type = "string",
        enum = list("completed", "errored", "canceled", "skipped", "dispatched",
                    "objects", "newer", "older", "described_as"),
        description = "Which state to query"
      ),
      time = list(
        type = "string",
        description = "Datetime string for newer/older actions (e.g. '2026-01-01' or '2026-01-01 12:00:00')"
      ),
      pattern = list(
        type = "string",
        description = "Description pattern for described_as action (regex)"
      ),
      include_metadata = list(
        type = "boolean",
        description = "If true, join with tar_meta to include name/error/time columns"
      )
    ),
    required = list("store", "action")
  )
)
```

- [ ] **Step 4: Add dispatch case**

Inside `dispatch_tool()`, add after the `tar_progress` dispatch line:

```r
tar_status = .rs$run(tool_tar_status,
  args = list(
    store = args$store,
    action = args$action,
    time = args$time,
    pattern = args$pattern,
    include_metadata = isTRUE(args$include_metadata)
  )),
```

- [ ] **Step 5: Smoke test**

Call `tar_status` with:
- `store=<store>, action="errored"` — should return list of errored target names
- `store=<store>, action="completed", include_metadata=true` — should return data frame with name/error/time
- `store=<store>, action="newer", time="2026-01-01"` — should return recent targets
- `action="newer"` without `time` — should return `success=false` with helpful error

- [ ] **Step 6: Commit**

```bash
git -C /home/ubuntu/.claude/plugins/marketplaces/pioneer-claude-marketplace \
  add plugins/targets-toolkit/mcp/targets-mcp/server.R
git -C /home/ubuntu/.claude/plugins/marketplaces/pioneer-claude-marketplace \
  commit -m "feat(targets-mcp): add tar_status tool for pipeline state queries"
```

---

### Task 3: Add `tar_inspect` tool

**Files:**
- Modify: `server.R`

- [ ] **Step 1: Verify R functions in isolation**

```bash
Rscript -e '
  library(targets)
  store <- "/mnt/data/analysis-results/karim_naguib/sclc/dco3/_targets"
  # Check sitrep
  print(head(tar_sitrep(store = store)))
  # Check network structure
  net <- tar_network(store = store)
  cat("network vertices:", ncol(net$vertices), "cols,", nrow(net$vertices), "rows\n")
  cat("network edges:", ncol(net$edges), "cols,", nrow(net$edges), "rows\n")
  cat("vertex cols:", paste(names(net$vertices), collapse=", "), "\n")
  cat("edge cols:", paste(names(net$edges), collapse=", "), "\n")
  # Check mermaid output
  m <- tar_mermaid(store = store)
  cat("mermaid lines:", length(m), "\n")
  cat(head(m, 5), sep="\n")
  # Check process
  print(tar_process(store = store))
  # IMPORTANT: verify tar_deps signature with a string argument
  # Pick any target name from the manifest
  first_target <- tar_manifest(store = store)$name[1]
  cat("Testing tar_deps with target:", first_target, "\n")
  d <- tar_deps(first_target)
  print(d)
'
```
Expected: output without errors. Note exact column names of sitrep, network$vertices, network$edges. Verify `tar_deps()` accepts a plain string (not NSE) in the installed version.

- [ ] **Step 2: Add `tool_tar_inspect` function to server.R**

```r
tool_tar_inspect <- function(store, action, names = NULL, format = "mermaid") {
  action <- match.arg(action, c("sitrep", "deps", "timestamp", "traceback",
                                 "pid", "process", "visualize"))
  switch(action,
    sitrep = {
      s <- targets::tar_sitrep(store = store)
      list(success = TRUE, data = as.list(s), row_count = nrow(s))
    },
    deps = {
      if (is.null(names)) stop("names required for deps action (single target name)")
      d <- targets::tar_deps(names)
      list(success = TRUE, deps = d)
    },
    timestamp = {
      name_vec <- if (!is.null(names)) trimws(strsplit(names, ",")[[1]]) else NULL
      if (is.null(name_vec)) {
        ts <- targets::tar_timestamp(store = store)
        list(success = TRUE, timestamps = list(format(ts)))
      } else {
        ts <- sapply(name_vec, function(n) {
          format(targets::tar_timestamp(name = n, store = store))
        })
        list(success = TRUE, timestamps = as.list(ts))
      }
    },
    traceback = {
      if (is.null(names)) stop("names required for traceback action (single target name)")
      tb <- targets::tar_traceback(name = names, store = store)
      list(success = TRUE, traceback = tb)
    },
    pid = {
      list(success = TRUE, pid = targets::tar_pid())
    },
    process = {
      p <- targets::tar_process(store = store)
      list(success = TRUE, data = as.list(p))
    },
    visualize = {
      format <- match.arg(format, c("mermaid", "network"))
      if (format == "mermaid") {
        m <- targets::tar_mermaid(store = store)
        list(success = TRUE, mermaid = paste(m, collapse = "\n"))
      } else {
        net <- targets::tar_network(store = store)
        list(
          success = TRUE,
          vertices = as.list(net$vertices),
          edges = as.list(net$edges)
        )
      }
    }
  )
}
```

- [ ] **Step 3: Add `tar_inspect` entry to TOOLS list**

```r
list(
  name = "tar_inspect",
  description = paste(
    "Detailed pipeline and target inspection. All actions are read-only.",
    "Actions: 'sitrep' - cue-by-cue status table showing exactly why each target is or isn't outdated;",
    "'deps' - code dependencies of a target (reads pipeline script, not store; requires names=single target);",
    "'timestamp' - build timestamp(s) for target(s) (names=comma-separated, optional);",
    "'traceback' - full error traceback for a failed target (requires names=single target);",
    "'pid' - PID of the most recent pipeline process;",
    "'process' - process metadata (PID, start time);",
    "'visualize' - dependency graph (format='mermaid' returns mermaid.js text;",
    "format='network' returns vertices/edges data frames)."
  ),
  inputSchema = list(
    type = "object",
    properties = list(
      store = list(type = "string", description = "Path to _targets/ store directory"),
      action = list(
        type = "string",
        enum = list("sitrep", "deps", "timestamp", "traceback", "pid", "process", "visualize"),
        description = "Inspection action to perform"
      ),
      names = list(
        type = "string",
        description = "Target name(s). Single name for deps/traceback. Comma-separated for timestamp."
      ),
      format = list(
        type = "string",
        enum = list("mermaid", "network"),
        description = "Visualization format for action='visualize': 'mermaid' (default) or 'network'"
      )
    ),
    required = list("store", "action")
  )
)
```

- [ ] **Step 4: Add dispatch case**

```r
tar_inspect = .rs$run(tool_tar_inspect,
  args = list(
    store = args$store,
    action = args$action,
    names = args$names,
    format = args$format %||% "mermaid"
  )),
```

- [ ] **Step 5: Smoke test**

Call `tar_inspect` with:
- `action="sitrep"` — should return table with cue columns
- `action="visualize", format="mermaid"` — should return mermaid text starting with `graph`
- `action="visualize", format="network"` — should return object with `vertices` and `edges`
- `action="pid"` — should return a PID integer
- `action="deps"` without `names` — should return `success=false` with error

- [ ] **Step 6: Commit**

```bash
git -C /home/ubuntu/.claude/plugins/marketplaces/pioneer-claude-marketplace \
  add plugins/targets-toolkit/mcp/targets-mcp/server.R
git -C /home/ubuntu/.claude/plugins/marketplaces/pioneer-claude-marketplace \
  commit -m "feat(targets-mcp): add tar_inspect tool for pipeline inspection"
```

---

### Task 4: Phase 1 version bump and docs (v1.5.0)

**Files:**
- Modify: `plugin.json`
- Modify: `marketplace.json`
- Modify: `AGENT.md`

- [ ] **Step 1: Read current versions**

```bash
cat /home/ubuntu/.claude/plugins/marketplaces/pioneer-claude-marketplace/plugins/targets-toolkit/.claude-plugin/plugin.json | grep version
grep targets-toolkit /home/ubuntu/.claude/plugins/marketplaces/pioneer-claude-marketplace/.claude-plugin/marketplace.json
```

- [ ] **Step 2: Bump version in plugin.json**

Change `"version": "1.4.4"` to `"version": "1.5.0"` in plugin.json.

- [ ] **Step 3: Bump version in marketplace.json**

Find the targets-toolkit entry and change its version to `"1.5.0"` to match.

- [ ] **Step 4: Update AGENT.md with new tool capabilities**

In the targets-explorer agent's AGENT.md, add a section listing the new tools:
- `tar_status` — query which targets are in a given state
- `tar_inspect` — sitrep, deps, timestamps, tracebacks, dependency graph
- `tar_progress` — now also supports `summary` and `branches` actions

- [ ] **Step 5: Commit**

```bash
cd /home/ubuntu/.claude/plugins/marketplaces/pioneer-claude-marketplace
git add plugins/targets-toolkit/.claude-plugin/plugin.json \
        .claude-plugin/marketplace.json \
        plugins/targets-toolkit/agents/
git commit -m "feat(targets-toolkit): release v1.5.0 — tar_status and tar_inspect tools"
```

---

## Phase 2 — Branching, Config & Debug (v1.6.0)

### Task 5: Enhance `tar_load` — add everything and globals actions

**Files:**
- Modify: `server.R`

- [ ] **Step 1: Verify R functions in isolation**

```bash
Rscript -e '
  library(targets)
  store <- "/mnt/data/analysis-results/karim_naguib/sclc/dco3/_targets"
  # tar_load_everything loads ALL objects - only test on small store
  # tar_load_globals does not need a store
  tar_load_globals()
  cat("globals loaded\n")
  ls()
'
```

- [ ] **Step 2: Update `tool_tar_load` function in server.R**

Find the existing `tool_tar_load` function and replace:

```r
tool_tar_load <- function(store, action = "targets", names = NULL) {
  action <- match.arg(action, c("targets", "everything", "globals"))
  switch(action,
    targets = {
      if (is.null(names)) stop("names required for action='targets'")
      name_vec <- trimws(strsplit(names, ",")[[1]])
      targets::tar_load(names = name_vec, store = store, envir = .GlobalEnv)
      list(success = TRUE, message = paste("Loaded:", paste(name_vec, collapse = ", ")))
    },
    everything = {
      targets::tar_load_everything(store = store, envir = .GlobalEnv)
      loaded <- ls(.GlobalEnv)
      list(success = TRUE, message = paste("Loaded all targets. Objects in session:", length(loaded)),
           names = loaded)
    },
    globals = {
      # tar_load_globals() does not accept a store argument
      targets::tar_load_globals()
      loaded <- ls(.GlobalEnv)
      list(success = TRUE, message = paste("Loaded globals. Objects in session:", length(loaded)),
           names = loaded)
    }
  )
}
```

- [ ] **Step 3: Update the TOOLS entry for `tar_load`**

Replace existing `tar_load` inputSchema to add `action`:

```r
list(
  name = "tar_load",
  description = paste(
    "Load target values into the r_session environment for use with eval_expr.",
    "Actions: 'targets' (default) - load specific named targets (requires names);",
    "'everything' - load all available targets from store;",
    "'globals' - load global objects (functions, data) from _targets.R script (no store needed).",
    "After loading, use list_objects to see what was loaded, and eval_expr to inspect values."
  ),
  inputSchema = list(
    type = "object",
    properties = list(
      store = list(type = "string",
                   description = "Path to _targets/ store (ignored for action='globals')"),
      action = list(
        type = "string",
        enum = list("targets", "everything", "globals"),
        description = "Loading mode: 'targets' (default, requires names), 'everything', or 'globals'"
      ),
      names = list(
        type = "string",
        description = "Comma-separated target names to load (required for action='targets')"
      )
    ),
    required = list("store")
  )
)
```

- [ ] **Step 4: Update dispatch case**

```r
tar_load = .rs$run(tool_tar_load,
  args = list(
    store = args$store,
    action = args$action %||% "targets",
    names = args$names
  )),
```

- [ ] **Step 5: Smoke test**

- `action="globals"` — should load functions into session; list_objects should show them
- `action="everything"` — should load all targets; check list_objects
- `action="targets"` without `names` — should return helpful error

- [ ] **Step 6: Commit**

```bash
git -C /home/ubuntu/.claude/plugins/marketplaces/pioneer-claude-marketplace \
  add plugins/targets-toolkit/mcp/targets-mcp/server.R
git -C /home/ubuntu/.claude/plugins/marketplaces/pioneer-claude-marketplace \
  commit -m "feat(targets-mcp): enhance tar_load with everything and globals actions"
```

---

### Task 6: Enhance `tar_meta` — add fields parameter

**Files:**
- Modify: `server.R`

- [ ] **Step 1: Verify fields parameter behavior**

```bash
Rscript -e '
  library(targets)
  store <- "/mnt/data/analysis-results/karim_naguib/sclc/dco3/_targets"
  # Check available fields
  all_meta <- tar_meta(store = store)
  cat("Available fields:", paste(names(all_meta), collapse=", "), "\n")
  # Check fields selection
  selected <- tar_meta(fields = c("name", "error", "time"), store = store)
  print(head(selected))
'
```

- [ ] **Step 2: Update `tool_tar_meta` function**

Find the existing `tool_tar_meta` function and replace:

```r
tool_tar_meta <- function(store, names = NULL, fields = NULL) {
  name_vec <- if (!is.null(names)) trimws(strsplit(names, ",")[[1]]) else NULL
  field_vec <- if (!is.null(fields)) trimws(strsplit(fields, ",")[[1]]) else NULL

  meta <- targets::tar_meta(
    names = name_vec,
    fields = field_vec,
    store = store
  )
  list(success = TRUE, data = as.list(meta), row_count = nrow(meta))
}
```

- [ ] **Step 3: Update the TOOLS entry for `tar_meta`**

Replace existing `tar_meta` inputSchema to add `fields`:

```r
list(
  name = "tar_meta",
  description = paste(
    "Read target metadata: name, type, format, bytes, seconds, error messages, warnings, timestamps.",
    "Use fields to select specific columns (e.g. fields='name,error,time' to focus on failures).",
    "Use names to filter to specific targets. Check the 'error' column to debug failures."
  ),
  inputSchema = list(
    type = "object",
    properties = list(
      store = list(type = "string", description = "Path to _targets/ store directory"),
      names = list(type = "string",
                   description = "Comma-separated target names to filter (optional, returns all if omitted)"),
      fields = list(type = "string",
                    description = paste("Comma-separated column names to return.",
                                        "Available: name, type, data, command, depend, seed, path,",
                                        "time, size, bytes, format, repository, iteration, parent,",
                                        "children, seconds, warnings, error, traceback, rep_workers"))
    ),
    required = list("store")
  )
)
```

- [ ] **Step 4: Update dispatch case**

```r
tar_meta = .rs$run(tool_tar_meta,
  args = list(store = args$store, names = args$names, fields = args$fields)),
```

- [ ] **Step 5: Smoke test**

- `fields="name,error,time"` — should return only 3 columns
- `names="some_target", fields="name,seconds,bytes"` — filtered subset
- No fields — should return all columns (existing behavior preserved)

- [ ] **Step 6: Commit**

```bash
git -C /home/ubuntu/.claude/plugins/marketplaces/pioneer-claude-marketplace \
  add plugins/targets-toolkit/mcp/targets-mcp/server.R
git -C /home/ubuntu/.claude/plugins/marketplaces/pioneer-claude-marketplace \
  commit -m "feat(targets-mcp): enhance tar_meta with fields parameter"
```

---

### Task 7: Add `tar_branches` tool

**Files:**
- Modify: `server.R`

- [ ] **Step 1: Verify R functions**

```bash
Rscript -e '
  library(targets)
  store <- "/mnt/data/analysis-results/karim_naguib/sclc/dco3/_targets"
  # Find a dynamic target to test with
  meta <- tar_meta(store = store)
  dynamic <- meta[!is.na(meta$parent), "name"]
  if (length(dynamic) > 0) {
    parent_name <- meta[!is.na(meta$parent), "parent"][1]
    cat("Testing with parent target:", parent_name, "\n")
    print(tar_branch_index(name = parent_name, store = store))
    print(head(tar_branch_names(name = parent_name, store = store)))
    print(head(tar_branches(name = parent_name, store = store)))
  } else {
    cat("No dynamic targets found in store\n")
  }
'
```

- [ ] **Step 2: Add `tool_tar_branches` function**

```r
tool_tar_branches <- function(store, action, name) {
  action <- match.arg(action, c("index", "names", "list"))
  switch(action,
    index = {
      idx <- targets::tar_branch_index(name = name, store = store)
      list(success = TRUE, indexes = idx, count = length(idx))
    },
    names = {
      branch_names <- targets::tar_branch_names(name = name, store = store)
      list(success = TRUE, names = branch_names, count = length(branch_names))
    },
    list = {
      branches_df <- targets::tar_branches(name = name, store = store)
      list(success = TRUE, data = as.list(branches_df), row_count = nrow(branches_df))
    }
  )
}
```

- [ ] **Step 3: Add `tar_branches` TOOLS entry**

```r
list(
  name = "tar_branches",
  description = paste(
    "Dynamic branching information for targets created with pattern=map(...) or similar.",
    "Requires a parent (stem) target name.",
    "Actions: 'index' - integer indexes of each branch;",
    "'names' - character names of each branch (e.g. 'my_target_abc123ef');",
    "'list' - data frame with branch names and their dependency branch names."
  ),
  inputSchema = list(
    type = "object",
    properties = list(
      store = list(type = "string", description = "Path to _targets/ store directory"),
      action = list(
        type = "string",
        enum = list("index", "names", "list"),
        description = "Branch query type"
      ),
      name = list(type = "string", description = "Parent (stem) target name")
    ),
    required = list("store", "action", "name")
  )
)
```

- [ ] **Step 4: Add dispatch case**

```r
tar_branches = .rs$run(tool_tar_branches,
  args = list(store = args$store, action = args$action, name = args$name)),
```

- [ ] **Step 5: Smoke test with a known dynamic target from the store**

- `action="names", name=<known_dynamic_target>` — returns branch name vector
- `action="list", name=<known_dynamic_target>` — returns data frame
- Missing required `name` — returns helpful error

- [ ] **Step 6: Commit**

```bash
git -C /home/ubuntu/.claude/plugins/marketplaces/pioneer-claude-marketplace \
  add plugins/targets-toolkit/mcp/targets-mcp/server.R
git -C /home/ubuntu/.claude/plugins/marketplaces/pioneer-claude-marketplace \
  commit -m "feat(targets-mcp): add tar_branches tool for dynamic branching info"
```

---

### Task 8: Add `tar_config` tool

**Files:**
- Modify: `server.R`

- [ ] **Step 1: Verify R functions**

```bash
Rscript -e '
  library(targets)
  library(yaml)
  # Test config functions (no store needed)
  print(tar_config_get())             # all settings
  print(tar_config_projects())         # project names
  yaml_path <- tar_config_yaml()
  cat("YAML path:", yaml_path, "\n")
  if (file.exists(yaml_path)) {
    print(yaml::read_yaml(yaml_path))
  }
  print(tar_envvars())
  print(tar_option_get("format"))
  # Path helpers
  cat("store path:", tar_path_store(), "\n")
  cat("script path:", tar_path_script(), "\n")
'
```
Note: Run from `/mnt/code` so `_targets.yaml` is found.

- [ ] **Step 2: Add `tool_tar_config` function**

```r
tool_tar_config <- function(store = NULL, action, name = NULL, path_type = NULL) {
  action <- match.arg(action, c("get", "projects", "yaml", "envvars",
                                 "option_get", "paths"))
  switch(action,
    get = {
      val <- targets::tar_config_get(name = name)
      list(success = TRUE, value = val)
    },
    projects = {
      projs <- targets::tar_config_projects()
      list(success = TRUE, projects = projs)
    },
    yaml = {
      yaml_path <- targets::tar_config_yaml()
      if (!file.exists(yaml_path)) {
        list(success = FALSE, error = paste("_targets.yaml not found at:", yaml_path))
      } else {
        contents <- yaml::read_yaml(yaml_path)
        list(success = TRUE, path = yaml_path, contents = contents)
      }
    },
    envvars = {
      ev <- targets::tar_envvars()
      list(success = TRUE, data = as.list(ev))
    },
    option_get = {
      val <- targets::tar_option_get(name = name)
      list(success = TRUE, option = name, value = val)
    },
    paths = {
      if (is.null(path_type)) stop("path_type required for paths action: 'store', 'script', or 'target'")
      path_type <- match.arg(path_type, c("store", "script", "target"))
      path_val <- switch(path_type,
        store  = targets::tar_path_store(),
        script = targets::tar_path_script(),
        target = {
          if (is.null(name)) stop("name required for path_type='target'")
          targets::tar_path_target(name = name)
        }
      )
      list(success = TRUE, path_type = path_type, path = path_val)
    }
  )
}
```

- [ ] **Step 3: Add `tar_config` TOOLS entry**

```r
list(
  name = "tar_config",
  description = paste(
    "Read configuration settings and path information. All actions are read-only.",
    "Actions: 'get' - read a named _targets.yaml setting (name optional, returns all if omitted);",
    "'projects' - list all configured project names;",
    "'yaml' - read and return the full _targets.yaml contents as parsed YAML;",
    "'envvars' - show all targets-related environment variables and their values;",
    "'option_get' - get a global target option value (e.g. name='format', 'memory', 'garbage_collection');",
    "'paths' - resolve a path (path_type='store'|'script'|'target'; name required for 'target')."
  ),
  inputSchema = list(
    type = "object",
    properties = list(
      store = list(type = "string",
                   description = "Path to _targets/ store (optional, used where relevant)"),
      action = list(
        type = "string",
        enum = list("get", "projects", "yaml", "envvars", "option_get", "paths"),
        description = "Configuration query to perform"
      ),
      name = list(
        type = "string",
        description = "Setting name for 'get'/'option_get'; target name for paths with path_type='target'"
      ),
      path_type = list(
        type = "string",
        enum = list("store", "script", "target"),
        description = "Path type for action='paths'"
      )
    ),
    required = list("action")
  )
)
```

- [ ] **Step 4: Add dispatch case**

```r
tar_config = .rs$run(tool_tar_config,
  args = list(
    store = args$store,
    action = args$action,
    name = args$name,
    path_type = args$path_type
  )),
```

- [ ] **Step 5: Smoke test**

- `action="projects"` — should list project names from `_targets.yaml`
- `action="yaml"` — should return parsed YAML contents
- `action="envvars"` — should return env var table
- `action="paths", path_type="store"` — should return store path string
- `action="paths"` without `path_type` — should return helpful error

- [ ] **Step 6: Commit**

```bash
git -C /home/ubuntu/.claude/plugins/marketplaces/pioneer-claude-marketplace \
  add plugins/targets-toolkit/mcp/targets-mcp/server.R
git -C /home/ubuntu/.claude/plugins/marketplaces/pioneer-claude-marketplace \
  commit -m "feat(targets-mcp): add tar_config tool for configuration queries"
```

---

### Task 9: Add `tar_validate` tool

**Files:**
- Modify: `server.R`

- [ ] **Step 1: Verify R function**

```bash
Rscript -e '
  library(targets)
  setwd("/mnt/code")  # so _targets.R is found
  tar_validate()
  cat("Validation passed\n")
' 2>&1
```
Expected: either silent success or validation errors printed.

- [ ] **Step 2: Add `tool_tar_validate` function**

```r
tool_tar_validate <- function(store, script = NULL) {
  result <- tryCatch({
    if (!is.null(script)) {
      targets::tar_validate(script = script, store = store)
    } else {
      targets::tar_validate(store = store)
    }
    list(success = TRUE, message = "Pipeline validation passed. No issues found.")
  }, error = function(e) {
    list(success = FALSE, error = conditionMessage(e),
         message = "Validation failed. Fix the errors above before running tar_make.")
  })
  result
}
```

- [ ] **Step 3: Add `tar_validate` TOOLS entry**

```r
list(
  name = "tar_validate",
  description = paste(
    "Validate the pipeline before running. Checks target definitions for correctness.",
    "Run before tar_make to catch configuration errors early.",
    "Returns success message if valid, or error details if validation fails.",
    "Use script parameter for non-default pipeline script paths (multi-project setups)."
  ),
  inputSchema = list(
    type = "object",
    properties = list(
      store = list(type = "string", description = "Path to _targets/ store directory"),
      script = list(
        type = "string",
        description = "Path to pipeline script (defaults to _targets.R in working directory)"
      )
    ),
    required = list("store")
  )
)
```

- [ ] **Step 4: Add dispatch case**

```r
tar_validate = .rs$run(tool_tar_validate,
  args = list(store = args$store, script = args$script)),
```

- [ ] **Step 5: Smoke test**

- Call with valid store — should return success message
- Call with a bogus script path — should return validation error
- Verify `success=false` cases have `isError=true` in MCP response

- [ ] **Step 6: Commit**

```bash
git -C /home/ubuntu/.claude/plugins/marketplaces/pioneer-claude-marketplace \
  add plugins/targets-toolkit/mcp/targets-mcp/server.R
git -C /home/ubuntu/.claude/plugins/marketplaces/pioneer-claude-marketplace \
  commit -m "feat(targets-mcp): add tar_validate tool"
```

---

### Task 10: Add `tar_debug` tool

**Files:**
- Modify: `server.R`

- [ ] **Step 1: Verify R functions**

```bash
Rscript -e '
  library(targets)
  store <- "/mnt/data/analysis-results/karim_naguib/sclc/dco3/_targets"
  # List workspaces (may be empty if no workspace_on_error targets)
  print(tar_workspaces(store = store))
  # Test seed creation
  seed <- tar_seed_create("some_target")
  cat("Seed for some_target:", seed, "\n")
  # Test tar_source with a real file
  # tar_source("r/util.R")  # only if that file exists
'
```

- [ ] **Step 2: Add `tool_tar_debug` function**

```r
tool_tar_debug <- function(store = NULL, action, name = NULL, script = NULL,
                            seed = NULL) {
  action <- match.arg(action, c("workspace_load", "workspace_list", "instructions",
                                 "source", "seed_create"))
  switch(action,
    workspace_load = {
      if (is.null(name)) stop("name required for workspace_load action")
      if (is.null(store)) stop("store required for workspace_load action")
      # WARNING: loads many objects into shared r_session — call reset_session after debugging
      # tar_workspace() uses tidyselect NSE internally; wrap string name with all_of()
      targets::tar_workspace(names = tidyselect::all_of(name), store = store)
      loaded <- ls(.GlobalEnv)
      list(
        success = TRUE,
        message = paste0(
          "Workspace loaded for target '", name, "'. ",
          length(loaded), " objects now in session. ",
          "WARNING: Run reset_session after debugging to clear contamination."
        ),
        objects = loaded
      )
    },
    workspace_list = {
      ws <- targets::tar_workspaces(store = store)
      list(success = TRUE, workspaces = ws, count = length(ws))
    },
    instructions = {
      if (is.null(name)) stop("name required for instructions action")
      out <- capture.output(targets::tar_debug_instructions(name = name))
      list(success = TRUE, instructions = paste(out, collapse = "\n"))
    },
    source = {
      if (is.null(script)) stop("script required for source action")
      targets::tar_source(files = script)
      list(success = TRUE, message = paste("Sourced:", script))
    },
    seed_create = {
      if (is.null(name)) stop("name required for seed_create action")
      s <- targets::tar_seed_create(name = name, seed = seed)
      list(success = TRUE, name = name, seed = s)
    }
  )
}
```

- [ ] **Step 3: Add `tar_debug` TOOLS entry**

```r
list(
  name = "tar_debug",
  description = paste(
    "Debugging and workspace tools for investigating failed targets.",
    "Actions: 'workspace_load' - load a saved workspace + seed for a failed target into the session",
    "(WARNING: contaminates r_session — run reset_session after; requires name + store);",
    "'workspace_list' - list all saved workspaces (requires workspace_on_error=TRUE in target def);",
    "'instructions' - print step-by-step debug instructions for a target (requires name);",
    "'source' - source a helper R script into the r_session (requires script path);",
    "'seed_create' - compute the deterministic RNG seed for a named target (requires name)."
  ),
  inputSchema = list(
    type = "object",
    properties = list(
      store = list(type = "string",
                   description = "Path to _targets/ store (required for workspace actions)"),
      action = list(
        type = "string",
        enum = list("workspace_load", "workspace_list", "instructions", "source", "seed_create"),
        description = "Debug action to perform"
      ),
      name = list(type = "string",
                  description = "Target name (required for workspace_load, instructions, seed_create)"),
      script = list(type = "string",
                    description = "Script file path for source action"),
      seed = list(type = "integer",
                  description = "Optional seed override for seed_create action")
    ),
    required = list("action")
  )
)
```

- [ ] **Step 4: Add dispatch case**

```r
tar_debug = .rs$run(tool_tar_debug,
  args = list(
    store = args$store,
    action = args$action,
    name = args$name,
    script = args$script,
    seed = args$seed
  )),
```

- [ ] **Step 5: Smoke test**

- `action="workspace_list"` with store — should return list (possibly empty)
- `action="seed_create", name="test_target"` — should return an integer
- `action="instructions", name="some_target"` — should return debug steps
- `action="workspace_load"` without `name` — should return helpful error

- [ ] **Step 6: Commit**

```bash
git -C /home/ubuntu/.claude/plugins/marketplaces/pioneer-claude-marketplace \
  add plugins/targets-toolkit/mcp/targets-mcp/server.R
git -C /home/ubuntu/.claude/plugins/marketplaces/pioneer-claude-marketplace \
  commit -m "feat(targets-mcp): add tar_debug tool"
```

---

### Task 11: Phase 2 version bump and docs (v1.6.0)

**Files:**
- Modify: `plugin.json`
- Modify: `marketplace.json`
- Modify: `AGENT.md`

- [ ] **Step 1: Bump both version files to `1.6.0`**

- [ ] **Step 2: Update AGENT.md**

Add to the agent capabilities section:
- `tar_branches` — branch indexes, names, and reconstruction for dynamic targets
- `tar_config` — `_targets.yaml` settings, env vars, path resolution
- `tar_validate` — pre-flight pipeline validation
- `tar_debug` — workspace loading, seed creation, sourcing helpers
- `tar_load` — now also supports `everything` and `globals` modes
- `tar_meta` — now supports `fields` parameter for column selection

- [ ] **Step 3: Commit**

```bash
cd /home/ubuntu/.claude/plugins/marketplaces/pioneer-claude-marketplace
git add plugins/targets-toolkit/.claude-plugin/plugin.json \
        .claude-plugin/marketplace.json \
        plugins/targets-toolkit/agents/
git commit -m "feat(targets-toolkit): release v1.6.0 — full read-only toolset complete"
```

---

## Phase 3 — Destructive Operations (v2.0.0)

### Task 12: Add `tar_destroy` tool and deprecate `tar_invalidate`

**Files:**
- Modify: `server.R`
- Modify: `plugin.json`
- Modify: `marketplace.json`

- [ ] **Step 1: Verify store validation and R functions**

```bash
Rscript -e '
  library(targets)
  store <- "/mnt/data/analysis-results/karim_naguib/sclc/dco3/_targets"
  # Verify prune_list (safe - read-only)
  print(tar_prune_list(store = store))
  # Verify store validation pattern
  cat("meta dir exists:", dir.exists(file.path(store, "meta")), "\n")
'
```

- [ ] **Step 2: Add helper function `validate_store_for_write`**

Add near the top of the tool implementations section (before `tool_tar_destroy`):

```r
# Helper: validate that a store path is a real targets store before any write operation
validate_store_for_write <- function(store) {
  if (!dir.exists(store)) {
    stop(paste("Store directory not found:", store))
  }
  if (!dir.exists(file.path(store, "meta"))) {
    stop(paste("Not a valid targets store (no meta/ directory):", store))
  }
  invisible(TRUE)
}
```

- [ ] **Step 3: Add `tool_tar_destroy` function**

```r
tool_tar_destroy <- function(store, action, names = NULL, scope = "all",
                              confirm = FALSE) {
  action <- match.arg(action, c("delete", "invalidate", "meta_delete",
                                 "prune", "prune_list", "destroy", "unblock"))

  # prune_list is a dry-run — no confirmation required
  if (action != "prune_list" && !isTRUE(confirm)) {
    stop(paste0(
      "DESTRUCTIVE action '", action, "' requires confirm=true. ",
      "This will permanently modify the targets store. ",
      "Set confirm=true to proceed."
    ))
  }

  # Validate store exists for all write operations
  if (action != "prune_list") {
    validate_store_for_write(store)
  }

  switch(action,
    delete = {
      if (is.null(names)) stop("names required for delete action")
      name_vec <- trimws(strsplit(names, ",")[[1]])
      targets::tar_delete(names = name_vec, store = store)
      list(success = TRUE, message = paste("Deleted output for:", paste(name_vec, collapse = ", ")))
    },
    invalidate = {
      if (is.null(names)) stop("names required for invalidate action")
      name_vec <- trimws(strsplit(names, ",")[[1]])
      targets::tar_invalidate(names = name_vec, store = store)
      list(success = TRUE, message = paste("Invalidated metadata for:", paste(name_vec, collapse = ", ")))
    },
    meta_delete = {
      if (is.null(names)) stop("names required for meta_delete action")
      name_vec <- trimws(strsplit(names, ",")[[1]])
      targets::tar_meta_delete(names = name_vec, store = store)
      list(success = TRUE, message = paste("Deleted metadata for:", paste(name_vec, collapse = ", ")))
    },
    prune = {
      targets::tar_prune(store = store)
      list(success = TRUE, message = "Pruned targets no longer in the pipeline.")
    },
    prune_list = {
      to_prune <- targets::tar_prune_list(store = store)
      list(success = TRUE, targets_to_prune = to_prune, count = length(to_prune),
           message = "Dry run — these targets would be removed by tar_destroy(action='prune').")
    },
    destroy = {
      if (!is.null(names) && nchar(trimws(names)) > 0) {
        stop("names is not applicable for the destroy action (destroys entire store scope). ",
             "Use action='delete' to remove specific targets.")
      }
      valid_scopes <- c("all", "cloud", "local", "meta", "process", "progress",
                        "objects", "scratch", "workspaces")
      scope <- match.arg(scope, valid_scopes)
      targets::tar_destroy(destroy = scope, store = store)
      list(success = TRUE,
           message = paste0("Store destroyed (scope='", scope, "'). ",
                            "The store at '", store, "' has been cleared."))
    },
    unblock = {
      targets::tar_unblock_process(store = store)
      list(success = TRUE,
           message = paste0(
             "Process lock file deleted. Store unblocked. ",
             "Only use this after a pipeline crash — ",
             "unblocking a running pipeline corrupts its state."
           ))
    }
  )
}
```

- [ ] **Step 4: Add `tar_destroy` TOOLS entry**

```r
list(
  name = "tar_destroy",
  description = paste(
    "DESTRUCTIVE — requires confirm=true. Store cleanup and invalidation operations.",
    "All actions except 'prune_list' permanently modify the store and require confirm=true.",
    "Actions: 'delete' - delete output values for named targets (requires names);",
    "'invalidate' - delete metadata for targets, forcing rerun (requires names);",
    "'meta_delete' - delete metadata records entirely (requires names);",
    "'prune' - remove targets no longer in the pipeline;",
    "'prune_list' - DRY RUN: list what prune would remove (safe, no confirm needed);",
    "'destroy' - delete store contents by scope (scope: 'all'=default/wipes store, 'meta', 'objects', etc.);",
    "'unblock' - delete process lock file after a pipeline crash (DANGER: corrupts a live pipeline)."
  ),
  inputSchema = list(
    type = "object",
    properties = list(
      store = list(type = "string", description = "Path to _targets/ store directory"),
      action = list(
        type = "string",
        enum = list("delete", "invalidate", "meta_delete", "prune", "prune_list", "destroy", "unblock"),
        description = "Destructive action to perform"
      ),
      names = list(
        type = "string",
        description = "Comma-separated target names (required for delete/invalidate/meta_delete)"
      ),
      scope = list(
        type = "string",
        enum = list("all", "cloud", "local", "meta", "process", "progress",
                    "objects", "scratch", "workspaces"),
        description = "Scope for destroy action (default 'all' wipes entire store)"
      ),
      confirm = list(
        type = "boolean",
        description = "Must be true for all destructive actions (all except prune_list)"
      )
    ),
    required = list("store", "action")
  )
)
```

- [ ] **Step 5: Add dispatch case for `tar_destroy`**

```r
tar_destroy = .rs$run(tool_tar_destroy,
  args = list(
    store = args$store,
    action = args$action,
    names = args$names,
    scope = args$scope %||% "all",
    confirm = isTRUE(args$confirm)
  )),
```

- [ ] **Step 6: Add deprecation warning to the existing `tar_invalidate` tool**

Find the `tar_invalidate` TOOLS entry and update its description to start with:
```
"[DEPRECATED] Use tar_destroy(action='invalidate') instead. Will be removed in v2.1.0. ..."
```

Also update `tool_tar_invalidate` to emit a warning in its return:
```r
tool_tar_invalidate <- function(store, names) {
  name_vec <- trimws(strsplit(names, ",")[[1]])
  targets::tar_invalidate(names = name_vec, store = store)
  list(
    success = TRUE,
    message = paste("Invalidated:", paste(name_vec, collapse = ", ")),
    deprecation_warning = "tar_invalidate is deprecated. Use tar_destroy(action='invalidate') instead. Will be removed in v2.1.0."
  )
}
```

- [ ] **Step 7: Smoke test**

- `action="prune_list"` — safe dry-run, returns list without confirm
- `action="delete"` without `confirm=true` — should return helpful error
- `action="destroy"` with `names` set — should return error about names not applicable
- `action="invalidate", names="some_target", confirm=true` — should succeed (test on a non-critical target)

- [ ] **Step 8: Bump versions to 2.0.0 in plugin.json and marketplace.json**

- [ ] **Step 9: Commit**

```bash
cd /home/ubuntu/.claude/plugins/marketplaces/pioneer-claude-marketplace
git add plugins/targets-toolkit/mcp/targets-mcp/server.R \
        plugins/targets-toolkit/.claude-plugin/plugin.json \
        .claude-plugin/marketplace.json
git commit -m "feat(targets-toolkit): v2.0.0 — add tar_destroy, deprecate tar_invalidate"
```

---

## Phase 4 — Pipeline Execution (v2.1.0)

### Task 13: Add `tar_make` tool and remove `tar_invalidate` alias

**Files:**
- Modify: `server.R`
- Modify: `plugin.json`
- Modify: `marketplace.json`

- [ ] **Step 1: Verify `callr::r_bg` pattern**

```bash
Rscript -e '
  library(callr)
  library(targets)
  store <- "/mnt/data/analysis-results/karim_naguib/sclc/dco3/_targets"

  # Test tar_exist_process for pipeline detection
  cat("Process exists:", tar_exist_process(store = store), "\n")

  # Test r_bg launch pattern (dry run — do not actually run tar_make)
  # Verify the log path convention
  log_path <- file.path(dirname(store), "tar_make.log")
  cat("Log would go to:", log_path, "\n")

  # Verify PID check pattern
  if (tar_exist_process(store = store)) {
    proc <- tar_process(store = store)
    pid <- proc$pid
    # Check if PID is alive (tools::SIGINT would not work; use system2)
    is_alive <- tryCatch({
      tools::pskill(pid, 0L)  # signal 0 = just check existence
      TRUE
    }, error = function(e) FALSE)
    cat("PID", pid, "is alive:", is_alive, "\n")
  }
'
```

- [ ] **Step 2: Add `server_tool_tar_make` function to server.R**

**Critical architecture note:** This function runs in the **server process** (not via `.rs$run()`), because `callr::r_bg()` must be launched from the server's own R process. Calling `callr::r_bg()` from inside a `callr::r_session` worker can deadlock due to inherited file descriptors. The pipeline-running check uses `.rs$run()` for the process check (since targets data is in the r_session), but the `r_bg()` launch happens in the server process.

Add `server_tool_tar_make` (NOT prefixed `tool_`) near the dispatch section, outside of any `tool_*` block:

```r
# Special case: runs in server process, not via .rs$run(), so callr::r_bg() is safe
server_tool_tar_make <- function(args) {
  store   <- args$store
  names   <- args$names
  confirm <- isTRUE(args$confirm)

  if (!confirm) {
    return(list(
      success = FALSE,
      error = paste0(
        "tar_make requires confirm=true. ",
        "This will execute the pipeline, which may run for a long time. ",
        "Monitor progress with tar_progress(store='", store, "')."
      )
    ))
  }

  # Check for already-running pipeline via r_session (where targets is loaded)
  proc_check <- tryCatch(
    .rs$run(function(store) {
      list(
        exists = targets::tar_exist_process(store = store),
        proc   = tryCatch(as.list(targets::tar_process(store = store)), error = function(e) NULL)
      )
    }, args = list(store = store)),
    error = function(e) list(exists = FALSE, proc = NULL)
  )

  if (isTRUE(proc_check$exists) && !is.null(proc_check$proc$pid)) {
    pid <- as.integer(proc_check$proc$pid)
    is_alive <- tryCatch({ tools::pskill(pid, 0L); TRUE }, error = function(e) FALSE)
    if (is_alive) {
      return(list(
        success = FALSE,
        error = paste0(
          "A pipeline is already running for this store (PID: ", pid, "). ",
          "Use tar_progress(store='", store, "') to monitor it, or ",
          "tar_destroy(action='unblock', confirm=true) only if the process has crashed."
        )
      ))
    }
  }

  name_vec <- if (!is.null(names) && nchar(trimws(names)) > 0) {
    trimws(strsplit(names, ",")[[1]])
  } else NULL

  log_path <- file.path(dirname(store), "tar_make.log")

  # Launch from server process — safe for callr::r_bg()
  bg <- callr::r_bg(
    func = function(name_vec, store, log_path) {
      library(targets)
      sink(log_path, append = FALSE, split = FALSE)
      on.exit(sink())
      if (is.null(name_vec)) {
        tar_make(store = store)
      } else {
        tar_make(names = name_vec, store = store)
      }
    },
    args = list(name_vec = name_vec, store = store, log_path = log_path),
    supervise = FALSE
  )

  list(
    success = TRUE,
    pid     = bg$get_pid(),
    log     = log_path,
    message = paste0(
      "Pipeline started (PID: ", bg$get_pid(), "). ",
      "Log: ", log_path, ". ",
      "Use tar_progress(store='", store, "') to monitor progress."
    )
  )
}
```

- [ ] **Step 3: Add `tar_make` TOOLS entry**

```r
list(
  name = "tar_make",
  description = paste(
    "Execute the pipeline. Requires confirm=true. Returns immediately — pipeline runs in background.",
    "Refuses to start if a pipeline is already running for this store.",
    "Monitor progress with tar_progress(store=...) after starting.",
    "Pipeline stdout/stderr are written to tar_make.log in the store's parent directory.",
    "Use names to build only specific targets and their dependencies (optional)."
  ),
  inputSchema = list(
    type = "object",
    properties = list(
      store = list(type = "string", description = "Path to _targets/ store directory"),
      names = list(
        type = "string",
        description = "Comma-separated target names to build (optional; builds all if omitted)"
      ),
      confirm = list(
        type = "boolean",
        description = "Must be true to start execution"
      )
    ),
    required = list("store", "confirm")
  )
)
```

- [ ] **Step 4: Add dispatch case for `tar_make` in `dispatch_tool()`**

This is the **only tool** that does NOT use `.rs$run()`. Add it as a special case before the main switch, or as a named case that calls `server_tool_tar_make` directly:

```r
# In dispatch_tool(), add BEFORE or alongside the switch():
if (tool_name == "tar_make") {
  return(server_tool_tar_make(args))
}
# (The rest of the switch handles all other tools via .rs$run())
```

- [ ] **Step 5: Remove the deprecated `tar_invalidate` tool**

  - Remove the `tar_invalidate` entry from the TOOLS list
  - Remove the `tool_tar_invalidate` function
  - Remove the `tar_invalidate` case from the dispatch switch

- [ ] **Step 6: Smoke test `tar_make`**

- `confirm=false` — should return helpful error with monitoring instructions
- `confirm=true` but pipeline already running — should return error with PID
- `confirm=true`, valid store, `names="some_fast_target"` — should start, return PID and log path
- After starting, call `tar_progress` to verify it's running
- Check log file exists and has output

- [ ] **Step 7: Smoke test `tar_invalidate` removal**

- Call `tar_invalidate` — should now return "Unknown tool: tar_invalidate" error
- Verify `tar_destroy(action="invalidate", ...)` still works as the replacement

- [ ] **Step 8: Bump versions to 2.1.0 in plugin.json and marketplace.json**

- [ ] **Step 9: Update AGENT.md and SKILL.md with final tool list**

Document the complete 19-tool inventory in both files.

- [ ] **Step 10: Commit**

```bash
cd /home/ubuntu/.claude/plugins/marketplaces/pioneer-claude-marketplace
git add plugins/targets-toolkit/mcp/targets-mcp/server.R \
        plugins/targets-toolkit/.claude-plugin/plugin.json \
        .claude-plugin/marketplace.json \
        plugins/targets-toolkit/agents/ \
        plugins/targets-toolkit/skills/
git commit -m "feat(targets-toolkit): v2.1.0 — add tar_make, remove deprecated tar_invalidate"
```

---

## Final Verification

After all phases are complete:

- [ ] **Confirm tool count**

```bash
grep '"name"' /home/ubuntu/.claude/plugins/marketplaces/pioneer-claude-marketplace/plugins/targets-toolkit/mcp/targets-mcp/server.R | grep -v "description\|properties\|store\|action\|names\|format\|confirm\|scope\|script\|seed\|time\|pattern\|path_type\|fields" | wc -l
```
Expected: 19

- [ ] **Confirm versions match**

```bash
grep version /home/ubuntu/.claude/plugins/marketplaces/pioneer-claude-marketplace/plugins/targets-toolkit/.claude-plugin/plugin.json
grep targets-toolkit /home/ubuntu/.claude/plugins/marketplaces/pioneer-claude-marketplace/.claude-plugin/marketplace.json | grep version
```
Expected: both show `2.1.0`

- [ ] **Restart MCP server and run full tool list**

After restarting the plugin, call `tools/list` and verify all 19 tools are present.
