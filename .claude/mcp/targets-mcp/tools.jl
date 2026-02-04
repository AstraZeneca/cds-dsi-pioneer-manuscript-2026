# MCP Tool Definitions for targets

# Server info for MCP
const SERVER_INFO = Dict(
    "name" => "targets-mcp",
    "version" => "0.1.0"
)

# Tool definitions following MCP schema
const TOOLS = [
    Dict(
        "name" => "tar_progress",
        "description" => "Read target progress from the most recent tar_make() run. Returns a data frame with columns: name, type (stem/pattern/branch), parent, branches count, and progress status. Progress values: 'dispatched' (sent to run, may be queued), 'completed', 'skipped', 'canceled', or 'errored'. Use to check pipeline status before loading targets.",
        "inputSchema" => Dict(
            "type" => "object",
            "properties" => Dict(
                "store" => Dict(
                    "type" => "string",
                    "description" => "Path to the _targets store directory. Required - do not rely on defaults."
                ),
                "project" => Dict(
                    "type" => "string",
                    "description" => "Project directory path (defaults to working directory)"
                )
            ),
            "required" => ["store"]
        )
    ),
    Dict(
        "name" => "tar_meta",
        "description" => "Read metadata for all recorded targets. Returns columns: name, type (function/object/stem/branch/map/cross), data hash, command hash, depend hash, seed, path, time (last modified), size hash, bytes, format, error, warnings, seconds (runtime). Use to debug failures (check 'error' column), understand performance (check 'seconds', 'bytes'), or verify target state.",
        "inputSchema" => Dict(
            "type" => "object",
            "properties" => Dict(
                "store" => Dict(
                    "type" => "string",
                    "description" => "Path to the _targets store directory. Required."
                ),
                "names" => Dict(
                    "type" => "string",
                    "description" => "Comma-separated target names to get metadata for (omit for all targets)"
                ),
                "project" => Dict(
                    "type" => "string",
                    "description" => "Project directory path"
                )
            ),
            "required" => ["store"]
        )
    ),
    Dict(
        "name" => "tar_outdated",
        "description" => "Check which targets are outdated and will rerun on tar_make(). A target is outdated if: its command changed, a dependency changed, it errored last time, or it never ran. Returns vector of target names. Empty vector means all targets are up-to-date.",
        "inputSchema" => Dict(
            "type" => "object",
            "properties" => Dict(
                "store" => Dict(
                    "type" => "string",
                    "description" => "Path to the _targets store directory. Required."
                ),
                "project" => Dict(
                    "type" => "string",
                    "description" => "Project directory path"
                )
            ),
            "required" => ["store"]
        )
    ),
    Dict(
        "name" => "tar_manifest",
        "description" => "List all targets defined in the pipeline with their deparsed commands. Shows the planned execution before running. Useful for understanding pipeline structure, reviewing target definitions, and seeing what each target computes.",
        "inputSchema" => Dict(
            "type" => "object",
            "properties" => Dict(
                "store" => Dict(
                    "type" => "string",
                    "description" => "Path to the _targets store directory. Required."
                ),
                "project" => Dict(
                    "type" => "string",
                    "description" => "Project directory path"
                )
            ),
            "required" => ["store"]
        )
    ),
    Dict(
        "name" => "tar_read",
        "description" => "Read a target's return value from _targets/objects/. For file targets (format='file'), returns the file paths instead of contents. The value is loaded into the session for further exploration with eval_expr. For dynamic branching patterns, use 'branches' to select specific branch indices. Large objects (>1MB) return a summary with hints for accessing specific parts.",
        "inputSchema" => Dict(
            "type" => "object",
            "properties" => Dict(
                "name" => Dict(
                    "type" => "string",
                    "description" => "Name of the target to read"
                ),
                "store" => Dict(
                    "type" => "string",
                    "description" => "Path to the _targets store directory. Required."
                ),
                "branches" => Dict(
                    "type" => "string",
                    "description" => "Comma-separated branch indices for dynamic branching patterns (e.g., '1,2,3')"
                ),
                "project" => Dict(
                    "type" => "string",
                    "description" => "Project directory path"
                )
            ),
            "required" => ["name", "store"]
        )
    ),
    Dict(
        "name" => "tar_load",
        "description" => "Load target return values into the session environment. Unlike tar_read which returns the value, tar_load assigns it to a variable with the target's name. For file targets (format='file'), loads the paths. Use 'branches' for dynamic branching patterns. After loading, use list_objects to see loaded targets and eval_expr to interact with them.",
        "inputSchema" => Dict(
            "type" => "object",
            "properties" => Dict(
                "names" => Dict(
                    "type" => "string",
                    "description" => "Comma-separated target names to load"
                ),
                "store" => Dict(
                    "type" => "string",
                    "description" => "Path to the _targets store directory. Required."
                ),
                "project" => Dict(
                    "type" => "string",
                    "description" => "Project directory path"
                )
            ),
            "required" => ["names", "store"]
        )
    ),
    Dict(
        "name" => "list_objects",
        "description" => "List all objects currently loaded in the session. Shows object names, classes, sizes, and whether they came from targets. Use after tar_load to see what's available.",
        "inputSchema" => Dict(
            "type" => "object",
            "properties" => Dict(
                "project" => Dict(
                    "type" => "string",
                    "description" => "Project directory path"
                )
            )
        )
    ),
    Dict(
        "name" => "eval_expr",
        "description" => "Evaluate an R expression in the session environment. Use this to inspect loaded objects, run summaries, extract specific fields, etc. Examples: 'head(my_data)', 'fit\$summary()', 'names(result)'",
        "inputSchema" => Dict(
            "type" => "object",
            "properties" => Dict(
                "expr" => Dict(
                    "type" => "string",
                    "description" => "R expression to evaluate"
                ),
                "project" => Dict(
                    "type" => "string",
                    "description" => "Project directory path"
                )
            ),
            "required" => ["expr"]
        )
    ),
    Dict(
        "name" => "reset_session",
        "description" => "Clear all loaded objects from the session. Use when you want to start fresh or free memory.",
        "inputSchema" => Dict(
            "type" => "object",
            "properties" => Dict(
                "project" => Dict(
                    "type" => "string",
                    "description" => "Project directory path"
                )
            )
        )
    ),
    Dict(
        "name" => "list_sessions",
        "description" => "List all active sessions with their loaded objects and store paths. Useful for understanding current state across projects.",
        "inputSchema" => Dict(
            "type" => "object",
            "properties" => Dict()
        )
    ),
    Dict(
        "name" => "tar_invalidate",
        "description" => "Delete metadata records to force targets to rerun on next tar_make(). Keeps the cached return values in _targets/objects/ (can be manually recovered) but tar_read/tar_load won't work until targets rerun. For patterns, invalidates all branches. Use when you need to force a rerun regardless of whether code/data changed.",
        "inputSchema" => Dict(
            "type" => "object",
            "properties" => Dict(
                "names" => Dict(
                    "type" => "string",
                    "description" => "Comma-separated target names to invalidate"
                ),
                "store" => Dict(
                    "type" => "string",
                    "description" => "Path to the _targets store directory. Required."
                ),
                "project" => Dict(
                    "type" => "string",
                    "description" => "Project directory path"
                )
            ),
            "required" => ["names", "store"]
        )
    ),
    Dict(
        "name" => "list_pipelines",
        "description" => "List all currently running targets pipelines. Shows controller PIDs, store paths, Stan model types, and CPU usage. Use this to see what's running before deciding to kill a pipeline.",
        "inputSchema" => Dict(
            "type" => "object",
            "properties" => Dict()
        )
    ),
    Dict(
        "name" => "kill_pipeline",
        "description" => "Kill a running targets pipeline. Terminates the controller process and any orphaned Stan processes. Provide either the controller PID or the store path.",
        "inputSchema" => Dict(
            "type" => "object",
            "properties" => Dict(
                "pid" => Dict(
                    "type" => "integer",
                    "description" => "Controller process ID to kill"
                ),
                "store" => Dict(
                    "type" => "string",
                    "description" => "Path to _targets store (used to find controller if pid not provided)"
                )
            )
        )
    ),
    Dict(
        "name" => "tar_list_variations",
        "description" => "Parse sclc_targets.R and extract tar_map variations (models, dcos, fit types). Returns structured data about available variations and target naming patterns. Use this to understand what target variations are available and construct valid target names.",
        "inputSchema" => Dict(
            "type" => "object",
            "properties" => Dict(
                "project" => Dict(
                    "type" => "string",
                    "description" => "Project directory path (defaults to /mnt/code)"
                )
            )
        )
    )
]

# Endpoints that use GET (read-only, no body parameters)
const GET_ENDPOINTS = Set(["list_pipelines", "list_sessions", "health"])
