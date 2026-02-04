#!/usr/bin/env Rscript
# MCP Wrapper for targets-mcp
# Translates MCP JSON-RPC (stdin/stdout) to HTTP calls to plumber backend

# CRITICAL: Handle initialize request IMMEDIATELY before any heavy loading
# This is necessary because Claude Code has a very short timeout for MCP servers

# Track if we already handled initialize
.INIT_HANDLED <- FALSE

# Open stdin connection ONCE at the start and keep it
.STDIN_CON <- file("stdin", "r", blocking = TRUE)

# Check if first message is initialize and handle it immediately
first_line <- readLines(.STDIN_CON, n = 1, warn = FALSE)
if (length(first_line) > 0 && grepl('"method"\\s*:\\s*"initialize"', first_line)) {
  # Fast path: respond to initialize immediately
  response <- '{"jsonrpc":"2.0","id":0,"result":{"protocolVersion":"2024-11-05","capabilities":{"tools":{}},"serverInfo":{"name":"targets-mcp","version":"0.1.0"}}}'
  cat(response, "\n", sep = "")
  flush(stdout())
  .INIT_HANDLED <- TRUE
} else {
  # Not an initialize request - we'll need to process it later
  .FIRST_LINE <- first_line
}

# Now do the heavy initialization
# Redirect stderr to log file for debugging
log_file <- "/tmp/targets-mcp-wrapper.log"
log_con <- file(log_file, open = "a")
sink(log_con, type = "message", append = TRUE)
cat(sprintf("\n=== MCP Wrapper started at %s ===\n", Sys.time()), file = stderr())
cat(sprintf("[%s] PID: %d\n", Sys.time(), Sys.getpid()), file = stderr())
cat(sprintf("[%s] Working directory: %s\n", Sys.time(), getwd()), file = stderr())
cat(sprintf("[%s] TAR_PROJECT: %s\n", Sys.time(), Sys.getenv("TAR_PROJECT", unset = "not set")), file = stderr())
cat(sprintf("[%s] Fast initialize response sent\n", Sys.time()), file = stderr())

# Temporarily redirect stdout to stderr to capture any stray messages
stdout_sink <- file(log_file, open = "a")
sink(stdout_sink, type = "output")

# Activate renv if it exists
renv_activate <- file.path("/mnt/code", "renv", "activate.R")
if (file.exists(renv_activate)) {
  cat(sprintf("[%s] Activating renv...\n", Sys.time()), file = stderr())
  # Suppress all output from renv activation
  suppressMessages(suppressWarnings(source(renv_activate, local = TRUE)))
  cat(sprintf("[%s] renv activated\n", Sys.time()), file = stderr())
} else {
  cat(sprintf("[%s] WARNING: renv not found at %s\n", Sys.time(), renv_activate), file = stderr())
}

# Disable conflicted package messages completely
options(conflicts.policy = list(warn = FALSE))

suppressPackageStartupMessages({
  library(jsonlite)
  library(httr)
})

# Restore stdout for JSON-RPC communication
sink(type = "output")
close(stdout_sink)

cat(sprintf("[%s] Packages loaded\n", Sys.time()), file = stderr())

# Configuration
PLUMBER_PORT <- 8484
PLUMBER_HOST <- "127.0.0.1"
PLUMBER_URL <- sprintf("http://%s:%d", PLUMBER_HOST, PLUMBER_PORT)

# Get script directory (same logic as start_server.R)
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    script_path <- sub("^--file=", "", file_arg)
    if (!startsWith(script_path, "/")) {
      script_path <- file.path(getwd(), script_path)
    }
    return(dirname(normalizePath(script_path, mustWork = FALSE)))
  }
  getwd()
}

# Tool definitions for MCP
TOOLS <- list(
  list(
    name = "tar_progress",
    description = "Read target progress from the most recent tar_make() run. Returns a data frame with columns: name, type (stem/pattern/branch), parent, branches count, and progress status. Progress values: 'dispatched' (sent to run, may be queued), 'completed', 'skipped', 'canceled', or 'errored'. Use to check pipeline status before loading targets.",
    inputSchema = list(
      type = "object",
      properties = list(
        store = list(
          type = "string",
          description = "Path to the _targets store directory. Required - do not rely on defaults."
        ),
        project = list(
          type = "string",
          description = "Project directory path (defaults to working directory)"
        )
      ),
      required = list("store")
    )
  ),
  list(
    name = "tar_meta",
    description = "Read metadata for all recorded targets. Returns columns: name, type (function/object/stem/branch/map/cross), data hash, command hash, depend hash, seed, path, time (last modified), size hash, bytes, format, error, warnings, seconds (runtime). Use to debug failures (check 'error' column), understand performance (check 'seconds', 'bytes'), or verify target state.",
    inputSchema = list(
      type = "object",
      properties = list(
        store = list(
          type = "string",
          description = "Path to the _targets store directory. Required."
        ),
        names = list(
          type = "string",
          description = "Comma-separated target names to get metadata for (omit for all targets)"
        ),
        project = list(
          type = "string",
          description = "Project directory path"
        )
      ),
      required = list("store")
    )
  ),
  list(
    name = "tar_outdated",
    description = "Check which targets are outdated and will rerun on tar_make(). A target is outdated if: its command changed, a dependency changed, it errored last time, or it never ran. Returns vector of target names. Empty vector means all targets are up-to-date.",
    inputSchema = list(
      type = "object",
      properties = list(
        store = list(
          type = "string",
          description = "Path to the _targets store directory. Required."
        ),
        project = list(
          type = "string",
          description = "Project directory path"
        )
      ),
      required = list("store")
    )
  ),
  list(
    name = "tar_manifest",
    description = "List all targets defined in the pipeline with their deparsed commands. Shows the planned execution before running. Useful for understanding pipeline structure, reviewing target definitions, and seeing what each target computes.",
    inputSchema = list(
      type = "object",
      properties = list(
        store = list(
          type = "string",
          description = "Path to the _targets store directory. Required."
        ),
        project = list(
          type = "string",
          description = "Project directory path"
        )
      ),
      required = list("store")
    )
  ),
  list(
    name = "tar_read",
    description = "Read a target's return value from _targets/objects/. For file targets (format='file'), returns the file paths instead of contents. The value is loaded into the session for further exploration with eval_expr. For dynamic branching patterns, use 'branches' to select specific branch indices. Large objects (>1MB) return a summary with hints for accessing specific parts.",
    inputSchema = list(
      type = "object",
      properties = list(
        name = list(
          type = "string",
          description = "Name of the target to read"
        ),
        store = list(
          type = "string",
          description = "Path to the _targets store directory. Required."
        ),
        branches = list(
          type = "string",
          description = "Comma-separated branch indices for dynamic branching patterns (e.g., '1,2,3')"
        ),
        project = list(
          type = "string",
          description = "Project directory path"
        )
      ),
      required = list("name", "store")
    )
  ),
  list(
    name = "tar_load",
    description = "Load target return values into the session environment. Unlike tar_read which returns the value, tar_load assigns it to a variable with the target's name. For file targets (format='file'), loads the paths. Use 'branches' for dynamic branching patterns. After loading, use list_objects to see loaded targets and eval_expr to interact with them.",
    inputSchema = list(
      type = "object",
      properties = list(
        names = list(
          type = "string",
          description = "Comma-separated target names to load"
        ),
        store = list(
          type = "string",
          description = "Path to the _targets store directory. Required."
        ),
        project = list(
          type = "string",
          description = "Project directory path"
        )
      ),
      required = list("names", "store")
    )
  ),
  list(
    name = "list_objects",
    description = "List all objects currently loaded in the session. Shows object names, classes, sizes, and whether they came from targets. Use after tar_load to see what's available.",
    inputSchema = list(
      type = "object",
      properties = list(
        project = list(
          type = "string",
          description = "Project directory path"
        )
      )
    )
  ),
  list(
    name = "eval_expr",
    description = "Evaluate an R expression in the session environment. Use this to inspect loaded objects, run summaries, extract specific fields, etc. Examples: 'head(my_data)', 'fit$summary()', 'names(result)'",
    inputSchema = list(
      type = "object",
      properties = list(
        expr = list(
          type = "string",
          description = "R expression to evaluate"
        ),
        project = list(
          type = "string",
          description = "Project directory path"
        )
      ),
      required = list("expr")
    )
  ),
  list(
    name = "reset_session",
    description = "Clear all loaded objects from the session. Use when you want to start fresh or free memory.",
    inputSchema = list(
      type = "object",
      properties = list(
        project = list(
          type = "string",
          description = "Project directory path"
        )
      )
    )
  ),
  list(
    name = "list_sessions",
    description = "List all active sessions with their loaded objects and store paths. Useful for understanding current state across projects.",
    inputSchema = list(
      type = "object",
      properties = list()
    )
  ),
  list(
    name = "tar_invalidate",
    description = "Delete metadata records to force targets to rerun on next tar_make(). Keeps the cached return values in _targets/objects/ (can be manually recovered) but tar_read/tar_load won't work until targets rerun. For patterns, invalidates all branches. Use when you need to force a rerun regardless of whether code/data changed.",
    inputSchema = list(
      type = "object",
      properties = list(
        names = list(
          type = "string",
          description = "Comma-separated target names to invalidate"
        ),
        store = list(
          type = "string",
          description = "Path to the _targets store directory. Required."
        ),
        project = list(
          type = "string",
          description = "Project directory path"
        )
      ),
      required = list("names", "store")
    )
  ),
  list(
    name = "list_pipelines",
    description = "List all currently running targets pipelines. Shows controller PIDs, store paths, Stan model types, and CPU usage. Use this to see what's running before deciding to kill a pipeline.",
    inputSchema = list(
      type = "object",
      properties = list()
    )
  ),
  list(
    name = "kill_pipeline",
    description = "Kill a running targets pipeline. Terminates the controller process and any orphaned Stan processes. Provide either the controller PID or the store path.",
    inputSchema = list(
      type = "object",
      properties = list(
        pid = list(
          type = "integer",
          description = "Controller process ID to kill"
        ),
        store = list(
          type = "string",
          description = "Path to _targets store (used to find controller if pid not provided)"
        )
      )
    )
  )
)

# Server info
SERVER_INFO <- list(
  name = "targets-mcp",
  version = "0.1.0"
)

# Find process using a specific port
find_process_on_port <- function(port) {
  # Try to find process using ps and grep
  tryCatch({
    # Look for R processes that might be running the server
    cmd <- sprintf("ps aux | grep -E 'start_server.*%d|plumber.*%d' | grep -v grep | awk '{print $2}'", port, port)
    result <- system(cmd, intern = TRUE, ignore.stderr = TRUE)
    if (length(result) > 0 && nzchar(result[1])) {
      return(as.integer(result[1]))
    }
  }, error = function(e) {
    # Ignore errors
  })
  return(NULL)
}

# Kill stale plumber process
kill_stale_process <- function(port) {
  pid <- find_process_on_port(port)
  if (!is.null(pid)) {
    tryCatch({
      system2("kill", args = as.character(pid), stdout = FALSE, stderr = FALSE)
      Sys.sleep(1)  # Give it time to die
      return(TRUE)
    }, error = function(e) {
      return(FALSE)
    })
  }
  return(FALSE)
}

# Start plumber server if not running, reuse if healthy
ensure_plumber_running <- function() {
  cat(sprintf("[%s] Checking if plumber server is running...\n", Sys.time()), file = stderr())

  # Check if server is responding and healthy
  tryCatch({
    response <- httr::GET(paste0(PLUMBER_URL, "/health"), httr::timeout(2))
    if (httr::status_code(response) == 200) {
      # Server is healthy, reuse it
      cat(sprintf("[%s] Server is healthy, reusing existing session\n", Sys.time()), file = stderr())
      return(TRUE)
    }
  }, error = function(e) {
    # Server not responding, may need to clean up
    cat(sprintf("[%s] Server not responding: %s\n", Sys.time(), conditionMessage(e)), file = stderr())
  })

  # Server not healthy - check if port is in use by stale process
  kill_stale_process(PLUMBER_PORT)

  # Start plumber in background
  script_dir <- get_script_dir()
  start_script <- file.path(script_dir, "start_server.R")

  # Use system2 to start in background
  system2(
    "Rscript",
    args = c(start_script, as.character(PLUMBER_PORT)),
    wait = FALSE,
    stdout = FALSE,
    stderr = FALSE
  )

  # Wait for server to start (with shorter intervals)
  max_attempts <- 20
  for (i in 1:max_attempts) {
    Sys.sleep(0.5)  # Check every 500ms instead of 1s
    tryCatch({
      response <- httr::GET(paste0(PLUMBER_URL, "/health"), httr::timeout(1))
      if (httr::status_code(response) == 200) {
        cat(sprintf("[%s] Server started successfully after %d attempts\n", Sys.time(), i), file = stderr())
        return(TRUE)
      }
    }, error = function(e) {
      # Keep waiting
      if (i %% 4 == 0) {  # Log every 2 seconds
        cat(sprintf("[%s] Still waiting for server (attempt %d/%d)...\n", Sys.time(), i, max_attempts), file = stderr())
      }
    })
  }

  cat(sprintf("[%s] FATAL: Failed to start plumber server after %d attempts\n", Sys.time(), max_attempts), file = stderr())
  stop("Failed to start plumber server")
}

# Endpoints that use GET (read-only, no parameters)
GET_ENDPOINTS <- c("list_pipelines", "list_sessions", "health")

# Call plumber endpoint
call_plumber <- function(endpoint, params = list()) {
  url <- paste0(PLUMBER_URL, "/", endpoint)

  # Use GET for read-only endpoints, POST for others
  if (endpoint %in% GET_ENDPOINTS) {
    response <- httr::GET(url)
  } else {
    response <- httr::POST(
      url,
      body = params,
      encode = "json",
      httr::content_type_json()
    )
  }

  content <- httr::content(response, as = "parsed", simplifyVector = TRUE)
  content
}

# Handle MCP initialize request
handle_initialize <- function(params) {
  # Use client's protocol version or default
  protocol_version <- if (!is.null(params$protocolVersion)) params$protocolVersion else "2024-11-05"
  cat(sprintf("[%s] Client protocol version: %s\n", Sys.time(), protocol_version), file = stderr())
  cat(sprintf("[%s] Client capabilities: %s\n", Sys.time(), toJSON(params$capabilities, auto_unbox = TRUE)), file = stderr())

  # Return proper MCP initialize response with all required capabilities
  list(
    protocolVersion = protocol_version,
    capabilities = list(
      tools = list(),
      prompts = list(),
      resources = list()
    ),
    serverInfo = SERVER_INFO
  )
}

# Handle tools/list request
handle_tools_list <- function() {
  list(tools = TOOLS)
}

# Handle tools/call request
handle_tools_call <- function(params) {
  tool_name <- params$name
  tool_args <- params$arguments

  if (is.null(tool_args)) {
    tool_args <- list()
  }

  # Map tool name to plumber endpoint
  endpoint <- tool_name

  # Call plumber
  result <- tryCatch({
    call_plumber(endpoint, tool_args)
  }, error = function(e) {
    list(success = FALSE, error = conditionMessage(e))
  })

  # Format result for MCP
  if (isTRUE(result$success)) {
    list(
      content = list(
        list(
          type = "text",
          text = toJSON(result, auto_unbox = TRUE, pretty = TRUE, null = "null")
        )
      )
    )
  } else {
    list(
      content = list(
        list(
          type = "text",
          text = paste("Error:", result$error %||% "Unknown error")
        )
      ),
      isError = TRUE
    )
  }
}

# Process a single JSON-RPC request
process_request <- function(request) {
  method <- request$method
  params <- request$params
  id <- request$id

  result <- tryCatch({
    switch(method,
      "initialize" = handle_initialize(params),
      "tools/list" = handle_tools_list(),
      "tools/call" = handle_tools_call(params),
      "notifications/initialized" = NULL,  # No response needed
      "ping" = list(),  # Empty response for ping
      stop(sprintf("Unknown method: %s", method))
    )
  }, error = function(e) {
    list(
      code = -32603,
      message = conditionMessage(e)
    )
  })

  # Don't send response for notifications
  if (is.null(id)) {
    return(NULL)
  }

  # Build response
  if (!is.null(result$code)) {
    # Error response
    list(
      jsonrpc = "2.0",
      id = id,
      error = result
    )
  } else {
    # Success response
    list(
      jsonrpc = "2.0",
      id = id,
      result = result
    )
  }
}

# Main loop - read from stdin, write to stdout
main <- function() {
  cat(sprintf("[%s] main() starting...\n", Sys.time()), file = stderr())

  # Ensure plumber is running
  tryCatch({
    ensure_plumber_running()
    cat(sprintf("[%s] Plumber server ready\n", Sys.time()), file = stderr())
  }, error = function(e) {
    cat(sprintf("[%s] FATAL: Failed to ensure plumber running: %s\n", Sys.time(), conditionMessage(e)), file = stderr())
    stop(e)
  })

  cat(sprintf("[%s] Entering main loop, waiting for stdin...\n", Sys.time()), file = stderr())

  # Use the global stdin connection we opened at startup
  con <- .STDIN_CON

  # If we have a buffered first line that wasn't initialize, process it first
  pending_line <- if (exists(".FIRST_LINE") && length(.FIRST_LINE) > 0) .FIRST_LINE else NULL

  while (TRUE) {
    # Try to read a line, handling errors gracefully
    if (!is.null(pending_line)) {
      line <- pending_line
      pending_line <- NULL
      cat(sprintf("[%s] Processing buffered line\n", Sys.time()), file = stderr())
    } else {
      line <- tryCatch({
        readLines(con, n = 1, warn = FALSE)
      }, error = function(e) {
        cat(sprintf("[%s] Error reading from stdin: %s\n", Sys.time(), conditionMessage(e)), file = stderr())
        character(0)
      })
    }

    cat(sprintf("[%s] Read from stdin: length=%d\n", Sys.time(), length(line)), file = stderr())

    if (length(line) == 0) {
      # EOF - this shouldn't happen in normal operation
      cat(sprintf("[%s] EOF received - connection closed by client\n", Sys.time()), file = stderr())
      break
    }

    if (nchar(trimws(line)) == 0) {
      cat(sprintf("[%s] Empty line, skipping\n", Sys.time()), file = stderr())
      next
    }

    cat(sprintf("[%s] Received: %s\n", Sys.time(), substr(line, 1, 100)), file = stderr())

    # Parse JSON-RPC request
    request <- tryCatch({
      fromJSON(line, simplifyVector = FALSE)
    }, error = function(e) {
      cat(sprintf("[%s] JSON parse error: %s\n", Sys.time(), conditionMessage(e)), file = stderr())
      NULL
    })

    if (is.null(request)) {
      # Invalid JSON, skip
      cat(sprintf("[%s] Invalid JSON, skipping\n", Sys.time()), file = stderr())
      next
    }

    cat(sprintf("[%s] Processing method: %s\n", Sys.time(), request$method), file = stderr())

    # Skip initialize if we already handled it at startup
    if (request$method == "initialize" && .INIT_HANDLED) {
      cat(sprintf("[%s] Skipping initialize (already handled at startup)\n", Sys.time()), file = stderr())
      next
    }

    # Process request
    response <- process_request(request)

    if (!is.null(response)) {
      # Write response to stdout
      response_json <- toJSON(response, auto_unbox = TRUE, null = "null")
      cat(sprintf("[%s] Sending response: %s\n", Sys.time(), substr(response_json, 1, 100)), file = stderr())
      cat(response_json, "\n", sep = "")
      flush(stdout())
      cat(sprintf("[%s] Response sent and flushed\n", Sys.time()), file = stderr())
    } else {
      cat(sprintf("[%s] No response to send (notification)\n", Sys.time()), file = stderr())
    }
  }
}

# Run main
tryCatch({
  main()
}, error = function(e) {
  cat(sprintf("[%s] FATAL ERROR: %s\n", Sys.time(), conditionMessage(e)), file = stderr())
  cat(sprintf("[%s] Traceback:\n", Sys.time()), file = stderr())
  print(traceback(), file = stderr())
  quit(status = 1)
})
