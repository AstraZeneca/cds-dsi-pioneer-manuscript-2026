# targets-mcp Plumber API
# Exposes {targets} functionality as HTTP endpoints

library(plumber)
library(targets)
library(jsonlite)

# Helper to get script directory
get_script_dir <- function() {
  # Method 1: Check for TARGETS_MCP_DIR env var FIRST (set by start_server.R)
  env_dir <- Sys.getenv("TARGETS_MCP_DIR", "")
  if (env_dir != "" && dir.exists(env_dir)) {
    return(env_dir)
  }

  # Method 2: Try to derive from this file's path if source'd
  for (i in seq_len(sys.nframe())) {
    ofile <- sys.frame(i)$ofile
    if (!is.null(ofile)) {
      return(dirname(normalizePath(ofile, mustWork = FALSE)))
    }
  }

  # Method 3: Fallback to known path (used during plumber parsing)
  known_path <- "/mnt/code/.claude/mcp/targets-mcp"
  if (dir.exists(known_path)) {
    return(known_path)
  }

  # Last resort: current directory
  getwd()
}

# Source session manager and process manager
script_dir <- get_script_dir()
source(file.path(script_dir, "sessions.R"))
source(file.path(script_dir, "process.R"))

# Default size limit for returned objects (1MB)
SIZE_LIMIT_BYTES <- 1024 * 1024

#* @apiTitle targets-mcp API
#* @apiDescription MCP server for {targets} package - provides stateful R sessions for target exploration

#* Health check
#* @get /health
function() {
  list(
    status = "ok",
    timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    sessions_active = length(ls(sessions, all.names = TRUE))
  )
}

#* Get pipeline progress
#* @param store Path to _targets store
#* @param project Project directory (for session, defaults to working directory)
#* @post /tar_progress
function(store = NULL, project = getwd()) {
  session <- get_session(project, store)
  store_path <- require_store(session)

  tryCatch({
    progress <- targets::tar_progress(store = store_path)

    # Convert to list for JSON serialization
    if (nrow(progress) == 0) {
      list(
        success = TRUE,
        data = list(),
        message = "No targets in progress"
      )
    } else {
      list(
        success = TRUE,
        data = as.list(as.data.frame(progress)),
        row_count = nrow(progress)
      )
    }
  }, error = function(e) {
    list(
      success = FALSE,
      error = conditionMessage(e)
    )
  })
}

#* Get target metadata
#* @param store Path to _targets store
#* @param project Project directory
#* @param names Target names to get metadata for (comma-separated, or NULL for all)
#* @param fields Fields to return (comma-separated, or NULL for default)
#* @post /tar_meta
function(store = NULL, project = getwd(), names = NULL, fields = NULL) {
  session <- get_session(project, store)
  store_path <- require_store(session)

  tryCatch({
    # Parse names if provided
    target_names <- if (!is.null(names) && names != "") {
      strsplit(names, ",")[[1]] |> trimws()
    } else {
      NULL
    }

    # Parse fields if provided
    target_fields <- if (!is.null(fields) && fields != "") {
      strsplit(fields, ",")[[1]] |> trimws() |> as.symbol() |> (\(x) substitute(c(x)))()
    } else {
      NULL
    }

    meta <- targets::tar_meta(
      names = target_names,
      store = store_path
    )

    if (nrow(meta) == 0) {
      list(
        success = TRUE,
        data = list(),
        message = "No metadata found"
      )
    } else {
      # Select useful columns for JSON output
      useful_cols <- intersect(
        c("name", "type", "format", "bytes", "seconds", "error", "warnings", "time"),
        names(meta)
      )
      meta_subset <- meta[, useful_cols, drop = FALSE]

      list(
        success = TRUE,
        data = as.list(as.data.frame(meta_subset)),
        row_count = nrow(meta_subset)
      )
    }
  }, error = function(e) {
    list(
      success = FALSE,
      error = conditionMessage(e)
    )
  })
}

#* Show outdated targets
#* @param store Path to _targets store
#* @param project Project directory
#* @post /tar_outdated
function(store = NULL, project = getwd()) {
  session <- get_session(project, store)
  store_path <- require_store(session)

  tryCatch({
    outdated <- targets::tar_outdated(store = store_path, callr_function = NULL)

    list(
      success = TRUE,
      data = outdated,
      count = length(outdated),
      message = if (length(outdated) == 0) "All targets up to date" else NULL
    )
  }, error = function(e) {
    list(
      success = FALSE,
      error = conditionMessage(e)
    )
  })
}

#* Get target manifest
#* @param store Path to _targets store
#* @param project Project directory
#* @post /tar_manifest
function(store = NULL, project = getwd()) {
  session <- get_session(project, store)
  store_path <- require_store(session)

  tryCatch({
    manifest <- targets::tar_manifest(store = store_path)

    if (nrow(manifest) == 0) {
      list(
        success = TRUE,
        data = list(),
        message = "No targets in manifest"
      )
    } else {
      # Return name and command columns
      useful_cols <- intersect(c("name", "command"), names(manifest))
      manifest_subset <- manifest[, useful_cols, drop = FALSE]

      list(
        success = TRUE,
        data = as.list(as.data.frame(manifest_subset)),
        row_count = nrow(manifest_subset)
      )
    }
  }, error = function(e) {
    list(
      success = FALSE,
      error = conditionMessage(e)
    )
  })
}

#* Read a target value
#* @param name Target name (required)
#* @param store Path to _targets store
#* @param project Project directory
#* @param branches Branch indices (comma-separated integers)
#* @post /tar_read
function(name, store = NULL, project = getwd(), branches = NULL) {
  if (missing(name) || is.null(name) || name == "") {
    return(list(success = FALSE, error = "Target name is required"))
  }

  session <- get_session(project, store)
  store_path <- require_store(session)

  tryCatch({
    # Parse branches if provided
    branch_indices <- if (!is.null(branches) && branches != "") {
      as.integer(strsplit(branches, ",")[[1]])
    } else {
      NULL
    }

    # Read the target
    value <- targets::tar_read(
      name = as.symbol(name),
      store = store_path,
      branches = branch_indices
    )

    # Store in session for later use
    assign(name, value, envir = session)
    session$.loaded_targets <- unique(c(session$.loaded_targets, name))

    # Check size and return appropriately
    serialized_size <- object.size(value)

    if (serialized_size > SIZE_LIMIT_BYTES) {
      # Large object - return summary
      summary_text <- capture.output(str(value, max.level = 2, give.attr = FALSE))

      list(
        success = TRUE,
        warning = sprintf(
          "Object '%s' is %s. Returning summary. Use eval_expr to access specific parts.",
          name,
          format(serialized_size, units = "auto")
        ),
        summary = paste(summary_text, collapse = "\n"),
        class = class(value),
        size = as.numeric(serialized_size),
        loaded_in_session = TRUE,
        hint = sprintf("Object is loaded in session. Use eval_expr with expressions like: %s$...", name)
      )
    } else {
      # Small enough to return directly
      list(
        success = TRUE,
        data = value,
        class = class(value),
        size = as.numeric(serialized_size),
        loaded_in_session = TRUE
      )
    }
  }, error = function(e) {
    list(
      success = FALSE,
      error = conditionMessage(e)
    )
  })
}

#* Load targets into session
#* @param names Target names (comma-separated, required)
#* @param store Path to _targets store
#* @param project Project directory
#* @post /tar_load
function(names, store = NULL, project = getwd()) {
  if (missing(names) || is.null(names) || names == "") {
    return(list(success = FALSE, error = "Target names are required"))
  }

  session <- get_session(project, store)
  store_path <- require_store(session)

  target_names <- strsplit(names, ",")[[1]] |> trimws()

  tryCatch({
    # Load each target into session
    targets::tar_load(
      names = target_names,
      store = store_path,
      envir = session
    )

    # Update loaded targets tracking
    session$.loaded_targets <- unique(c(session$.loaded_targets, target_names))

    list(
      success = TRUE,
      loaded = target_names,
      message = sprintf("Loaded %d target(s) into session", length(target_names))
    )
  }, error = function(e) {
    list(
      success = FALSE,
      error = conditionMessage(e)
    )
  })
}

#* List objects in session
#* @param project Project directory
#* @post /list_objects
function(project = getwd()) {
  session <- get_session(project)

  # Get all user objects (not metadata)
  all_objects <- ls(session, all.names = FALSE)

  object_info <- lapply(all_objects, function(obj_name) {
    obj <- get(obj_name, envir = session)
    list(
      name = obj_name,
      class = class(obj),
      size = as.numeric(object.size(obj)),
      is_target = obj_name %in% session$.loaded_targets
    )
  })

  list(
    success = TRUE,
    objects = object_info,
    count = length(all_objects),
    session_info = session_info(session)
  )
}

#* Evaluate R expression in session
#* @param expr R expression to evaluate (required)
#* @param project Project directory
#* @post /eval_expr
function(expr, project = getwd()) {
  if (missing(expr) || is.null(expr) || expr == "") {
    return(list(success = FALSE, error = "Expression is required"))
  }

  session <- get_session(project)

  tryCatch({
    # Parse and evaluate expression
    parsed <- parse(text = expr)
    result <- eval(parsed, envir = session)

    # Capture output if it's a print/summary call
    output <- capture.output(result)

    # Check size
    serialized_size <- object.size(result)

    if (serialized_size > SIZE_LIMIT_BYTES) {
      list(
        success = TRUE,
        warning = sprintf("Result is %s. Returning summary.", format(serialized_size, units = "auto")),
        output = paste(output, collapse = "\n"),
        class = class(result),
        size = as.numeric(serialized_size)
      )
    } else {
      list(
        success = TRUE,
        data = result,
        output = if (length(output) > 0) paste(output, collapse = "\n") else NULL,
        class = class(result),
        size = as.numeric(serialized_size)
      )
    }
  }, error = function(e) {
    list(
      success = FALSE,
      error = conditionMessage(e)
    )
  })
}

#* Reset session (clear loaded objects)
#* @param project Project directory
#* @post /reset_session
function(project = getwd()) {
  result <- reset_session(project)

  if (result) {
    list(
      success = TRUE,
      message = "Session reset successfully"
    )
  } else {
    list(
      success = FALSE,
      error = "Session not found"
    )
  }
}

#* List all active sessions
#* @get /list_sessions
function() {
  sessions_list <- list_sessions()

  list(
    success = TRUE,
    sessions = sessions_list,
    count = length(sessions_list)
  )
}

#* Invalidate targets
#* @param names Target names (comma-separated, required)
#* @param store Path to _targets store
#* @param project Project directory
#* @post /tar_invalidate
function(names, store = NULL, project = getwd()) {
  if (missing(names) || is.null(names) || names == "") {
    return(list(success = FALSE, error = "Target names are required"))
  }

  session <- get_session(project, store)
  store_path <- require_store(session)

  target_names <- strsplit(names, ",")[[1]] |> trimws()

  tryCatch({
    targets::tar_invalidate(
      names = target_names,
      store = store_path
    )

    list(
      success = TRUE,
      invalidated = target_names,
      message = sprintf("Invalidated %d target(s)", length(target_names))
    )
  }, error = function(e) {
    list(
      success = FALSE,
      error = conditionMessage(e)
    )
  })
}

#* List all running pipelines
#* @serializer unboxedJSON
#* @get /list_pipelines
function() {
  tryCatch({
    list_pipelines()
  }, error = function(e) {
    list(
      success = FALSE,
      error = conditionMessage(e)
    )
  })
}

#* Kill a running pipeline
#* @param pid Controller process ID (optional if store provided)
#* @param store Path to _targets store (optional if pid provided)
#* @post /kill_pipeline
function(pid = NULL, store = NULL) {
  # Convert pid to integer if provided
  if (!is.null(pid) && pid != "") {
    pid <- as.integer(pid)
  } else {
    pid <- NULL
  }

  tryCatch({
    kill_pipeline(pid = pid, store = store)
  }, error = function(e) {
    list(
      success = FALSE,
      error = conditionMessage(e)
    )
  })
}

#* List target variations from sclc_targets.R
#* Sources the targets file and traverses the resulting list structure
#* @param project Project directory (defaults to /mnt/code)
#* @post /tar_list_variations
function(project = "/mnt/code") {
  targets_file <- file.path(project, "targets", "sclc_targets.R")

  if (!file.exists(targets_file)) {
    return(list(
      success = FALSE,
      error = sprintf("sclc_targets.R not found at %s", targets_file)
    ))
  }

  tryCatch({
    # Load required packages
    suppressPackageStartupMessages({
      require(targets)
      require(tarchetypes)
      require(tidyverse)
      require(crew)
    })

    # Set up environment
    old_wd <- getwd()
    on.exit(setwd(old_wd))
    setwd(project)

    # Source .Rprofile for init_project() and other helpers
    if (file.exists(".Rprofile")) {
      source(".Rprofile", local = FALSE)
    }

    # Set required env var
    Sys.setenv(TAR_PROJECT = "sclc")

    # Source the targets file - this creates sclc_targets
    source(targets_file, local = FALSE)

    # Get all target names recursively
    get_all_target_names <- function(targets_list) {
      names_list <- character()
      for (item in targets_list) {
        if (inherits(item, "tar_target")) {
          names_list <- c(names_list, item$settings$name)
        } else if (is.list(item)) {
          names_list <- c(names_list, get_all_target_names(item))
        }
      }
      names_list
    }

    # Build tree structure from targets list
    build_tree <- function(targets_list, depth = 0, max_depth = 4) {
      result <- list()
      for (i in seq_along(targets_list)) {
        item <- targets_list[[i]]
        item_name <- names(targets_list)[i]
        if (is.null(item_name)) item_name <- ""

        if (inherits(item, "tar_target")) {
          result <- c(result, list(list(
            type = "target",
            name = item$settings$name
          )))
        } else if (is.list(item) && length(item) > 0) {
          children <- if (depth < max_depth) build_tree(item, depth + 1, max_depth) else list()
          result <- c(result, list(list(
            type = "group",
            name = item_name,
            count = length(item),
            children = children
          )))
        }
      }
      result
    }

    all_names <- get_all_target_names(sclc_targets)

    # Extract variations from target names by finding common patterns
    extract_variations <- function(names) {
      # Known dco values (extracted from target names ending in _apr or _aug)
      dcos <- unique(sub(".*_([^_]+)$", "\\1", names[grepl("_(apr|aug)$", names)]))

      # Known models - extract from patterns like tumor_ssls_lfo_{model}_{dco}
      # Model names can have underscores (e.g., no_oe, no_ctdna, ctdna_only)
      lfo_names <- names[grepl("^tumor_ssls_lfo_.*_(apr|aug)$", names)]
      # Remove prefix and suffix to get model name
      models <- unique(sub("^tumor_ssls_lfo_(.*)_(apr|aug)$", "\\1", lfo_names))
      # Filter out names that still contain lfo-related suffixes
      models <- models[!grepl("^(clean|endpoints)", models)]

      # Fit types
      fit_types <- c("prior", "posterior")

      list(
        dco = dcos,
        model = models,
        fit_type = fit_types
      )
    }

    variations <- extract_variations(all_names)

    # Get sample targets for key patterns (just a few for validation)
    sample_targets <- list(
      tumor_ssls_lfo = head(all_names[grepl("^tumor_ssls_lfo_[^_]+_(apr|aug)$", all_names)], 3),
      tumor_ssls_res = head(all_names[grepl("^tumor_ssls_res_", all_names)], 3),
      prior_tumor_ssls = head(all_names[grepl("^prior_tumor_ssls", all_names)], 3),
      tumor_ssls = head(all_names[grepl("^tumor_ssls_[^_]+_(apr|aug)$", all_names)], 3)
    )

    list(
      success = TRUE,
      total_targets = length(all_names),
      variations = variations,
      sample_targets = sample_targets
    )
  }, error = function(e) {
    list(
      success = FALSE,
      error = conditionMessage(e),
      traceback = paste(capture.output(traceback()), collapse = "\n")
    )
  })
}
