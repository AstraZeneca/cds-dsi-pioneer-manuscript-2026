# Session Manager for targets-mcp
# Manages per-project R environments for stateful target exploration

# Global session storage
sessions <- new.env(parent = emptyenv())

# Session timeout in seconds (30 minutes)
SESSION_TIMEOUT <- 30 * 60

#' Find the targets store for a project
#'
#' @param project_path Path to the project directory
#' @return Path to the _targets store, or NULL if not found
find_store <- function(project_path) {
  project_path <- normalizePath(project_path, mustWork = FALSE)

  # 1. Check _targets.yaml in project directory
  yaml_path <- file.path(project_path, "_targets.yaml")
  if (file.exists(yaml_path)) {
    yaml_content <- yaml::read_yaml(yaml_path)
    # Check for store in main config or project-specific config
    if (!is.null(yaml_content$store)) {
      store <- yaml_content$store
      if (!startsWith(store, "/")) {
        store <- file.path(project_path, store)
      }
      if (dir.exists(store)) {
        return(normalizePath(store))
      }
    }
  }

  # 2. Check default _targets/ in project directory
  default_store <- file.path(project_path, "_targets")
  if (dir.exists(default_store)) {
    return(normalizePath(default_store))
  }

  # Not found - return NULL (caller should handle error)
  # Note: We intentionally don't have hardcoded fallback paths.
  # The store must be explicitly provided or discoverable from _targets.yaml
  NULL
}

#' Get or create a session for a project/store combination
#'
#' @param project_path Path to the project directory
#' @param store Optional explicit store path (overrides discovery)
#' @return Session environment
get_session <- function(project_path, store = NULL) {
  project_path <- normalizePath(project_path, mustWork = FALSE)

  # Use store path as part of key if provided (allows multiple stores per project)
  store_for_key <- if (!is.null(store) && nzchar(store)) {
    normalizePath(store, mustWork = FALSE)
  } else {
    find_store(project_path)
  }

  # Key includes both project and store to allow multiple stores
  key <- paste0(project_path, "::", store_for_key)

  if (!exists(key, envir = sessions)) {
    # Create new session
    session <- new.env(parent = globalenv())
    session$.project <- project_path
    session$.store <- store_for_key
    session$.created <- Sys.time()
    session$.last_used <- Sys.time()

    # Track loaded objects (names only, for list_objects)
    session$.loaded_targets <- character(0)

    sessions[[key]] <- session
  }

  # Update last used time
  sessions[[key]]$.last_used <- Sys.time()

  sessions[[key]]
}

#' Get session info as a list (for JSON serialization)
#'
#' @param session Session environment
#' @return List with session info
session_info <- function(session) {
  # Get user-defined objects (exclude . prefixed metadata)
  all_names <- ls(session, all.names = FALSE)

  list(
    project = session$.project,
    store = session$.store,
    objects_loaded = session$.loaded_targets,
    all_objects = all_names,
    created = format(session$.created, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    last_used = format(session$.last_used, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
  )
}

#' List all active sessions
#'
#' @return List of session info objects
list_sessions <- function() {
  keys <- ls(sessions, all.names = TRUE)
  lapply(keys, function(key) {
    session_info(sessions[[key]])
  })
}

#' Reset a session (clear all loaded objects)
#'
#' @param project_path Path to the project directory
#' @return TRUE if session was reset, FALSE if not found
reset_session <- function(project_path) {
  project_path <- normalizePath(project_path, mustWork = FALSE)
  key <- project_path


  if (!exists(key, envir = sessions)) {
    return(FALSE)
  }

  session <- sessions[[key]]

  # Remove all user objects (keep metadata)
  user_objects <- ls(session, all.names = FALSE)
  rm(list = user_objects, envir = session)

  # Reset loaded targets tracking
  session$.loaded_targets <- character(0)
  session$.last_used <- Sys.time()

  TRUE
}

#' Destroy a session completely
#'
#' @param project_path Path to the project directory
#' @return TRUE if session was destroyed, FALSE if not found
destroy_session <- function(project_path) {
  project_path <- normalizePath(project_path, mustWork = FALSE)
  key <- project_path

  if (!exists(key, envir = sessions)) {
    return(FALSE)
  }

  rm(list = key, envir = sessions)
  TRUE
}

#' Clean up stale sessions (idle > timeout)
#'
#' @param timeout_seconds Timeout in seconds (default: SESSION_TIMEOUT)
#' @return Number of sessions cleaned up
cleanup_stale_sessions <- function(timeout_seconds = SESSION_TIMEOUT) {
  keys <- ls(sessions, all.names = TRUE)
  now <- Sys.time()
  cleaned <- 0

  for (key in keys) {
    session <- sessions[[key]]
    idle_time <- as.numeric(difftime(now, session$.last_used, units = "secs"))

    if (idle_time > timeout_seconds) {
      rm(list = key, envir = sessions)
      cleaned <- cleaned + 1
    }
  }

  cleaned
}

#' Require a store path for a session
#'
#' @param session Session environment
#' @return Store path or stop with error
require_store <- function(session) {
  if (is.null(session$.store) || !dir.exists(session$.store)) {
    stop(sprintf(
      "No targets store found for project '%s'. Please provide an explicit store path.",
      session$.project
    ))
  }
  session$.store
}
