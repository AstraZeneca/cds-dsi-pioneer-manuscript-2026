# Process Manager for targets-mcp
# Discovers running targets pipelines and manages their lifecycle

#' Find all running Stan model processes
#'
#' @return Data frame with PID, PPID, command, and model type
find_stan_processes <- function() {
  # Get all processes matching Stan model executables
  cmd <- "ps aux | grep -E '(sf-ssls-lfo|sf-ssm-log-space)' | grep -v grep"
  output <- tryCatch({
    system(cmd, intern = TRUE)
  }, error = function(e) {
    character(0)
  })

  if (length(output) == 0) {
    return(data.frame(
      pid = integer(0),
      ppid = integer(0),
      cpu = numeric(0),
      cmd = character(0),
      model_type = character(0),
      output_file = character(0),
      stringsAsFactors = FALSE
    ))
  }

  # Parse ps output
  results <- lapply(output, function(line) {
    parts <- strsplit(trimws(line), "\\s+")[[1]]
    if (length(parts) < 11) return(NULL)

    pid <- as.integer(parts[2])
    cpu <- as.numeric(parts[3])
    cmd <- paste(parts[11:length(parts)], collapse = " ")

    # Extract model type
    model_type <- if (grepl("sf-ssls-lfo", cmd)) {
      "sf-ssls-lfo"
    } else if (grepl("sf-ssm-log-space", cmd)) {
      "sf-ssm-log-space"
    } else {
      "unknown"
    }

    # Extract output file path
    output_match <- regmatches(cmd, regexpr("output file=[^ ]+", cmd))
    output_file <- if (length(output_match) > 0) {
      sub("output file=", "", output_match)
    } else {
      NA_character_
    }

    list(pid = pid, cpu = cpu, cmd = cmd, model_type = model_type, output_file = output_file)
  })

  results <- Filter(Negate(is.null), results)
  if (length(results) == 0) {
    return(data.frame(
      pid = integer(0),
      ppid = integer(0),
      cpu = numeric(0),
      cmd = character(0),
      model_type = character(0),
      output_file = character(0),
      stringsAsFactors = FALSE
    ))
  }

  # Get parent PIDs
  pids <- sapply(results, `[[`, "pid")
  ppid_cmd <- sprintf("ps -o pid,ppid -p %s --no-headers", paste(pids, collapse = ","))
  ppid_output <- tryCatch({
    system(ppid_cmd, intern = TRUE)
  }, error = function(e) {
    character(0)
  })

  ppid_map <- list()
  for (line in ppid_output) {
    parts <- strsplit(trimws(line), "\\s+")[[1]]
    if (length(parts) >= 2) {
      ppid_map[[parts[1]]] <- as.integer(parts[2])
    }
  }

  # Build result data frame
  df <- data.frame(
    pid = sapply(results, `[[`, "pid"),
    ppid = sapply(results, function(r) ppid_map[[as.character(r$pid)]] %||% NA_integer_),
    cpu = sapply(results, `[[`, "cpu"),
    cmd = sapply(results, `[[`, "cmd"),
    model_type = sapply(results, `[[`, "model_type"),
    output_file = sapply(results, `[[`, "output_file"),
    stringsAsFactors = FALSE
  )

  df
}

#' Get process tree for a PID
#'
#' @param pid Process ID
#' @return Character string describing the process tree
get_process_tree <- function(pid) {
  cmd <- sprintf("pstree -p -s %d 2>/dev/null", pid)
  tryCatch({
    output <- system(cmd, intern = TRUE)
    paste(output, collapse = "\n")
  }, error = function(e) {
    NA_character_
  })
}

#' Find the controller process for a Stan process
#'
#' @param stan_pid PID of a Stan process
#' @return List with controller info or NULL
find_controller <- function(stan_pid) {
  tree <- get_process_tree(stan_pid)
  if (is.na(tree)) return(NULL)

  # Parse the tree to find the R controller
  # Typical pattern: bash(1)---R(CONTROLLER)---R(WORKER)---sf-ssls-lfo(STAN)
  # or: bash(1)---R(CONTROLLER)---sf-ssls-lfo(STAN)

  # Extract all R process PIDs from the tree
  r_pids <- as.integer(regmatches(tree, gregexpr("R\\((\\d+)\\)", tree))[[1]] |>
    sub("R\\(", "", x = _) |>
    sub("\\)", "", x = _))

  if (length(r_pids) == 0) return(NULL)

  # The first R process is typically the controller
  controller_pid <- r_pids[1]

  # Get controller details
  cmd <- sprintf("ps -o pid,ppid,lstart,cmd -p %d --no-headers", controller_pid)
  output <- tryCatch({
    system(cmd, intern = TRUE)
  }, error = function(e) {
    return(NULL)
  })

  if (length(output) == 0) return(NULL)

  # Parse the output (format: PID PPID Mon Day HH:MM:SS YYYY CMD...)
  parts <- strsplit(trimws(output), "\\s+")[[1]]
  if (length(parts) < 7) return(NULL)

  list(
    pid = as.integer(parts[1]),
    ppid = as.integer(parts[2]),
    started = paste(parts[3:7], collapse = " "),
    cmd = paste(parts[8:length(parts)], collapse = " ")
  )
}

#' Find which targets store a controller is managing
#'
#' @param controller_pid PID of the controller process
#' @return Path to the _targets store or NULL
find_store_for_controller <- function(controller_pid) {
  # Check common store locations for matching PID
  store_locations <- c(
    "/mnt/data/analysis-results/karim_naguib/sclc/test4/_targets",
    "/mnt/data/analysis-results/karim_naguib/sclc/test3/_targets",
    "/mnt/data/analysis-results/karim_naguib/sclc/test2/_targets",
    "/mnt/data/analysis-results/karim_naguib/sclc/test/_targets"
  )

  for (store in store_locations) {
    process_file <- file.path(store, "meta", "process")
    if (file.exists(process_file)) {
      content <- tryCatch({
        readLines(process_file)
      }, error = function(e) {
        character(0)
      })

      # Parse pipe-delimited format
      for (line in content) {
        if (startsWith(line, "pid|")) {
          stored_pid <- as.integer(sub("pid\\|", "", line))
          if (!is.na(stored_pid) && stored_pid == controller_pid) {
            return(store)
          }
        }
      }
    }
  }

  NULL
}

#' List all running pipelines with full details
#'
#' @return List of pipeline info objects
list_pipelines <- function() {
  # Find all Stan processes
  stan_procs <- find_stan_processes()

  if (nrow(stan_procs) == 0) {
    return(list(
      success = TRUE,
      pipelines = list(),
      message = "No running pipelines found"
    ))
  }

  # Group by parent PID to find unique controllers
  controllers <- list()

  for (i in seq_len(nrow(stan_procs))) {
    proc <- stan_procs[i, ]

    # Find controller
    controller <- find_controller(proc$pid)
    if (is.null(controller)) next

    controller_key <- as.character(controller$pid)

    if (!controller_key %in% names(controllers)) {
      # Find store
      store <- find_store_for_controller(controller$pid)

      # Derive store name from output path if not found via metadata
      store_name <- if (!is.null(store)) {
        # Extract test name from path like .../sclc/test3/_targets
        m <- regmatches(store, regexpr("test\\d*", store))
        if (length(m) > 0) m else basename(dirname(store))
      } else if (!is.na(proc$output_file)) {
        m <- regmatches(proc$output_file, regexpr("test\\d*", proc$output_file))
        if (length(m) > 0) m else "unknown"
      } else {
        "unknown"
      }

      controllers[[controller_key]] <- list(
        controller_pid = controller$pid,
        controller_started = controller$started,
        store = store,
        store_name = store_name,
        model_type = proc$model_type,
        stan_processes = list(),
        total_cpu = 0
      )
    }

    # Add this Stan process
    controllers[[controller_key]]$stan_processes <- c(
      controllers[[controller_key]]$stan_processes,
      list(list(
        pid = proc$pid,
        cpu = proc$cpu,
        chain_id = regmatches(proc$cmd, regexpr("id=\\d+", proc$cmd)) |>
          sub("id=", "", x = _) |>
          as.integer()
      ))
    )
    controllers[[controller_key]]$total_cpu <- controllers[[controller_key]]$total_cpu + proc$cpu
  }

  list(
    success = TRUE,
    pipelines = unname(controllers),
    count = length(controllers)
  )
}

#' Kill a pipeline by controller PID or store path
#'
#' @param pid Controller PID (optional)
#' @param store Store path (optional, used to find controller if pid not provided)
#' @return Result of the kill operation
kill_pipeline <- function(pid = NULL, store = NULL) {
  # Find controller PID
  controller_pid <- pid

  if (is.null(controller_pid) && !is.null(store)) {
    # Find controller from store metadata
    process_file <- file.path(store, "meta", "process")
    if (file.exists(process_file)) {
      content <- readLines(process_file)
      for (line in content) {
        if (startsWith(line, "pid|")) {
          controller_pid <- as.integer(sub("pid\\|", "", line))
          break
        }
      }
    }
  }

  if (is.null(controller_pid)) {
    return(list(
      success = FALSE,
      error = "Could not determine controller PID. Provide either pid or store path."
    ))
  }

  # Check if process exists
  check_cmd <- sprintf("ps -p %d --no-headers 2>/dev/null", controller_pid)
  exists <- length(system(check_cmd, intern = TRUE)) > 0

  if (!exists) {
    return(list(
      success = FALSE,
      error = sprintf("Process %d does not exist", controller_pid)
    ))
  }

  # Find child Stan processes first (in case they become orphaned)
  stan_procs <- find_stan_processes()
  related_stan <- stan_procs[stan_procs$ppid == controller_pid |
                             sapply(stan_procs$pid, function(p) {
                               ctrl <- find_controller(p)
                               !is.null(ctrl) && ctrl$pid == controller_pid
                             }), ]

  # Kill controller
  kill_result <- system(sprintf("kill %d 2>&1", controller_pid), intern = TRUE)

  # Wait briefly and check for orphaned Stan processes
  Sys.sleep(1)

  orphaned <- list()
  for (pid in related_stan$pid) {
    check <- sprintf("ps -p %d --no-headers 2>/dev/null", pid)
    if (length(system(check, intern = TRUE)) > 0) {
      # Process still running, kill it
      system(sprintf("kill %d 2>&1", pid))
      orphaned <- c(orphaned, pid)
    }
  }

  # Verify everything is dead
  Sys.sleep(1)
  final_check <- sprintf("ps -p %d --no-headers 2>/dev/null", controller_pid)
  controller_dead <- length(system(final_check, intern = TRUE)) == 0

  list(
    success = controller_dead,
    killed_controller = controller_pid,
    killed_orphans = unlist(orphaned),
    message = if (controller_dead) {
      sprintf("Successfully killed pipeline (controller PID %d)", controller_pid)
    } else {
      "Controller may still be running. Use kill -9 if needed."
    }
  )
}
