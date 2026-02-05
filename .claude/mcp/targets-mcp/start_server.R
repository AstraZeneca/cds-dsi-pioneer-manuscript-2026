#!/usr/bin/env Rscript
# Start the targets-mcp plumber server
# Usage: Rscript start_server.R [port]

# Activate renv if it exists
renv_activate <- file.path("/mnt/code", "renv", "activate.R")
if (file.exists(renv_activate)) {
  cat("Activating renv...\n")
  suppressMessages(suppressWarnings(source(renv_activate, local = TRUE)))
}

# Disable conflicted package messages
options(conflicts.policy = list(warn = FALSE))

suppressPackageStartupMessages(library(plumber))

# Get port from args or use default
args <- commandArgs(trailingOnly = TRUE)
port <- if (length(args) > 0) as.integer(args[1]) else 8484

# Get the directory where this script lives
get_script_dir <- function() {
  # Try multiple methods to find the script path
  # Method 1: commandArgs for Rscript
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    script_path <- sub("^--file=", "", file_arg)
    # Handle both absolute and relative paths
    if (!startsWith(script_path, "/")) {
      script_path <- file.path(getwd(), script_path)
    }
    return(dirname(normalizePath(script_path, mustWork = FALSE)))
  }

  # Method 2: sys.frame for source()
  for (i in seq_len(sys.nframe())) {
    ofile <- sys.frame(i)$ofile
    if (!is.null(ofile)) {
      return(dirname(normalizePath(ofile, mustWork = FALSE)))
    }
  }

  # Fallback: current directory
  getwd()
}

script_dir <- get_script_dir()
api_file <- file.path(script_dir, "api.R")

# Set env var so api.R can find sessions.R
Sys.setenv(TARGETS_MCP_DIR = script_dir)

cat(sprintf("Starting targets-mcp server on port %d...\n", port))
cat(sprintf("API file: %s\n", api_file))
cat(sprintf("TARGETS_MCP_DIR: %s\n", script_dir))

# Create and run plumber
pr <- plumber::plumb(api_file)
pr$run(host = "127.0.0.1", port = port)
