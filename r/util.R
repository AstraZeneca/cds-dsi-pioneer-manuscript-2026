# nolint start: object_usage_linter

#' Utility Functions for Stan Model Sampling and Analysis
#'
#' This file contains utility functions for conducting Bayesian analysis of
#' oncology clinical trials, with particular focus on tumor growth modeling
#' and survival analysis using Stan.
#'
#' Key functionality includes:
#' - Stan model sampling with automatic file saving
#' - Kaplan-Meier survival analysis helpers
#' - RECIST response evaluation and trajectory analysis
#' - Prior specification utilities for tumor growth and PFS models
#' - Data preprocessing and visualization helpers
#' - File format handling for rvar objects
#'
#' Dependencies:
#' @importFrom posterior rfun is_rvar
#' @importFrom survival survfit2 Surv
#' @importFrom broom tidy
#' @importFrom purrr map_dfr map2_chr accumulate
#' @importFrom dplyr mutate select bind_cols group_by case_when
#' @importFrom tibble is_tibble as_tibble
#' @importFrom qs2 qs_save qs_read
#' @importFrom scales label_number
#' @importFrom targets tar_combine_raw tar_select_targets tar_format
#'

#' Sample from a Stan model and optionally save the output
#'
#' @param model A Stan model object
#' @param ... Additional arguments passed to the sample method
#' @param output_dir Directory to save output files
#' @param output_basename Base name for output files
#' @param timestamp Boolean, whether to include a timestamp in file names
#' @param no_save Boolean, if TRUE, don't save any output files
#'
#' @return A fitted Stan model object
sample_and_save <- function(
  exe_file,
  ...,
  output_dir,
  output_basename = NULL,
  timestamp = TRUE,
  no_save = FALSE,
  save_profiles = TRUE,
  sampler_fun = c("sample", "pathfinder", "variational")
) {
  sampler_fun <- arg_match(sampler_fun)

  fs::dir_create(output_dir, recurse = TRUE)

  # Load the compiled model from exe_file
  model <- cmdstan_model(exe_file = exe_file)

  # Ensure the compiled Stan executable has execute permissions
  if (fs::file_exists(exe_file)) {
    fs::file_chmod(exe_file, "u+x")
  }

  if (!no_save && !timestamp) {
    # fit <- model$sample(..., output_dir = output_dir, output_basename = output_basename)
    fit <- exec(
      model[[sampler_fun]],
      output_dir = output_dir,
      output_basename = output_basename,
      ...
    )
  } else {
    # fit <- model$sample(...)
    fit <- exec(model[[sampler_fun]], !!!list(...))

    if (!no_save) {
      fit$save_output_files(
        dir = output_dir,
        basename = output_basename,
        random = FALSE,
        timestamp = timestamp
      )
    }
  }

  if (!no_save && save_profiles) {
    tryCatch(
      {
        fit$save_profile_files(
          dir = output_dir,
          basename = output_basename,
          random = FALSE,
          timestamp = timestamp
        )
      },
      error = function(e) {
        warning("Failed to save profile files: ", conditionMessage(e))
        # Optionally return a sentinel or do nothing
        NULL
      }
    )
  }

  return(fit)
}

#' Re-run generated quantities block against existing posterior draws
#'
#' Uses CmdStanR's generate_quantities() to re-execute only the GQ block of a
#' compiled Stan model against an existing fit's parameter draws, without
#' re-running MCMC sampling. This is used when only GQ code has changed.
#'
#' @param exe_file Path to compiled Stan executable (with updated GQ code)
#' @param fit CmdStanMCMC fit object providing existing posterior draws
#' @param data Stan data list
#' @param output_dir Directory to save GQ CSV output files
#' @param parallel_chains Number of chains to run in parallel (default 4)
#' @return CmdStanGQ fit object
generate_quantities_from_fit <- function(
  exe_file,
  fit,
  data,
  output_dir,
  parallel_chains = 4
) {
  model <- cmdstan_model(exe_file = exe_file)

  if (fs::file_exists(exe_file)) {
    fs::file_chmod(exe_file, "u+x")
  }

  fs::dir_create(output_dir, recurse = TRUE)

  model$generate_quantities(
    fitted_params = fit,
    data = data,
    output_dir = output_dir,
    parallel_chains = parallel_chains
  )
}

#' Compute CIF draws from existing fit CSVs by applying the state-3 correction in R
#'
#' TEMPORARY workaround: the full model GQ has an expensive patient-state generation
#' step that makes `generate_quantities_from_fit()` unreliable (37+ min/target, chain
#' crashes). Instead, we read the pre-computed per-patient endpoint arrays directly
#' from the fit CSVs and recompute only the CIF in R.
#'
#' The only correction needed is: for state-3 patients where T_SF < T_dropout
#' (case 2), flip `spop_ms_right_censored` from 1 → 0. The timing (`spop_ms_pfs`)
#' is already correct in the existing draws because
#' `spop_ms_pfs = spop_pfs = min(T_SF, T_dropout) = T_SF` for these patients.
#'
#' @param fit CmdStanMCMC fit object (existing posterior draws)
#' @param stan_data Stan data list (needs n_patients, n_trials, max_all_t,
#'   trial_patient_pos)
#' @return posterior::draws_df with spop/sample_cif_01/02/03[trial,time] variables
compute_cif_from_draws <- function(fit, stan_data) {
  n_patients <- stan_data$n_patients
  n_trials   <- stan_data$n_trials
  # max_all_t and trial_patient_pos are Stan transformed_data — reconstruct here.
  # Stan: max_all_t = max(max(t_patient_visits) + 1, extend_max_all_t)
  max_all_t  <- max(max(stan_data$t_patient_visits) + 1L, stan_data$extend_max_all_t)
  T_len      <- max_all_t + 1L
  n_trial_patients <- tabulate(stan_data$patient_trial, nbins = n_trials)
  tpp <- c(1L, cumsum(n_trial_patients) + 1L)  # length n_trials + 1

  pvars <- c(
    "spop_ms_pfs", "spop_ms_right_censored",
    "spop_right_censored", "spop_target_right_censored", "spop_target_pfs",
    "spop_pfs", "spop_os", "spop_os_censored",
    "sample_ms_pfs", "sample_ms_right_censored",
    "sample_os", "sample_os_censored"
  )

  d_mat <- posterior::as_draws_matrix(fit$draws(variables = pvars))
  n_draws <- nrow(d_mat)

  # Pre-build column index maps (1..n_patients per variable)
  col_idx <- function(varname) {
    base::match(paste0(varname, "[", seq_len(n_patients), "]"), colnames(d_mat))
  }
  ci <- list(
    sms_pfs   = col_idx("spop_ms_pfs"),
    sms_rc    = col_idx("spop_ms_right_censored"),
    s_rc      = col_idx("spop_right_censored"),
    st_rc     = col_idx("spop_target_right_censored"),
    st_pfs    = col_idx("spop_target_pfs"),
    s_pfs     = col_idx("spop_pfs"),
    s_os      = col_idx("spop_os"),
    s_osc     = col_idx("spop_os_censored"),
    sam_pfs   = col_idx("sample_ms_pfs"),
    sam_rc    = col_idx("sample_ms_right_censored"),
    sam_os    = col_idx("sample_os"),
    sam_osc   = col_idx("sample_os_censored")
  )

  # Output CIF arrays [n_draws, n_trials, T_len]
  cif_arrs <- list(
    spop_cif_01   = array(0, c(n_draws, n_trials, T_len)),
    spop_cif_02   = array(0, c(n_draws, n_trials, T_len)),
    spop_cif_03   = array(0, c(n_draws, n_trials, T_len)),
    sample_cif_01 = array(0, c(n_draws, n_trials, T_len)),
    sample_cif_02 = array(0, c(n_draws, n_trials, T_len)),
    sample_cif_03 = array(0, c(n_draws, n_trials, T_len))
  )

  for (i in seq_len(n_draws)) {
    row <- d_mat[i, ]
    get_int <- function(idx) as.integer(round(row[idx]))

    sms_pfs   <- get_int(ci$sms_pfs)
    sms_rc    <- get_int(ci$sms_rc)
    s_rc      <- get_int(ci$s_rc)
    st_rc     <- get_int(ci$st_rc)
    st_pfs    <- get_int(ci$st_pfs)
    s_pfs     <- get_int(ci$s_pfs)
    s_os      <- get_int(ci$s_os)
    s_osc     <- get_int(ci$s_osc)
    sam_pfs   <- get_int(ci$sam_pfs)
    sam_rc    <- get_int(ci$sam_rc)
    sam_os    <- get_int(ci$sam_os)
    sam_osc   <- get_int(ci$sam_osc)

    # Mirror the Stan GQ correction (_endpoints_generated_quantities.stan:336-347).
    # This overwrites spop_ms_pfs/spop_ms_right_censored before the CIF module runs.
    # Case 1: PFS event (combined) — use combined PFS time (min of target + ms)
    case1 <- s_rc == 0L
    sms_pfs[case1] <- s_pfs[case1]
    sms_rc[case1]  <- 0L
    # Case 2: dropout where SF progression precedes dropout — reclassify as 0→1
    case2 <- s_rc == 1L & st_rc == 0L & st_pfs <= sms_pfs
    sms_pfs[case2] <- st_pfs[case2]
    sms_rc[case2]  <- 0L
    # Case 3 (else): leave unchanged — dropout or fully censored
    # Note: Stan GQ applies NO correction to sample_* variables.

    for (s in seq_len(n_trials)) {
      tr   <- seq(tpp[s], tpp[s + 1L] - 1L)
      n_tr <- length(tr)
      if (n_tr == 0L) next

      # ── spop (unconditional posterior predictive) ──────────────────────────
      pfs_e <- sms_rc[tr] == 0L
      dd    <- pfs_e & s_osc[tr] == 0L & sms_pfs[tr] == s_os[tr]
      prog  <- pfs_e & !dd
      drop_ <- !pfs_e & sms_pfs[tr] <= max_all_t

      cif_arrs$spop_cif_01[i, s, ] <- cumsum(tabulate(sms_pfs[tr][prog],  nbins = T_len)) / n_tr
      cif_arrs$spop_cif_02[i, s, ] <- cumsum(tabulate(sms_pfs[tr][dd],    nbins = T_len)) / n_tr
      cif_arrs$spop_cif_03[i, s, ] <- cumsum(tabulate(sms_pfs[tr][drop_], nbins = T_len)) / n_tr

      # ── sample (conditional on observed data, no correction — matches Stan GQ)
      pfs_e_s <- sam_rc[tr] == 0L
      dd_s    <- pfs_e_s & sam_osc[tr] == 0L & sam_pfs[tr] == sam_os[tr]
      prog_s  <- pfs_e_s & !dd_s
      drop_s  <- !pfs_e_s & sam_pfs[tr] <= max_all_t

      cif_arrs$sample_cif_01[i, s, ] <- cumsum(tabulate(sam_pfs[tr][prog_s], nbins = T_len)) / n_tr
      cif_arrs$sample_cif_02[i, s, ] <- cumsum(tabulate(sam_pfs[tr][dd_s],   nbins = T_len)) / n_tr
      cif_arrs$sample_cif_03[i, s, ] <- cumsum(tabulate(sam_pfs[tr][drop_s], nbins = T_len)) / n_tr
    }
  }

  # Return a draws_matrix with column names matching Stan convention:
  # spop_cif_01[trial,time], spop_cif_02[trial,time], etc.
  # We build a plain numeric matrix first, then convert via as_draws_matrix().
  # DO NOT return a tibble/data.frame — posterior::as_draws_matrix() on a
  # data.frame routes through as_draws_df(), which reconsolidates bracket-named
  # columns into list-typed array variables, breaking downstream as.numeric().
  n_cif_cols <- length(cif_arrs) * n_trials * T_len
  result_mat <- matrix(0, nrow = n_draws, ncol = n_cif_cols)
  col_names  <- character(n_cif_cols)
  k <- 0L
  for (nm in names(cif_arrs)) {
    arr <- cif_arrs[[nm]]
    for (s in seq_len(n_trials)) {
      for (t in seq_len(T_len)) {
        k <- k + 1L
        result_mat[, k] <- arr[, s, t]
        col_names[k] <- paste0(nm, "[", s, ",", t, "]")
      }
    }
  }
  colnames(result_mat) <- col_names
  posterior::as_draws_matrix(result_mat)
}

#' Select draws from CmdStanR fit using tidyselect patterns
#'
#' Uses cmdstanr::read_cmdstan_csv with variable selection to read only the
#' needed variables directly from CSV files. This bypasses any caching in
#' the fit object and ensures minimal memory usage.
#'
#' @param fit A CmdStanMCMC fit object
#' @param ... Tidyselect expressions to filter variables (e.g., ends_with("_pop"))
#' @return A draws_array object with selected variables
#' @export
select_draws <- function(fit, ...) {
  # Get all variable names from the fit (base names without indices)
  all_vars <- fit$metadata()$stan_variables

  if (...length() == 0) {
    # No selection - read all variables
    selected_vars <- NULL
  } else {
    # Use substitute to capture the unevaluated selection expression
    # This ensures tidyselect helpers are found in eval_select's environment
    selection <- substitute(c(...))
    selected_idx <- tidyselect::eval_select(
      selection,
      data = rlang::set_names(all_vars, all_vars)
    )
    selected_vars <- all_vars[selected_idx]

    if (length(selected_vars) == 0) {
      stop("No variables matched the selection criteria")
    }
  }

  # Read directly from CSV files with variable selection
  # This bypasses fit object caching for better memory efficiency
  cmdstanr::read_cmdstan_csv(
    files = fit$output_files(),
    variables = selected_vars
  )$post_warmup_draws
}

#' Extract and compile decorated Stan functions
#'
#' @param stan_file Path to .stan file with decorated functions
#' @param includes Optional vector of #include statements
#' @return cmdstanr model object with exported functions
export_stan_functions <- function(stan_file, includes = NULL) {
  # Read the Stan file
  content <- readLines(stan_file) %>% paste(collapse = "\n")

  # Extract functions marked with @stan_export anywhere in their documentation
  # Find @stan_export markers and extract the following function

  # Split content into lines for easier processing
  lines <- strsplit(content, "\n")[[1]]
  export_lines <- which(str_detect(lines, "@stan_export"))

  if (length(export_lines) == 0) {
    stop("No functions marked with @stan_export found in ", stan_file)
  }

  exported_functions <- c()

  for (export_line in export_lines) {
    # Find the start of the function signature after @stan_export
    func_start <- NULL

    # Look for the start of a function signature (could be multi-line)
    for (i in (export_line + 1):length(lines)) {
      line <- trimws(lines[i])

      # Skip empty lines and comment lines
      if (
        line == "" ||
          str_detect(line, "^\\s*//") ||
          str_detect(line, "^\\s*/\\*")
      ) {
        next
      }

      # Check if this line starts a function (return type or tuple)
      if (
        str_detect(
          line,
          "^\\s*(?:real|int|vector|matrix|row_vector|array|void|tuple)"
        )
      ) {
        func_start <- i
        break
      }
    }

    if (is.null(func_start)) {
      next
    }

    # Find the complete function by tracking braces
    brace_count <- 0
    func_end <- NULL
    found_opening_brace <- FALSE

    for (i in func_start:length(lines)) {
      line <- lines[i]

      # Count opening and closing braces
      open_braces <- str_count(line, "\\{")
      close_braces <- str_count(line, "\\}")

      if (open_braces > 0) {
        found_opening_brace <- TRUE
      }

      brace_count <- brace_count + open_braces - close_braces

      if (found_opening_brace && brace_count == 0) {
        func_end <- i
        break
      }
    }

    if (!is.null(func_end)) {
      func_text <- paste(lines[func_start:func_end], collapse = "\n")
      exported_functions <- c(exported_functions, func_text)
    }
  }

  if (length(exported_functions) == 0) {
    stop("No valid functions found after @stan_export markers in ", stan_file)
  }

  functions_code <- exported_functions

  # Create Stan program with functions block
  stan_program <- paste0(
    "functions {\n",
    if (!is.null(includes)) paste0("  ", includes, collapse = "\n"),
    "\n",
    paste(functions_code, collapse = "\n\n"),
    "\n",
    "}\n\n",
    "data {}\n",
    "parameters {}\n",
    "model {}\n"
  )

  # Compile directly from string and expose functions
  model <- cmdstan_model(
    stan_file = write_stan_file(stan_program),
    force_recompile = TRUE
  )
  model$expose_functions()
  return(model)
}

# Compute hash of Stan model source files
# Returns a hash string that changes when any source file content changes
compute_stan_source_hash <- function(model_file, include_files = NULL) {
  all_source_files <- c(model_file, include_files)
  all_source_files |>
    sort() |>
    map(read_lines) |>
    digest::digest(algo = "md5")
}

build_model <- function(
  model_file,
  include_files = NULL,
  dir = NULL,
  compile_cores = parallel::detectCores(),
  ...
) {
  # Compute hash of all source file contents to detect changes
  source_hash <- compute_stan_source_hash(model_file, include_files)

  # Determine expected executable path
  model_name <- tools::file_path_sans_ext(fs::path_file(model_file))
  exe_dir <- dir %||% fs::path_dir(model_file)
  exe_path <- fs::path(exe_dir, model_name)
  hash_file <- str_c(exe_path, ".source_hash")

  # Compare with stored hash - delete binary if sources changed
  skip_compile <- FALSE
  if (fs::file_exists(exe_path)) {
    if (fs::file_exists(hash_file)) {
      stored_hash <- read_lines(hash_file, n_max = 1)
      if (!identical(stored_hash, source_hash)) {
        message("Source hash changed - forcing recompilation")
        fs::file_delete(exe_path)
      } else {
        message("Source hash matches - skipping recompilation")
        skip_compile <- TRUE
        # Ensure executable has proper permissions when reusing from shared storage
        # Fixes "Permission denied" errors across Domino jobs/workspaces
        Sys.chmod(exe_path, mode = "0755")
      }
    } else {
      message("Hash file not found - forcing recompilation")
      fs::file_delete(exe_path)
    }
  }

  # Set MAKEFLAGS for parallel compilation
  withr::local_envvar(MAKEFLAGS = str_c("-j", compile_cores))

  model <- cmdstan_model(
    model_file,
    cpp_options = lst(
      stan_threads = TRUE,
      "CXXFLAGS += -O3",
      "CXXFLAGS += --march=native" # Optimize for local CPU architecture
    ),
    # stanc_options = list("O1"),      # Stan compiler optimizations
    dir = dir,
    force_recompile = !skip_compile,
    ...
  )

  # Save source hash for future comparisons
  write_lines(source_hash, hash_file)

  # Track the executable by including its hash in the return value
  exe_path <- model$exe_file()

  # Ensure executable has proper permissions for shared storage (Domino)
  # Fixes "Permission denied" errors across jobs/workspaces
  Sys.chmod(exe_path, mode = "0755")

  exe_hash <- digest::digest(file = exe_path, algo = "md5")

  # Store the hash as an attribute so targets tracks it
  attr(model, "exe_hash") <- exe_hash

  return(model)
}

remove_incomplete_cases <- function(data, incomplete) {
  if (is_empty(incomplete)) {
    return(data)
  } else {
    return(slice(data, -incomplete))
  }
}

get_conditioning_subgroups <- function(data, cond, other_cond) {
  map(cond$cond_group_expr, \(x) {
    transmute(data, cond = !!x & !!other_cond) |> pull(cond) |> which()
  })
}

#' Convert Kaplan-Meier estimates to a tibble (data frame) format
#'
#' @param trt_data Analysis data
#' @param key Identifier for the data group (e.g., treatment arm)
#' @param pfs_var Name of variable were PFS is stored in the data
#'
#' @return tibble object with Kaplan-Meier results.
km_to_tibble <- function(
  trt_data,
  key,
  pfs_sym,
  censored_sym,
  probs = c(0.5, 0.8)
) {
  survfit_objs <- rlang::inject(
    lst(
      lb = ggsurvfit::survfit2(
        Surv(!!pfs_sym + 1 - !!censored_sym, 1 - !!censored_sym) ~ 1,
        trt_data
      ),
      ub = ggsurvfit::survfit2(
        Surv(
          !!pfs_sym + interval_censored + 1 - !!censored_sym,
          1 - !!censored_sym
        ) ~ 1,
        trt_data
      ),
    )
  )

  survfit_objs |>
    imap_dfr(\(fit, btype) {
      # Get tidy KM estimates
      km_data <- broom::tidy(fit) |>
        select(t = time, s = estimate, n = n.risk, c = n.censor, e = n.event)

      # Get quantiles with confidence intervals from survfit object
      quant_result <- quantile(fit, probs = probs)
      quantiles <- tibble(
        quantile = probs,
        pfs = quant_result$quantile,
        pfs_lower = quant_result$lower,
        pfs_upper = quant_result$upper
      )

      tibble(
        btype = btype,
        km_data = list(km_data),
        quantiles = list(quantiles)
      )
    }) |>
    bind_cols(key)
}

#' Extract Kaplan-Meier estimates with median PFS
#'
#' @param analysis_data The full analysis dataset
#' @param pfs_var Name of the PFS time variable
#' @param censored_var Name of the censoring indicator variable
#' @param ... Additional grouping variables (e.g., treatment arm, biomarker subgroups)
#'
#' @return A tibble with KM estimates and median PFS (with confidence intervals) for each group
get_km_res <- function(
  analysis_data,
  pfs_var,
  censored_var,
  ...,
  probs = c(0.5, 0.8)
) {
  pfs_sym <- rlang::ensym(pfs_var)
  censored_sym <- rlang::ensym(censored_var)

  analysis_data |>
    group_by(trial, ...) |>
    group_map(
      \(trt_data, key) {
        km_to_tibble(trt_data, key, pfs_sym, censored_sym, probs = probs)
      },
      .keep = TRUE
    ) |>
    bind_rows()
}

#' Adjust PFS Data Based on Calendar Day Cutoff
#'
#' This function artificially cuts off analysis data at a specific calendar day,
#' adjusting PFS times and censoring indicators as if the data collection had
#' stopped at that cutoff date. This is useful for retrospective analyses and
#' simulating data availability at earlier timepoints.
#'
#' @param analysis_data A tibble containing the analysis dataset with nested visit_data
#' @param cutoff_calendar_day The calendar day at which to cut off the data
#' @param pfs_var Name of the PFS variable (unquoted, default: pfs)
#' @param week_var Name of the week variable in visit_data (unquoted, default: week)
#' @param calendar_day_var Name of the calendar day variable in visit_data (unquoted, default: calendar_day)
#' @param require_post_baseline Logical; if TRUE, excludes patients who have no post-baseline
#'   visits (week > 0) at the cutoff. This filters out retroactively-added patients who weren't
#'   in the original DCO. Default is FALSE for backwards compatibility.
#'
#' @return A modified tibble with adjusted PFS, right_censored, and interval_censored values
#'
#' @details
#' For each patient, the function:
#' 1. Filters visits to only those occurring on or before the cutoff calendar day
#' 2. Optionally removes patients with no post-baseline visits (if require_post_baseline = TRUE)
#' 3. If the patient's original PFS event occurred after the cutoff, sets them as right-censored
#'    at the time of their last visit before/at the cutoff
#' 4. Recalculates interval censoring based on the gap between the last observed visit
#'    and the cutoff date
#' 5. Preserves the original PFS if it occurred before the cutoff
#'
#' @examples
#' \dontrun{
#' # Simulate a data cutoff 6 months into the trial
#' cutoff_data <- apply_calendar_cutoff(
#'   analysis_data,
#'   cutoff_calendar_day = 180
#' )
#'
#' # Use with custom variable names
#' cutoff_data <- apply_calendar_cutoff(
#'   analysis_data,
#'   cutoff_calendar_day = 180,
#'   pfs_var = my_pfs,
#'   week_var = my_week
#' )
#' }
apply_calendar_cutoff <- function(
  analysis_data,
  cutoff_calendar_day,
  pfs_var = pfs,
  week_var = week,
  calendar_day_var = visit_calendar_day,
  require_post_baseline = FALSE
) {
  analysis_data |>
    mutate(
      # Filter visit data to only include visits at or before cutoff
      visit_data = map(visit_data, \(vd) {
        vd |> filter({{ calendar_day_var }} <= cutoff_calendar_day)
      })
    ) |>
    # Remove patients who have no visits before/at the cutoff
    filter(map_int(visit_data, nrow) > 0) |>
    # Optionally require at least one post-baseline visit (week > 0)
    # This filters out patients who only have baseline/screening data,

    # which typically indicates retroactively-added patients not in original DCO
    filter(
      !require_post_baseline |
        map_lgl(visit_data, \(vd) any(pull(vd, {{ week_var }}) > 0))
    ) |>
    mutate(
      # Determine last observed week before cutoff
      last_obs_week = map_int(visit_data, \(vd) {
        vd |> pull({{ week_var }}) |> max(na.rm = TRUE)
      }),

      # Determine last observed calendar day before cutoff
      last_obs_calendar_day = map_int(visit_data, \(vd) {
        vd |> pull({{ calendar_day_var }}) |> max(na.rm = TRUE)
      }),

      # Adjust PFS: if event was after cutoff, censor at last observation
      cutoff_pfs = if_else(
        last_obs_week < {{ pfs_var }},
        last_obs_week,
        {{ pfs_var }}
      ),

      # Adjust right censoring: becomes censored if original event was after cutoff
      cutoff_right_censored = if_else(
        last_obs_week < {{ pfs_var }},
        1L,
        right_censored
      ),

      # Adjust interval censoring: calculate gap from last visit to cutoff
      # Only applies if we're now right-censored due to cutoff
      cutoff_interval_censored = if_else(
        last_obs_week < {{ pfs_var }},
        pmax(
          0L,
          as.integer(cutoff_calendar_day - last_obs_calendar_day) %/% 7L
        ),
        if_else(
          cutoff_right_censored == 1L,
          pmax(
            0L,
            as.integer(cutoff_calendar_day - last_obs_calendar_day) %/% 7L
          ),
          interval_censored
        )
      )
    ) |>
    # Replace original variables with cutoff versions
    mutate(
      {{ pfs_var }} := cutoff_pfs,
      right_censored = cutoff_right_censored,
      interval_censored = cutoff_interval_censored
    ) |>
    # Clean up temporary columns
    select(!c(starts_with("cutoff_"), last_obs_week, last_obs_calendar_day))
}

add_confirmed_resp_priors <- function(stan_data, priors) {
  stan_data |>
    list_assign(!!!priors)
}

add_pfs_crcr_priors <- function(
  stan_data,
  crcr_priors,
  tumor_priors,
  pfs_priors
) {
  add_confirmed_resp_priors(stan_data, crcr_priors) |>
    list_assign(!!!tumor_priors, !!!pfs_priors)
}

#' Generate a histogram of time-to-events for a single draw
#'
#' This function creates a histogram of time-to-event data for a single draw from a
#' posterior distribution. It uses R's base hist() function but returns only the counts,
#' not the full histogram object.
#'
#' @param pred A numeric vector of predicted time-to-event values
#' @param breaks A numeric vector specifying the breakpoints between histogram cells
#' @param ... Additional arguments passed to hist()
#'
#' @return A numeric vector of counts for each histogram bin
#'
sample_hist <- function(pred, breaks, freq = TRUE, ...) {
  # hist() is a base R function to generate histograms from data and provided breaks.
  hist(
    pmax(pmin(pred, max(breaks)), min(breaks)),
    breaks = breaks,
    plot = FALSE,
    ...
  )[[if (freq) "counts" else "density"]]
}

# This function is used to treated_pfs_analysis_dataallow us to generate a distribution of histograms
rvar_sample_hist <- posterior::rfun(sample_hist, rvar_dots = FALSE)
rvar_weighted_mean <- posterior::rfun(weighted.mean, rvar_args = "x")
rvar_plogis <- posterior::rfun(plogis, rvar_args = "q")

#' Name coefficient indices with meaningful labels
#'
#' This function takes a data frame with coefficient indices and adds meaningful labels
#' to these coefficients based on their index and the provided Stan data. It also
#' optionally adds trial labels.
#'
#' @param data A data frame containing the coefficient indices to be named
#' @param coef_idx_col The name of the column in 'data' that contains the coefficient indices
#' @param trial_col The name of the column in 'data' that contains trial identifiers (optional)
#' @param stan_data A list containing Stan data, including 'covar_design_matrix' and 'patient_trial'
#'
#' @return A modified data frame with additional columns:
#'   - 'covar': A factor column with meaningful names for each coefficient
#'   - 'trial': A factor column with trial labels (if trial_col is provided)
name_coef_indices <- function(data, coef_idx_col, trial_col, stan_data) {
  data |>
    mutate(
      covar = case_when(
        {{ coef_idx_col }} == 1 ~ "baseline sum of tumor sizes",
        {{ coef_idx_col }} == 2 ~ "first post-treatment sum of tumor sizes",
        {{ coef_idx_col }} - 2 <= ncol(stan_data$covar_design_matrix) ~
          colnames(stan_data$covar_design_matrix)[pmax(
            1,
            {{ coef_idx_col }} - 2
          )] |>
          str_replace(r"{factor\((.+),\sordered\s=\sFALSE\)}", "\\1 "),
        TRUE ~ "confirmed response"
      ) |>
        as_factor(),
      trial = if (!is_null(trial_col)) {
        factor({{ trial_col }}, labels = levels(stan_data$patient_trial))
      },
    )
}

get_fake_stan_data_list <- function(prior_res, origin_stan_data, n = 5) {
  get_all_confirmed_response(prior_res) |>
    select(starts_with("rep_")) |>
    unnest_rvars() |>
    filter(.draw <= n) |>
    rename_with(\(n) str_remove(n, "^rep_")) |>
    select(!c(.chain, .iteration)) |>
    group_by(.draw) |>
    group_map(\(d, k, ...) {
      list_assign(
        origin_stan_data,
        !!!d,
        draw = first(k$.draw),
        confirmed_response_interval_censored = rep(0, nrow(d))
      )
    })
}

#' Determine RECIST 1.1 Response
#'
#' This function calculates the RECIST 1.1 response category based on measurements
#' of target lesions, and optionally considers non-target lesions and new lesions.
#'
#' @param baseline_sld Baseline sum of longest diameter
#' @param current_sld Current sum of longest diameter
#' @param nadir_sld The smallest sum of measurements observed so far (default is NULL, which means baseline is used as nadir)
#' @param include_non_target Logical indicating whether to consider non-target lesions (default is FALSE)
#' @param non_target_response A character string: "CR", "SD", or "PD" (only used if include_non_target = TRUE)
#' @param new_lesions Logical indicating whether new lesions have appeared (default is FALSE)
#'
#' @return A character string indicating the RECIST 1.1 response category
#'
determine_recist_response <- function(
  baseline_sld,
  current_sld,
  nadir_sld = NULL,
  include_non_target = FALSE,
  non_target_response = "NON-CR/NON-PD",
  new_lesions = FALSE
) {
  # Input validation
  if (
    !is.numeric(baseline_sld) ||
      !is.numeric(current_sld) ||
      baseline_sld < 0 ||
      current_sld < 0
  ) {
    stop("baseline_sld and current_sld must be non-negative numeric values")
  }
  if (!is.null(nadir_sld) && (!is.numeric(nadir_sld) || nadir_sld < 0)) {
    stop("nadir_sld must be a non-negative numeric value or NULL")
  }
  if (!is.logical(include_non_target) || !is.logical(new_lesions)) {
    stop("include_non_target and new_lesions must be logical values")
  }
  if (
    !is.na(non_target_response) &&
      !non_target_response %in% c("CR", "NON-CR/NON-PD", "PD", "NE")
  ) {
    stop(
      "non_target_response must be 'CR', 'NON-CR/NON-PD', 'PD', or NA, got: ",
      non_target_response
    )
  }

  # If nadir sum not provided, use baseline as nadir
  if (is_null(nadir_sld)) {
    nadir_sld <- min(baseline_sld, current_sld)
  } else {
    nadir_sld <- min(nadir_sld, current_sld) # Update nadir if current sum is smaller
  }

  change_from_baseline <- (current_sld - baseline_sld) / baseline_sld
  absolute_diff_from_nadir <- current_sld - nadir_sld
  change_from_nadir <- ifelse(
    nadir_sld > 0,
    absolute_diff_from_nadir / nadir_sld,
    NA
  )

  # RECIST 1.1: Order of checks: CR, PD, PR, SD
  if (!include_non_target && !new_lesions) {
    case_when(
      current_sld == 0 ~ "CR",
      # PD: ≥20% increase from nadir AND ≥5mm absolute increase
      nadir_sld > 0 &
        current_sld > nadir_sld &
        change_from_nadir >= 0.2 &
        absolute_diff_from_nadir >= 5 ~ "PD",
      change_from_baseline <= -0.3 ~ "PR",
      .default = "SD"
    )
  } else {
    case_when(
      # If including non-target lesions or new lesions
      (nadir_sld > 0 &
        current_sld > nadir_sld &
        change_from_nadir >= 0.2 &
        absolute_diff_from_nadir >= 5) ||
        (!is.na(non_target_response) && non_target_response == "PD") ||
        new_lesions ~ "PD",
      current_sld == 0 &&
        (non_target_response == "CR" || is.na(non_target_response)) ~ "CR",
      change_from_baseline <= -0.3 ~ "PR",
      .default = "SD"
    )
  }
}

determine_trajectory_recist_target_response <- function(sld) {
  nadir <- accumulate(sld, min)

  map2_chr(sld[-1], nadir[-1], \(curr_sld, curr_nadir) {
    determine_recist_response(first(sld), curr_sld, curr_nadir)
  })
}

weeks_to_months <- function(weeks) weeks * 7 * 12 / 365.25
label_weeks_to_months <- scales::label_number(scale = weeks_to_months(1))
months_to_weeks <- function(months) months / weeks_to_months(1)

lognormal_sd <- function(mu = 0, sigma) {
  # Calculate standard deviation of lognormal variable X
  # where log(X) ~ N(mu, sigma)
  #
  # Args:
  #   mu: mean parameter of the normal distribution in log space
  #   sigma: standard deviation parameter of the normal distribution in log space
  #
  # Returns:
  #   standard deviation of the lognormal random variable X

  sqrt((exp(sigma^2) - 1) * exp(2 * mu + sigma^2))
}

tar_bind_rows <- function(target_name, mapped, names, ...) {
  tar_combine_raw(
    deparse(substitute(target_name)),
    tar_select_targets(mapped, names),
    command = expression(bind_rows(!!!.x)),
    ...
  )
}

#' Safe QS2 Format for RVar Objects in Targets Pipeline
#'
#' A custom targets format that safely handles tibbles containing posterior::rvar
#' objects by removing problematic cache attributes before serialization and
#' ensuring proper tibble conversion on read.
#'
#' @details
#' This format addresses issues with serializing rvar objects that contain
#' cached attributes which can cause problems during the save/load process.
#' The format:
#' 1. Detects tibbles containing rvar columns
#' 2. Removes the "cache" attribute from rvar objects before saving
#' 3. Uses qs2 for efficient serialization
#' 4. Ensures objects are returned as tibbles on read
#'
#' The marshal/unmarshal functions are pass-through (identity functions)
#' since the main processing happens in write/read.
#'
#' @section Usage:
#' Use this format in targets pipelines when working with posterior samples
#' stored as rvar objects:
#' ```
#' tar_target(
#'   name = my_posterior_data,
#'   command = analyze_posterior(),
#'   format = rvar_safe_qs2_format
#' )
#' ```
#'
#' @section Performance:
#' - Uses qs2 for fast serialization of large objects
#' - Minimal overhead for non-rvar objects
#' - Only processes rvar columns when detected
#'
#' @return A targets format object with custom read/write methods
#'
#' @seealso
#' - [targets::tar_format()] for creating custom formats
#' - [posterior::rvar()] for random variable objects
#' - [qs2::qs_save()] and [qs2::qs_read()] for serialization
#'
#' @examples
#' \dontrun{
#' # In a _targets.R file
#' library(targets)
#' library(posterior)
#'
#' tar_pipeline(
#'   tar_target(
#'     posterior_results,
#'     my_stan_analysis(),
#'     format = rvar_safe_qs2_format
#'   )
#' )
#' }
rvar_safe_qs2_format <- tar_format(
  write = function(object, path) {
    if (
      tibble::is_tibble(object) &&
        any(purrr::map_lgl(object, posterior::is_rvar))
    ) {
      object <- as.data.frame(object) |>
        dplyr::mutate(across(where(posterior::is_rvar), \(r) {
          attr(r, "cache") <- NULL
          r
        }))
    }

    qs2::qs_save(object, path)
  },

  read = function(path) {
    object <- qs2::qs_read(path)

    if (is.data.frame(object) && !tibble::is_tibble(object)) {
      object <- tibble::as_tibble(object)
    }

    return(object)
  }
)

cmdstanr_format <- tar_format(
  read = function(path) {
    cmdstanr::as_cmdstan_fit(
      readr::read_rds(path)$csv_files,
      check_diagnostics = FALSE,
      format = "draws_list"
    )
    # readr::read_rds(path)$file |>
    #   readr::read_rds()
  },
  write = function(object, path) {
    # obj_file <- stringr::str_c(path, "_cmdstanr_object.rds")
    # object$save_object(path)

    # Calculate hash of all CSV files combined
    csv_files <- object$output_files()
    # csv_hash <- digest::digest(purrr::map(csv_files, \(f) digest::digest(file = f)), algo = "xxhash64")  # Fast hash algorithm
    csv_hash <- purrr::map(csv_files, \(f) {
      digest::digest(file = f, algo = "xxhash64")
    }) # Fast hash algorithm

    # Save both fit object and hash
    readr::write_rds(
      tibble::lst(
        # fit = object,
        # file = obj_file,
        csv_files,
        csv_hash # This changes when CSV content changes
      ),
      path
    )
  }
)

tar_cmdstan_sample <- function(
  name,
  model,
  stan_data,
  init_factory = \(...) \(...) NULL,
  ...
) {
  tar_target(
    name,
    sample_and_save(
      model,
      stan_data,
      init = init_factory(stan_data),
      timestamp = FALSE,
      format = cmdstanr_format,
      ...
    ),
  )
}

determine_visit_data_response <- function(visit_data) {
  visit_data |>
    group_by(usubjid) |>
    mutate(
      det_response = c(
        NA,
        determine_trajectory_recist_target_response(mmsumdiam)
      )
    ) |>
    ungroup() |>
    mutate(
      across(c(response, det_response), \(r) {
        ordered(r, levels = c("CR", "PR", "SD", "PD"))
      }),
    )
}

determine_pfs <- function(visit_data, pfs_confirm_visits = 1) {
  visit_data |>
    summarize(
      det_pfs_idx = find_consecutive(det_response, "PD", pfs_confirm_visits), # |> coalesce(max(week)),
      det_pfs = case_when(
        is.na(det_pfs_idx) ~ max(week),
        det_pfs_idx > 1 ~ max(1, week[det_pfs_idx - 1]),
        TRUE ~ 1
      ),
      det_right_censored = is.na(det_pfs_idx),
      det_interval_censored = if_else(
        det_right_censored,
        0,
        coalesce(week[det_pfs_idx] - det_pfs - 1, 0)
      )
    )
}

# Helper function to parse parameter specifications and find matching columns
get_param_names_from_dots <- function(dots, all_cols) {
  if (length(dots) == 0) {
    # If no parameters specified, return all non-metadata columns
    return(setdiff(all_cols, c(".chain", ".iteration", ".draw")))
  }

  # Extract the base parameter names from the expressions
  param_names <- character()

  param_names <- map_chr(dots, function(dot) {
    expr_str <- rlang::as_label(dot)
    # Extract base name (e.g., "beta" from "beta[i]")
    base_name <- stringr::str_extract(expr_str, "^[^\\[]+")
    # Find all columns that match this base name
    str_subset(all_cols, str_glue(r"{^{base_name}(\[|$)}"))
  })

  return(unique(param_names))
}

get_draws <- function(fit, ..., recover_data = NULL) {
  d <- enquos(...) |>
    map_chr(as_label) |>
    str_extract(r"{^[^\[]+}") |>
    fit$draws()

  if (!is_null(recover_data)) {
    d |> recover_types(recover_data)
  } else {
    d
  }
}

lite_spread_rvars <- function(fit, ..., ndraws = NULL, recover_data = NULL) {
  get_draws(fit, ..., recover_data = recover_data) |>
    tidybayes::spread_rvars(..., ndraws = ndraws)
}

lite_gather_rvars <- function(
  fit,
  ...,
  ndraws = NULL,
  recover_data = NULL,
  calc_rhat = FALSE,
  calc_ess = FALSE
) {
  d <- get_draws(fit, ..., recover_data = recover_data) |>
    tidybayes::gather_rvars(..., ndraws = ndraws)

  if (calc_rhat) {
    d <- d |> mutate(rh = posterior::rhat(.value))
  }

  if (calc_ess) {
    d <- d |>
      mutate(
        ess_b = posterior::ess_bulk(.value),
        ess_t = posterior::ess_tail(.value)
      )
  }

  return(d)
}

find_consecutive <- function(vec, x, n = 1) {
  vctrs::vec_unrep(vec) |>
    mutate(
      end_index = cumsum(times),
      start_index = end_index - times + 1
    ) |>
    filter(key == x, times >= n) |>
    pull(start_index) |>
    first() %||%
    NA
}

find_stan_includes <- function(stan_file, base_dir = NULL) {
  # Set base directory - find the stan root directory
  # Walk up from the main file to find a directory named "stan"
  if (is.null(base_dir)) {
    abs_stan_file <- normalizePath(stan_file, mustWork = TRUE)
    current_dir <- dirname(abs_stan_file)

    # Walk up to find the "stan" directory
    while (current_dir != dirname(current_dir)) {
      # Not at filesystem root
      if (basename(current_dir) == "stan") {
        base_dir <- current_dir
        break
      }
      current_dir <- dirname(current_dir)
    }

    # If we didn't find a "stan" directory, fall back to the file's directory
    if (is.null(base_dir)) {
      base_dir <- dirname(abs_stan_file)
    }
  }

  # Initialize list to store all found files
  all_files <- character(0)
  processed_files <- character(0)

  # Recursive function to process a single file
  process_file <- function(file_path) {
    # Convert to absolute path
    abs_path <- normalizePath(file_path, mustWork = TRUE)

    # Skip if already processed (prevents infinite loops)
    if (abs_path %in% processed_files) {
      return()
    }

    # Mark as processed
    processed_files <<- c(processed_files, abs_path)

    # Read the file
    if (!file.exists(abs_path)) {
      warning(paste("File not found:", abs_path))
      return()
    }

    lines <- readLines(abs_path, warn = FALSE)

    # Find #include statements - improved pattern that properly handles quotes/brackets
    # Pattern captures the path between quotes or angle brackets
    include_pattern <- '^\\s*#include\\s+[\"<]([^\"<>]+)[\">]'
    include_matches <- grep(include_pattern, lines, value = TRUE)

    if (length(include_matches) > 0) {
      # Extract file paths from include statements
      # This will capture the path between the opening and closing quote/bracket
      included_files <- sub(
        '^\\s*#include\\s+[\"<]([^\"<>]+)[\">].*$',
        '\\1',
        include_matches
      )
      # Trim any whitespace from extracted paths
      included_files <- trimws(included_files)

      for (inc_file in included_files) {
        # Handle relative paths
        if (!file.path(inc_file) == inc_file || !startsWith(inc_file, "/")) {
          # Relative path - resolve relative to current file's directory
          inc_path <- file.path(dirname(abs_path), inc_file)
        } else {
          # Absolute path
          inc_path <- inc_file
        }

        # Try to normalize the path
        tryCatch(
          {
            inc_path <- normalizePath(inc_path, mustWork = TRUE)

            # Add to results if not already there
            if (!inc_path %in% all_files) {
              all_files <<- c(all_files, inc_path)
            }

            # Recursively process the included file
            process_file(inc_path)
          },
          error = function(e) {
            # If file doesn't exist, try relative to base_dir (stan root)
            alt_path <- file.path(base_dir, inc_file)
            tryCatch(
              {
                alt_path <- normalizePath(alt_path, mustWork = TRUE)

                if (!alt_path %in% all_files) {
                  all_files <<- c(all_files, alt_path)
                }

                process_file(alt_path)
              },
              error = function(e2) {
                warning(paste(
                  "Could not find included file:",
                  inc_file,
                  "from",
                  abs_path
                ))
              }
            )
          }
        )
      }
    }
  }

  # Start processing from the main file
  process_file(stan_file)

  # Return sorted list of unique file paths
  return(sort(unique(all_files)))
}

# Example usage:
# included_files <- find_stan_includes("model.stan")
# print(included_files)

# nolint end: object_usage_linter

#' Find the first occurrence of n_succ consecutive elements from 'what' in 'all'
#'
#' @param all Integer vector to search within
#' @param what Integer vector of values to search for (can contain duplicates)
#' @param n_succ Number of consecutive elements to find (default 1)
#' @return Starting index (1-based) of the first sequence of n_succ consecutive elements where each element is in 'what'. Returns 0 if not found.
#' @examples
#' find_first(c(1,2,3,2,2,4), c(2,3), 3) # returns 2
find_first <- function(all, what, n_succ = 1) {
  n <- length(all)
  n_what <- length(what)
  if (n < n_succ) {
    return(0)
  }
  if (n_what == 0) {
    stop("find_first: 'what' vector cannot be empty")
  }
  sorted_what <- sort(what)
  for (i in seq_len(n - n_succ + 1)) {
    matches <- 0
    for (j in seq_len(n_succ)) {
      search_val <- all[i + j - 1]
      # Use binary search for efficiency
      found_match <- search_val %in% sorted_what
      if (found_match) {
        matches <- matches + 1
      } else {
        break
      }
    }
    if (matches == n_succ) return(i)
  }
  return(0)
}

# nolint end: object_usage_linter
