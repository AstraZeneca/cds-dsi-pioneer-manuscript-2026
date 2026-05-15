# Multi-Level Hierarchy Helper Functions
#
# Functions to create data structures for Stan's multi-level hierarchy system.
# Supports arbitrary N-level hierarchies (e.g., trial -> region -> site -> patient).

# level_intercept_mode: controls whether a hierarchy level has no effect,
# a fixed effect (FE, SD fixed as hyperparameter), or a random effect (RE,
# SD estimated from data via NCP). Consistent across tr, frac, init, and ms
# modules. See stan/modules/tr/flags.stan and stan/modules/multistate/flags.stan.
# Integer values must match the LEVEL_MODE_* constants in
# stan/_hierarchy_transformed_data.stan.
level_intercept_mode <- c(
  none = 0L,   # No intercept at this level
  fe = 1L,     # Fixed effect: SD is a data hyperparameter (no pooling)
  re = 2L,     # Random effect: SD estimated via NCP (hierarchical pooling)
  re_gp = 3L,  # Random effect + full GP residual (ms module only, for now)
  re_cp = 4L   # Random effect, centered: direct ~normal(0, sd) sampling
)

#' Create multi-level hierarchy structure for Stan
#'
#' @param patient_assignments Named list of vectors, one per level (excluding patient level).
#'   Names are level names (e.g., "trial", "region"), values are integer vectors
#'   mapping each patient to their group at that level.
#' @param n_patients Number of patients (used to create patient level)
#' @return List with hierarchy structure for Stan:
#'   - n_levels: Total number of levels including patients
#'   - n_groups_per_level: Array of group counts per level
#'   - patient_level_groups: Matrix mapping patients to groups at each level
#'   - level_names: Character vector of level names (for documentation)
create_hierarchy_structure <- function(patient_assignments, n_patients) {

  n_levels <- length(patient_assignments) + 1L  # +1 for patient level

  # Preserve level names for documentation
  level_names <- c(names(patient_assignments), "patient")

  # Convert factors to integers (factors are common from data frames)
  patient_assignments <- map(patient_assignments, as.integer)

  # Build patient_level_groups matrix
  patient_level_groups <- do.call(cbind, patient_assignments)
  patient_level_groups <- cbind(patient_level_groups, seq_len(n_patients))
  colnames(patient_level_groups) <- level_names

  # Calculate n_groups_per_level (named for clarity)
  n_groups_per_level <- c(
    map_int(patient_assignments, max),
    n_patients
  )
  names(n_groups_per_level) <- level_names

  # Note: level_names stored as attribute to avoid passing to Stan
  result <- list(
    n_levels = n_levels,
    n_groups_per_level = as.array(n_groups_per_level),
    patient_level_groups = patient_level_groups
  )
  attr(result, "level_names") <- level_names
  result
}

#' Validate an SD sub-hierarchy mode tibble against a level stack
#'
#' Checks each row of `sd_modes` against the current implementation's positional
#' (nested) rule: `sub_level` must appear positionally before `location_level` in
#' the `level_stack`. Empty tibbles pass trivially.
#'
#' @param sd_modes Tibble with columns `location_level`, `sub_level`, `mode`.
#'   An empty tibble (zero rows) is valid and represents the all-NONE configuration.
#' @param level_stack Character vector naming the levels in positional order,
#'   e.g. `c("trial", "arm", "patient")`.
#' @param allowed_modes Character vector of permitted mode strings.
#'   Defaults to `c("none", "fe", "re", "re_cp")`. `"gp"` is NOT in this default
#'   because the sub-hierarchy does not support RE_GP.
#'
#' @return Invisibly returns `sd_modes` if all rows are valid; stops otherwise.
validate_sd_modes <- function(sd_modes,
                              level_stack,
                              allowed_modes = c("none", "fe", "re", "re_cp")) {
  stopifnot(is.data.frame(sd_modes))
  if (nrow(sd_modes) == 0) return(invisible(sd_modes))

  required_cols <- c("location_level", "sub_level", "mode")
  missing <- setdiff(required_cols, names(sd_modes))
  if (length(missing) > 0) {
    stop("validate_sd_modes: missing columns ", paste(missing, collapse = ", "))
  }

  level_pos <- setNames(seq_along(level_stack), level_stack)

  for (i in seq_len(nrow(sd_modes))) {
    row <- sd_modes[i, ]
    loc <- row$location_level
    sub <- row$sub_level
    mode <- row$mode

    if (!(loc %in% level_stack)) {
      stop(sprintf(
        "validate_sd_modes: row %d (%s, %s): unknown location_level %s; not in level_stack = (%s)",
        i, loc, sub, loc, paste(level_stack, collapse = ", ")
      ))
    }
    if (!(sub %in% level_stack)) {
      stop(sprintf(
        "validate_sd_modes: row %d (%s, %s): unknown sub_level %s; not in level_stack = (%s)",
        i, loc, sub, sub, paste(level_stack, collapse = ", ")
      ))
    }
    if (level_pos[sub] >= level_pos[loc]) {
      stop(sprintf(
        "validate_sd_modes: row %d (%s, %s): sub_level must appear positionally before location_level in the (nested) level_stack",
        i, loc, sub
      ))
    }
    if (!(mode %in% allowed_modes)) {
      stop(sprintf(
        "validate_sd_modes: row %d (%s, %s): mode '%s' not in allowed_modes = (%s)",
        i, loc, sub, mode, paste(allowed_modes, collapse = ", ")
      ))
    }
  }

  invisible(sd_modes)
}
