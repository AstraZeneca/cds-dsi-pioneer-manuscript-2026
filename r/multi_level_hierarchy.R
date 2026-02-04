# Multi-Level Hierarchy Helper Functions
#
# Functions to create data structures for Stan's multi-level hierarchy system.
# Supports arbitrary N-level hierarchies (e.g., trial -> region -> site -> patient).

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
