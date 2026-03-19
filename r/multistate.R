# nolint start: object_usage_linter

# =============================================================================
# MULTISTATE PATIENT ROUTING
# =============================================================================
#
# Classification-first approach: validate inputs → classify into 6 patterns →
# derive all Stan fields as pure functions of pattern + ms_mode → validate output.
#
# The 6 patterns record what *happened* to the patient (objective).
# ms_mode controls how the model *treats* each patient (subjective).
#
# Pattern → Stan state mapping:
#
#   Pattern          | none | pfs | illness_death | full
#   -----------------|------|-----|---------------|-----
#   admin_censored   |  0   |  0  |       0       |  0
#   true_dropout     |  0   |  0  |       0       |  3
#   progressed_alive |  0   |  1  |       1       |  1
#   progressed_died  |  0   |  1  |       2       |  2
#   died_on_trial    |  0   |  1  |       2       |  2
#   died_off_trial   |  0   |  0  |       0       |  3


# =============================================================================
# INPUT VALIDATION
# =============================================================================

#' Validate inputs before multistate patient classification
#'
#' Checks that all required columns have correct types, no NAs where they must
#' not appear, and that no impossible boolean combinations are present.
#' Called automatically by \code{classify_ms_patients()}.
#'
#' Invariants checked (I1–I8):
#' \describe{
#'   \item{I1}{original_right_censored is logical, no NAs}
#'   \item{I2}{progression_before_death is logical; NAs allowed only for right-censored patients}
#'   \item{I3}{death is logical, no NAs}
#'   \item{I4}{death_week >= 0 for dead patients, no NAs}
#'   \item{I5}{patient_max_t > 0, no NAs}
#'   \item{I6}{potential_followup >= patient_max_t, no NAs}
#'   \item{I7}{Uncensored patients must have progression or death}
#'   \item{I8}{right_censored implies !progression_before_death}
#' }
#'
#' @param patient_data Tibble with columns: original_right_censored,
#'   progression_before_death, death, death_week, patient_max_t,
#'   potential_followup
#' @return invisible(TRUE) on success, or aborts with informative message
validate_ms_inputs <- function(patient_data) {
  if (!is.logical(patient_data$original_right_censored) ||
      anyNA(patient_data$original_right_censored)) {
    cli::cli_abort(
      "I1: {.field original_right_censored} must be logical with no NAs."
    )
  }

  if (!is.logical(patient_data$progression_before_death)) {
    cli::cli_abort(
      "I2: {.field progression_before_death} must be logical."
    )
  }
  # NAs are allowed for right-censored patients (NA means "not observed" = effectively FALSE).
  # Uncensored patients must have a definitive value.
  uncensored_na_pbd <- !patient_data$original_right_censored &
    is.na(patient_data$progression_before_death)
  if (any(uncensored_na_pbd)) {
    n <- sum(uncensored_na_pbd)
    cli::cli_abort(
      "I2: {n} uncensored patient{?s} ha{?s/ve} NA {.field progression_before_death}."
    )
  }

  if (!is.logical(patient_data$death) || anyNA(patient_data$death)) {
    cli::cli_abort("I3: {.field death} must be logical with no NAs.")
  }

  dead_death_weeks <- patient_data$death_week[patient_data$death]
  if (anyNA(dead_death_weeks) || any(dead_death_weeks < 0)) {
    cli::cli_abort(
      "I4: {.field death_week} must be non-negative with no NAs for dead patients."
    )
  }

  if (anyNA(patient_data$patient_max_t) || any(patient_data$patient_max_t <= 0)) {
    cli::cli_abort(
      "I5: {.field patient_max_t} must be positive with no NAs."
    )
  }

  if (anyNA(patient_data$potential_followup) ||
      any(patient_data$potential_followup < patient_data$patient_max_t)) {
    cli::cli_abort(
      "I6: {.field potential_followup} must be >= {.field patient_max_t} with no NAs."
    )
  }

  impossible <- !patient_data$original_right_censored &
    !patient_data$progression_before_death &
    !patient_data$death
  if (any(impossible)) {
    n <- sum(impossible)
    cli::cli_abort(
      "I7: {n} patient{?s} {?is/are} uncensored with no progression and no death \\
      (impossible combo: !rc & !pbd & !death)."
    )
  }

  # I8: right_censored + progression_before_death = TRUE is a data anomaly
  # (PD recorded after censoring date). These patients are still correctly
  # routed as admin_censored/true_dropout since patterns A/B match first.
  # Warn rather than abort so processing continues.
  i8_flag <- patient_data$original_right_censored &
    !is.na(patient_data$progression_before_death) &
    patient_data$progression_before_death
  if (any(i8_flag)) {
    n <- sum(i8_flag)
    cli::cli_warn(
      "I8: {n} patient{?s} {?is/are} right-censored AND have \\
      {.field progression_before_death} == TRUE (data anomaly; will be routed as censored)."
    )
  }

  invisible(TRUE)
}


# =============================================================================
# CLASSIFICATION
# =============================================================================

#' Classify patients into 6 mutually exclusive multistate patterns
#'
#' Each patient is assigned to exactly one of 6 patterns based on their observed
#' clinical outcome. The classification is purely objective — it records what
#' happened to the patient, independent of how the model will treat them.
#'
#' Patterns:
#' \describe{
#'   \item{admin_censored}{Right-censored close to DCO (potential follow-up nearly exhausted)}
#'   \item{true_dropout}{Right-censored well before DCO (lost to follow-up)}
#'   \item{progressed_alive}{Had RECIST progression, still alive at last follow-up}
#'   \item{progressed_died}{Had RECIST progression, then died}
#'   \item{died_on_trial}{Died without progression, last visit within buffer of death}
#'   \item{died_off_trial}{Died without progression, well after last on-trial visit}
#' }
#'
#' @param patient_data Tibble with columns: original_right_censored,
#'   progression_before_death, death, death_week, patient_max_t,
#'   potential_followup
#' @param admin_censor_buffer Integer number of weeks; patients whose remaining
#'   potential follow-up (or gap between last visit and death) is within this
#'   buffer are classified as admin_censored / died_on_trial respectively.
#'   Default 6L.
#' @return Factor with 6 levels (pattern labels), one entry per patient
classify_ms_patients <- function(patient_data, admin_censor_buffer = 6L) {
  validate_ms_inputs(patient_data)

  factor(
    dplyr::case_when(
      # A: Admin-censored — right-censored, follow-up consumed by DCO
      patient_data$original_right_censored &
        (patient_data$potential_followup - patient_data$patient_max_t <= admin_censor_buffer)
        ~ "admin_censored",

      # B: True dropout — right-censored well before DCO
      patient_data$original_right_censored
        ~ "true_dropout",

      # C: Progressed, still alive at last follow-up
      patient_data$progression_before_death & !patient_data$death
        ~ "progressed_alive",

      # D: Progressed, then died
      patient_data$progression_before_death & patient_data$death
        ~ "progressed_died",

      # E: Died on-trial (death within buffer weeks of last visit)
      patient_data$death & !patient_data$progression_before_death &
        (patient_data$death_week - patient_data$patient_max_t <= admin_censor_buffer)
        ~ "died_on_trial",

      # F: Died off-trial (died well after last on-trial visit)
      patient_data$death & !patient_data$progression_before_death
        ~ "died_off_trial"
    ),
    levels = c(
      "admin_censored", "true_dropout",
      "progressed_alive", "progressed_died",
      "died_on_trial", "died_off_trial"
    )
  )
}


# =============================================================================
# ROUTING VALIDATION
# =============================================================================

#' Validate routing output after multistate field derivation
#'
#' Checks that all derived Stan fields are internally consistent for the given
#' ms_mode. Called automatically by \code{derive_ms_fields()}.
#'
#' Invariants checked (R1–R9):
#' \describe{
#'   \item{R1}{Every patient has exactly one non-NA ms_pattern}
#'   \item{R2}{ms_final_state in \{0, ..., ms_max_state\} for the ms_mode}
#'   \item{R3}{State 0 → ms_censored_01 == 1}
#'   \item{R4}{State 1 → ms_time_01 > 0}
#'   \item{R5}{State 2 → ms_time_01 >= 0 (0 = 0→2 direct, >0 = 0→1→2)}
#'   \item{R6}{State 3 → ms_time_03 == patient_max_t}
#'   \item{R7}{ms_time_12 > 0 only for progressed-then-died patients}
#'   \item{R8}{No state 3 when ms_mode in \{"none", "pfs", "illness_death"\}}
#'   \item{R9}{ms_time_12 >= 1 for state 2 patients on 0→1→2 path}
#' }
#'
#' @param analysis_data Tibble with ms_pattern and patient_max_t columns
#' @param stan_fields Named list output of \code{derive_ms_fields()}
#' @param ms_mode Character, one of "none", "pfs", "illness_death", "full"
#' @return invisible(TRUE) on success, or aborts with informative message
validate_ms_routing <- function(analysis_data, stan_fields, ms_mode) {
  ms_mode <- match.arg(ms_mode, c("none", "pfs", "illness_death", "full"))

  ms_max_state <- switch(ms_mode,
    "none"          = 0L,
    "pfs"           = 1L,
    "illness_death" = 2L,
    "full"          = 3L
  )

  state <- stan_fields$ms_final_state

  if (anyNA(analysis_data$ms_pattern)) {
    n <- sum(is.na(analysis_data$ms_pattern))
    cli::cli_abort("R1: {n} patient{?s} ha{?s/ve} NA ms_pattern.")
  }

  if (any(state < 0L | state > ms_max_state)) {
    bad <- sum(state < 0L | state > ms_max_state)
    cli::cli_abort(
      "R2: {bad} patient{?s} ha{?s/ve} ms_final_state outside \\
      [0, {ms_max_state}] for ms_mode '{ms_mode}'."
    )
  }

  state0 <- state == 0L
  if (any(stan_fields$ms_censored_01[state0] != 1L)) {
    n <- sum(stan_fields$ms_censored_01[state0] != 1L)
    cli::cli_abort("R3: {n} state-0 patient{?s} ha{?s/ve} ms_censored_01 != 1.")
  }

  state1 <- state == 1L
  if (any(state1) && any(stan_fields$ms_time_01[state1] <= 0L)) {
    n <- sum(stan_fields$ms_time_01[state1] <= 0L)
    cli::cli_abort("R4: {n} state-1 patient{?s} ha{?s/ve} ms_time_01 <= 0.")
  }

  state2 <- state == 2L
  if (any(state2) && any(stan_fields$ms_time_01[state2] < 0L)) {
    cli::cli_abort("R5: State-2 patients must have ms_time_01 >= 0.")
  }

  state3 <- state == 3L
  if (any(state3)) {
    mismatch <- stan_fields$ms_time_03[state3] != analysis_data$patient_max_t[state3]
    if (any(mismatch)) {
      n <- sum(mismatch)
      cli::cli_abort("R6: {n} state-3 patient{?s}: ms_time_03 != patient_max_t.")
    }
  }

  # R7: ms_time_12 > 0 only for progressed patients (both C progressed_alive and D progressed_died
  # have valid non-zero sojourn times — D's is an event, C's is censored)
  progressed_pats <- analysis_data$ms_pattern %in% c("progressed_alive", "progressed_died")
  bad_time_12 <- stan_fields$ms_time_12 > 0L & !progressed_pats
  if (any(bad_time_12)) {
    n <- sum(bad_time_12)
    cli::cli_abort(
      "R7: {n} patient{?s} ha{?s/ve} ms_time_12 > 0 but have not progressed."
    )
  }

  # R8: No state 3 when ms_mode in {"none", "pfs", "illness_death"}
  # (already caught by R2 since max_state < 3, but kept as an explicit semantic check)
  if (ms_mode %in% c("none", "pfs", "illness_death") && any(state3)) {
    n <- sum(state3)
    cli::cli_abort(
      "R8: {n} patient{?s} ha{?s/ve} state 3 in ms_mode '{ms_mode}' (not allowed)."
    )
  }

  prog_then_died <- state2 & stan_fields$ms_time_01 > 0L
  if (any(prog_then_died) && any(stan_fields$ms_time_12[prog_then_died] < 1L)) {
    n <- sum(stan_fields$ms_time_12[prog_then_died] < 1L)
    cli::cli_abort(
      "R9: {n} progressed-then-died patient{?s} ha{?s/ve} ms_time_12 < 1."
    )
  }

  invisible(TRUE)
}


# =============================================================================
# FIELD DERIVATION
# =============================================================================

#' Derive all Stan multistate data fields from patient patterns and ms_mode
#'
#' Maps the 6 clinical patterns (from \code{classify_ms_patients()}) to Stan
#' likelihood fields according to the ms_mode. This is the single source of
#' truth for multistate Stan data — replaces both \code{update_dropout_state()}
#' and \code{apply_pfs_ms_mode_if()}.
#'
#' Key time-field semantics:
#' \describe{
#'   \item{ms_time_01}{Time to 0→1 transition or censoring. Uses
#'     \code{analysis_data$pfs} which is already detection-adjusted
#'     (pfs_original + interval_censored + 1 for events). Sentinel 0 = direct
#'     death (0→2 path).}
#'   \item{ms_time_12}{Post-progression sojourn. Computed as
#'     \code{pmax(1, death_week - pfs_adjusted)} — subtracting adjusted pfs
#'     accounts for the detection window shifting progression time forward.}
#'   \item{ms_time_02}{Direct death time. Uses raw \code{death_week} (no
#'     detection adjustment — deaths are observed exactly).}
#' }
#'
#' @param analysis_data Tibble from \code{prepare_analysis_data()} containing
#'   \code{ms_pattern} (from \code{classify_ms_patients()}), \code{pfs}
#'   (detection-adjusted), \code{patient_max_t}, \code{death_week},
#'   \code{ms_prog_deterministic}, and \code{interval_censored}.
#' @param ms_mode Character, one of "none", "pfs", "illness_death", "full"
#' @return Named list of all Stan multistate data fields, ready to merge into
#'   a Stan data list via \code{c()}.
derive_ms_fields <- function(
  analysis_data,
  ms_mode = c("none", "pfs", "illness_death", "full")
) {
  ms_mode <- match.arg(ms_mode)
  pat <- analysis_data$ms_pattern

  # -------------------------------------------------------------------------
  # 0→1 transition
  # -------------------------------------------------------------------------
  # In "pfs" mode: C, D, E are events (composite PFS = progression OR on-trial death)
  # In "illness_death"/"full" mode: C, D only (E goes via 0→2 direct death path)
  is_01_event <- dplyr::case_when(
    ms_mode == "none" ~ FALSE,
    ms_mode == "pfs"  ~ pat %in% c("progressed_alive", "progressed_died", "died_on_trial"),
    TRUE              ~ pat %in% c("progressed_alive", "progressed_died")
  )

  # analysis_data$pfs is already detection-adjusted for event patients.
  # E in illness_death/full mode gets sentinel 0 (signals 0→2 path to Stan).
  ms_time_01 <- dplyr::case_when(
    is_01_event ~ as.integer(analysis_data$pfs),
    pat == "died_on_trial" & ms_mode %in% c("illness_death", "full") ~ 0L,
    TRUE ~ as.integer(analysis_data$patient_max_t)
  )
  ms_censored_01 <- as.integer(!is_01_event)

  # -------------------------------------------------------------------------
  # 0→2 transition (direct death without prior progression)
  # -------------------------------------------------------------------------
  # Only E (died_on_trial) contributes a 0→2 event, and only in illness_death/full.
  # Deaths are observed exactly — no detection-window adjustment.
  is_02_event <- ms_mode %in% c("illness_death", "full") & pat == "died_on_trial"
  ms_time_02 <- dplyr::if_else(
    is_02_event,
    as.integer(analysis_data$death_week),
    as.integer(analysis_data$patient_max_t)
  )
  ms_censored_02 <- as.integer(!is_02_event)

  # -------------------------------------------------------------------------
  # 1→2 transition (post-progression sojourn)
  # -------------------------------------------------------------------------
  # pfs is already detection-adjusted, so death_week - pfs = sojourn after
  # the detection window. pmax(1L, ...) prevents negative sojourn when
  # death occurred in the same week as (or before) the adjusted progression time.
  ms_time_12 <- dplyr::case_when(
    pat == "progressed_died" ~
      pmax(1L, as.integer(analysis_data$death_week - analysis_data$pfs)),
    pat == "progressed_alive" ~
      pmax(1L, as.integer(analysis_data$patient_max_t - analysis_data$pfs)),
    TRUE ~ 0L
  )
  ms_censored_12 <- as.integer(pat != "progressed_died")

  # -------------------------------------------------------------------------
  # 0→3 and 3→2 transitions (dropout / off-trial)
  # -------------------------------------------------------------------------
  ms_time_03 <- as.integer(analysis_data$patient_max_t)
  ms_time_32 <- dplyr::if_else(
    pat == "died_off_trial",
    as.integer(analysis_data$death_week - analysis_data$patient_max_t),
    0L
  )
  ms_censored_32 <- as.integer(pat != "died_off_trial")

  # -------------------------------------------------------------------------
  # Final state (pattern × mode lookup)
  # -------------------------------------------------------------------------
  ms_final_state <- dplyr::case_when(
    ms_mode == "none"                                                              ~ 0L,
    ms_mode == "pfs"                                                               ~ as.integer(is_01_event),
    ms_mode %in% c("illness_death", "full") & pat == "progressed_alive"           ~ 1L,
    ms_mode %in% c("illness_death", "full") & pat %in% c("progressed_died", "died_on_trial") ~ 2L,
    ms_mode == "full" & pat %in% c("true_dropout", "died_off_trial")              ~ 3L,
    TRUE                                                                           ~ 0L
  )

  # -------------------------------------------------------------------------
  # GP grid sizing: max sojourn times + buffer for extrapolation
  # -------------------------------------------------------------------------
  ms_max_sojourn_t <- max(1L, max(ms_time_12, na.rm = TRUE) + 10L)
  ms_max_sojourn_t_32 <- max(1L, max(ms_time_32, na.rm = TRUE) + 10L)

  stan_fields <- list(
    ms_final_state        = ms_final_state,
    ms_time_01            = ms_time_01,
    ms_time_02            = ms_time_02,
    ms_time_12            = ms_time_12,
    ms_censored_01        = ms_censored_01,
    ms_censored_02        = ms_censored_02,
    ms_censored_12        = ms_censored_12,
    ms_time_03            = ms_time_03,
    ms_time_32            = ms_time_32,
    ms_censored_32        = ms_censored_32,
    ms_prog_deterministic = analysis_data$ms_prog_deterministic,
    ms_max_sojourn_t      = ms_max_sojourn_t,
    ms_max_sojourn_t_32   = ms_max_sojourn_t_32,
    ms_gp_grid_step       = 4L
  )

  validate_ms_routing(analysis_data, stan_fields, ms_mode)

  stan_fields
}

# nolint end: object_usage_linter
