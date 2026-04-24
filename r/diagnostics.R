# r/diagnostics.R
# Post-sampling diagnostics for posterior fits.
# Covers two orthogonal concerns:
#   1. Convergence: R-hat and ESS (check_convergence)
#   2. Sampling quality: NUTS chain diagnostics (summarize_nuts)
#
# Both functions are designed to consume already-extracted targets
# (e.g. *_draws_pop, *_draws_patient_params, *_nuts_param) so no
# additional CSV reads are required at diagnostic time.

# check_convergence -----------------------------------------------------------

#' Tiered R-hat and ESS check for hierarchical model draws
#'
#' Structural parameters (population-level, SDs, covariate coefficients) are
#' summarised individually. Patient-level effects — too numerous to inspect
#' one-by-one — are reduced to aggregate worst-case statistics.
#'
#' @param draws_pop     draws_array of structural parameters
#' @param draws_patient draws_array of patient-level parameters
#' @return Named list:
#'   $pop     — tibble with variable, rhat, ess_bulk, ess_tail, sorted by rhat desc
#'   $patient — single-row tibble: group, n_params, max_rhat, n_rhat_bad,
#'              min_ess_bulk, min_ess_tail
check_convergence <- function(draws_pop, draws_patient) {
  pop_summary <- draws_pop |>
    posterior::summarize_draws(
      rhat     = posterior::rhat,
      ess_bulk = posterior::ess_bulk,
      ess_tail = posterior::ess_tail
    ) |>
    dplyr::arrange(dplyr::desc(.data$rhat))

  patient_per_var <- draws_patient |>
    posterior::summarize_draws(
      rhat     = posterior::rhat,
      ess_bulk = posterior::ess_bulk,
      ess_tail = posterior::ess_tail
    )

  patient_summary <- tibble::tibble(
    group        = "patient_effects",
    n_params     = nrow(patient_per_var),
    max_rhat     = max(patient_per_var$rhat, na.rm = TRUE),
    n_rhat_bad   = as.integer(sum(patient_per_var$rhat > 1.01, na.rm = TRUE)),
    min_ess_bulk = min(patient_per_var$ess_bulk, na.rm = TRUE),
    min_ess_tail = min(patient_per_var$ess_tail, na.rm = TRUE)
  )

  list(pop = pop_summary, patient = patient_summary)
}

# summarize_nuts --------------------------------------------------------------

#' Per-chain NUTS sampler quality summary (sampling phase only)
#'
#' Expects the output of `bayesplot::nuts_params(fit)`, which by default
#' contains only post-warmup (sampling) iterations.
#'
#' E-BFMI (Energy Bayesian Fraction of Missing Information) is computed as
#' var(diff(energy__)) / var(energy__) per chain. Values below 0.3 suggest
#' the kinetic energy is not exploring the posterior geometry well.
#'
#' @param nuts_params_df Data frame from `bayesplot::nuts_params(fit)`
#'   with columns: Chain, Iteration, Parameter, Value
#' @param max_treedepth Maximum tree depth threshold (default 10, matching
#'   Stan's default `max_depth` setting). Hits are counted as
#'   treedepth__ >= max_treedepth.
#' @return Tibble with one row per chain:
#'   Chain, n_iter, n_divergent, n_max_treedepth, mean_accept,
#'   median_stepsize, total_leapfrog, ebfmi
summarize_nuts <- function(nuts_params_df, max_treedepth = 10) {
  nuts_params_df |>
    tidyr::pivot_wider(names_from = Parameter, values_from = Value) |>
    dplyr::summarise(
      .by = Chain,
      n_iter           = dplyr::n(),
      n_divergent      = as.integer(sum(divergent__)),
      n_max_treedepth  = as.integer(sum(treedepth__ >= max_treedepth)),
      mean_accept      = mean(accept_stat__),
      median_stepsize  = stats::median(stepsize__),
      total_leapfrog   = sum(n_leapfrog__),
      ebfmi            = stats::var(diff(energy__)) / stats::var(energy__)
    ) |>
    dplyr::arrange(Chain)
}

# energy_correlations --------------------------------------------------------

#' Rank parameters by correlation with HMC energy
#'
#' Low E-BFMI indicates the Hamiltonian kinetic energy is not exploring the
#' posterior geometry well. Parameters whose values correlate strongly with
#' the per-iteration energy are the dimensions driving that mis-match, and
#' are therefore the candidate culprits for geometry-induced pathologies
#' (heavy tails, funnels).
#'
#' @param draws         posterior::draws_array (or anything as_draws_df accepts)
#' @param nuts_param_df bayesplot::nuts_params output aligned to the same draws
#' @return Tibble ordered by |cor(energy, theta)| descending:
#'   variable, cor_energy, abs_cor, ess_bulk, rhat
energy_correlations <- function(draws, nuts_param_df) {
  np_wide <- nuts_param_df |>
    tidyr::pivot_wider(names_from = Parameter, values_from = Value) |>
    dplyr::arrange(Chain, Iteration)

  draws_df <- posterior::as_draws_df(draws) |>
    dplyr::arrange(.chain, .iteration)

  stopifnot(nrow(draws_df) == nrow(np_wide))

  theta_mat <- draws_df |>
    dplyr::select(-.chain, -.iteration, -.draw) |>
    as.matrix()

  # Some columns (derived fully deterministic quantities) may be constant;
  # cor() emits NaN for those — just drop them from the ranking.
  cors <- suppressWarnings(stats::cor(np_wide$energy__, theta_mat)[1, ])

  summary_df <- posterior::summarize_draws(
    draws,
    ess_bulk = posterior::ess_bulk,
    rhat     = posterior::rhat
  )

  tibble::tibble(
    variable   = colnames(theta_mat),
    cor_energy = unname(cors),
    abs_cor    = abs(unname(cors))
  ) |>
    dplyr::filter(!is.na(.data$abs_cor)) |>
    dplyr::left_join(summary_df, by = "variable") |>
    dplyr::arrange(dplyr::desc(.data$abs_cor))
}

# group_ess_summary ----------------------------------------------------------

#' Aggregate ESS / R-hat by parameter group prefix
#'
#' Expands the scalar `check_convergence()$patient` summary: instead of one
#' row for "all patient effects", produce one row per group prefix (e.g.
#' `frac_log_growth_patient`, `init_raw_level_intercept`). Useful for
#' pinpointing which hierarchical layer is under-mixing.
#'
#' @param draws posterior::draws_array
#' @return Tibble ordered by min_ess_bulk ascending:
#'   group, n_params, min_ess_bulk, median_ess_bulk, max_rhat, n_rhat_bad
group_ess_summary <- function(draws) {
  posterior::summarize_draws(
    draws,
    ess_bulk = posterior::ess_bulk,
    rhat     = posterior::rhat
  ) |>
    dplyr::mutate(
      group = sub("\\[.*\\]$", "", .data$variable)
    ) |>
    dplyr::summarise(
      .by             = group,
      n_params        = dplyr::n(),
      min_ess_bulk    = min(.data$ess_bulk, na.rm = TRUE),
      median_ess_bulk = stats::median(.data$ess_bulk, na.rm = TRUE),
      max_rhat        = max(.data$rhat, na.rm = TRUE),
      n_rhat_bad      = as.integer(sum(.data$rhat > 1.01, na.rm = TRUE))
    ) |>
    dplyr::arrange(.data$min_ess_bulk)
}
