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
#'   Chain, n_iter, n_divergent, n_max_treedepth, mean_accept, ebfmi
summarize_nuts <- function(nuts_params_df, max_treedepth = 10) {
  nuts_params_df |>
    tidyr::pivot_wider(names_from = Parameter, values_from = Value) |>
    dplyr::summarise(
      .by = Chain,
      n_iter           = dplyr::n(),
      n_divergent      = as.integer(sum(divergent__)),
      n_max_treedepth  = as.integer(sum(treedepth__ >= max_treedepth)),
      mean_accept      = mean(accept_stat__),
      ebfmi            = stats::var(diff(energy__)) / stats::var(energy__)
    ) |>
    dplyr::arrange(Chain)
}
