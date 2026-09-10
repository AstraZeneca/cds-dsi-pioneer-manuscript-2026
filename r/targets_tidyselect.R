#' Tidyselect helpers for targets dependency graphs
#'
#' Custom tidyselect functions that allow selecting targets based on
#' their position in the dependency graph.

#' Get downstream dependents from a dependency edge list
#'
#' @param edges A data frame with 'from' and 'to' columns
#' @param start_names Character vector of starting target names
#' @return Character vector of all downstream target names (excluding start_names)
get_downstream <- function(edges, start_names) {
  if (nrow(edges) == 0 || length(start_names) == 0) {
    return(character())
  }

  graph <- igraph::graph_from_data_frame(edges, directed = TRUE)

  start_names |>
    intersect(igraph::V(graph)$name) |>
    purrr::map(~ igraph::subcomponent(graph, .x, mode = "out")$name) |>
    purrr::reduce(union, .init = character()) |>
    setdiff(start_names)
}

#' Get upstream dependencies from a dependency edge list
#'
#' @param edges A data frame with 'from' and 'to' columns
#' @param start_names Character vector of starting target names
#' @return Character vector of all upstream target names (including start_names)
get_upstream <- function(edges, start_names) {
  if (nrow(edges) == 0 || length(start_names) == 0) {
    return(start_names)
  }

  graph <- igraph::graph_from_data_frame(edges, directed = TRUE)

  start_names |>
    intersect(igraph::V(graph)$name) |>
    purrr::map(~ igraph::subcomponent(graph, .x, mode = "in")$name) |>
    purrr::reduce(union, .init = character()) |>
    union(start_names)
}

#' Select targets and all their downstream dependents (subtrees)
#'
#' @param expr A tidyselect expression (e.g., starts_with("model_"), contains("fit"))
#' @return Integer vector of column positions for use in tidyselect context
#' @export
#'
#' @examples
#' \dontrun{
#' # Single target and its dependents
#' tar_make(names = depends_on(my_data))
#'
#' # Pattern matching
#' tar_make(names = depends_on(starts_with("model_")))
#'
#' # Multiple patterns
#' tar_make(names = depends_on(contains("fit") | contains("posterior")))
#' }
depends_on <- function(expr) {
  expr_quo <- rlang::enquo(expr)
  vars <- tidyselect::peek_vars(fn = "depends_on")

  selected_idx <- tidyselect::eval_select(
    expr_quo,
    data = rlang::set_names(seq_along(vars), vars)
  )
  selected_names <- vars[selected_idx]

  if (length(selected_names) == 0) {
    return(integer(0))
  }

  edges <- targets::tar_network(targets_only = TRUE)$edges

  all_targets <- get_downstream(edges, selected_names)

  which(vars %in% all_targets)
}

#' Select targets and all their upstream dependencies
#'
#' @param expr A tidyselect expression (e.g., ends_with("_summary"), matches("^report"))
#' @return Integer vector of column positions for use in tidyselect context
#' @export
#'
#' @examples
#' \dontrun{
#' # Single target and its dependencies
#' tar_make(names = upstream_of(final_report))
#'
#' # All dependencies of summary targets
#' tar_outdated(names = upstream_of(ends_with("_summary")))
#' }
upstream_of <- function(expr) {
  expr_quo <- rlang::enquo(expr)
  vars <- tidyselect::peek_vars(fn = "upstream_of")

  selected_idx <- tidyselect::eval_select(
    expr_quo,
    data = rlang::set_names(seq_along(vars), vars)
  )
  selected_names <- vars[selected_idx]

  if (length(selected_names) == 0) {
    return(integer(0))
  }

  edges <- targets::tar_network(targets_only = TRUE)$edges

  all_targets <- get_upstream(edges, selected_names)

  which(vars %in% all_targets)
}
