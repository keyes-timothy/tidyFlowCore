# tidyr_verbs.R
# This file contains implementations of tidyr verbs for flowCore's flowFrame
# and flowSet S4 classes.


#' Nest a \code{\link[flowCore]{flowFrame}} into a \code{\link[flowCore]{flowSet}}
#'
#' @param .data A \code{\link[flowCore]{flowFrame}}
#'
#' @param ... Columns to nest; these will appear in the inner \code{\link[flowCore]{flowFrame}}s
#' comprising the output \code{\link[flowCore]{flowSet}}.
#' Specified using name-variable pairs of the form new_col = c(col1, col2, col3).
#' The right hand side can be any valid tidyselect expression.
#' If not supplied, then ... is derived as all columns not selected by .by.
#'
#' @param .by Columns to nest by; these will be stored in the
#' \code{\link[Biobase]{pData}} of the output \code{\link[flowCore]{flowSet}}.
#' .by can be used in place of or in conjunction with columns supplied through ....
#' If not supplied, then .by is derived as all columns not selected by ....
#'
#' @param .key Unused.
#'
#' @param .names_sep Unused.
#'
#' @returns A \code{\link[flowCore]{flowSet}} wherein cells are grouped into
#' constituent \code{\link[flowCore]{flowFrame}}s based on which columns are used
#' to nest.
#'
#' @importFrom rlang enquos
#' @importFrom rlang !!
#' @importFrom rlang !!!
#'
#' @export
#'
#' @examples
#'
#' my_flowframe <-
#'   simulate_cytometry_data()$flowframe |>
#'   dplyr::mutate(
#'     random_group =
#'       sample(
#'         c("a", "b"),
#'         size = nrow(simulate_cytometry_data()$flowframe),
#'         replace = TRUE
#'       )
#'   )
#' my_flowframe |>
#'   tidyr::nest(.by = random_group)
#'
nest.flowFrame <- function(.data, ..., .by = NULL, .key = NULL, .names_sep = NULL) {

  .cols <- rlang::enquos(...)
  if (length(.cols) == 0) {
    result <-
      .data |>
      group_by.flowFrame({{ .by }})
  } else if (length(.cols) == 1) {
    names(.cols) <- NULL
    .cols <- .cols[[1]]
    result <-
      .data |>
      group_by.flowFrame(!! .cols)
  } else {
    names(.cols) <- NULL
    result <-
      .data |>
      group_by.flowFrame(!!! .cols)
  }

  return(result)
}

#' Unnest a \code{\link[flowCore]{flowSet}} into a single
#' \code{\link[flowCore]{flowFrame}}
#'
#' @param data A \code{\link[flowCore]{flowSet}}
#'
#' @param cols Columns in \code{\link[Biobase]{pData}} to unnest.
#'
#' @param ... Unused.
#'
#' @param keep_empty Unused.
#'
#' @param ptype Unused.
#'
#' @param names_sep Unused.
#'
#' @param names_repair Unused.
#'
#' @returns A \code{\link[flowCore]{flowFrame}} or \code{\link[flowCore]{flowSet}}
#' depending on the degree of unnest-ing. Note that unnest-ing and
#' ungrouping a \code{\link[flowCore]{flowSet}} are equivalent.
#'
#' @export
#'
#' @examples
#' my_flowset <- simulate_cytometry_data()$flowset
#'
#' my_flowset |>
#'   tidyr::unnest(cols = c(patient, cell_type))
#'
#' my_flowset |>
#'   tidyr::unnest(cols = patient)
#'
unnest.flowSet <-
  function (
    data,
    cols,
    ...,
    keep_empty = FALSE,
    ptype = NULL,
    names_sep = NULL,
    names_repair = "check_unique"
  ) {
    result <-
      data |>
      ungroup.flowSet({{ cols }})

    return(result)
  }
