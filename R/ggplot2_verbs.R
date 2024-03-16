# ggplot2_verbs.R

#' Create a new ggplot.
#'
#' @param data Default dataset to use for plot in the form of a
#' \code{\link[flowCore]{flowFrame}}. If not specified, must be supplied in each
#'  layer added to the plot.
#'
#' @param mapping Default list of aesthetic mappings to use for plot.
#' If not specified, must be supplied in each layer added to the plot. Note that
#' variable names used for aesthetic mappings come from the
#' \code{\link[flowCore]{featureNames}} of the input \code{\link[flowCore]{flowFrame}}.
#'
#' @param ... Other arguments passed on to methods. Not currently used.
#'
#' @param environment Deprecated. Used prior to tidy evaluation.
#'
#' @returns A \code{\link[ggplot2]{ggplot}}
#'
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 ggplot
#'
#' @export
#'
#' @examples
#' simulations <- simulate_cytometry_data()
#' test_flowframe <- simulations$flowframe
#'
#' flowframe_plot <-
#'   test_flowframe |>
#'   ggplot2::ggplot(ggplot2::aes(x = feature_1, y = feature_2)) +
#'   ggplot2::geom_point()
#'
ggplot.flowFrame <-
  function(data = NULL, mapping = ggplot2::aes(), ..., environment = parent.frame()) {
    tof_tibble <-
      as_tof_tbl(data, .name_method = "featureNames")
    result <-
      ggplot2::ggplot(data = tof_tibble, mapping = mapping, ..., environment = environment)
    return(result)
  }

#' Create a new ggplot.
#'
#' @param data Default dataset to use for plot in the form of a
#' \code{\link[flowCore]{flowSet}}. If not specified, must be supplied in each
#'  layer added to the plot. Note that any metadata stored in
#'  \code{\link[Biobase]{pData}} will be merged into the underlying
#'  flowCore-tibble abstraction and will thus be available for plotting.
#'
#' @param mapping Default list of aesthetic mappings to use for plot.
#' If not specified, must be supplied in each layer added to the plot. Note that
#' variable names used for aesthetic mappings come from the
#' \code{\link[flowCore]{featureNames}} of the input
#' \code{\link[flowCore]{flowSet}}'s constituent
#' \code{\link[flowCore]{flowFrame}}s.
#'
#' @param ... Other arguments passed on to methods. Not currently used.
#'
#' @param environment Deprecated. Used prior to tidy evaluation.
#'
#' @returns A \code{\link[ggplot2]{ggplot}}
#'
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 ggplot
#'
#' @export
#'
#' @examples
#' simulations <- simulate_cytometry_data()
#' test_flowset <- simulations$flowset
#'
#' flowset_plot <-
#'   test_flowset |>
#'   ggplot2::ggplot(ggplot2::aes(x = feature_1, y = feature_2)) +
#'   ggplot2::geom_point()
#'
#' flowset_plot_with_metadata <-
#'   test_flowset |>
#'   # note that `patient` below comes from the flowSet's metadata (pData)
#'   ggplot2::ggplot(ggplot2::aes(x = feature_1, y = feature_2, color = patient)) +
#'   ggplot2::geom_point()
ggplot.flowSet <-
  function(data = NULL, mapping = ggplot2::aes(), ..., environment = parent.frame()) {
    tof_tibble <-
      as_tof_tbl(data, .name_method = "featureNames", include_metadata = TRUE)
    result <-
      ggplot2::ggplot(
        data = tof_tibble,
        mapping = mapping,
        ...,
        environment = environment
      )
    return(result)
  }
