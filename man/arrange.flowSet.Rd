% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dplyr_verbs.R
\name{arrange.flowSet}
\alias{arrange.flowSet}
\title{Order rows using column values}
\usage{
\method{arrange}{flowSet}(.data, ..., .by_group = FALSE)
}
\arguments{
\item{.data}{A \code{\link[flowCore]{flowSet}}}

\item{...}{Variables, or functions of variables, to arrange by.}

\item{.by_group}{Unused.}
}
\value{
An object of the same type as .data. The output has the following
properties:
* All rows appear in the output, but (usually) in a different place.
* Columns are not modified.
* The \code{\link[flowCore]{flowSet}}'s
\code{\link[Biobase]{pData}} will be preserved.
}
\description{
Order rows using column values
}
\examples{
my_flowset <- simulate_cytometry_data()$flowset

my_flowset |>
  dplyr::arrange(feature_1)


}
