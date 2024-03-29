% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dplyr_verbs.R
\name{arrange.flowFrame}
\alias{arrange.flowFrame}
\title{Order rows using column values}
\usage{
\method{arrange}{flowFrame}(.data, ..., .by_group = FALSE)
}
\arguments{
\item{.data}{A \code{\link[flowCore]{flowFrame}}}

\item{...}{Variables, or functions of variables, to arrange by.}

\item{.by_group}{Unused.}
}
\value{
An object of the same type as .data. The output has the following
properties:
* All rows appear in the output, but (usually) in a different place.
* Columns are not modified.
* The \code{\link[flowCore]{flowFrame}}'s
\code{\link[flowCore]{identifier}} will be preserved.
}
\description{
Order rows using column values
}
\examples{

my_flowframe <- simulate_cytometry_data()$flowframe

 my_flowframe |>
   dplyr::arrange(feature_1)

}
