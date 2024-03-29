% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dplyr_verbs.R
\name{mutate.flowSet}
\alias{mutate.flowSet}
\title{Create, modify, and delete columns.}
\usage{
\method{mutate}{flowSet}(.data, ...)
}
\arguments{
\item{.data}{A \code{\link[flowCore]{flowSet}}}

\item{...}{Name-value pairs. The name (the left side of the equals sign)
gives the name of the column in the output. The right side of the equation
performs computations using the names of each channel according to
 \code{\link[flowCore]{featureNames}}. Supports tidyselection.}
}
\value{
A \code{\link[flowCore]{flowSet}}. The output has the following properties:
* Columns from .data will be preserved according to the .keep argument.
* Existing columns that are modified by ... will always be returned in their original location.
* New columns created through ... will be placed according to the .before and .after arguments.
* The number of rows is not affected.
* Columns given the value NULL will be removed.
}
\description{
Create, modify, and delete columns.
}
\examples{
my_flowset <-
  simulate_cytometry_data()$flowset

my_flowset |>
  dplyr::mutate(new_feature = feature_1 + feature_2)

}
