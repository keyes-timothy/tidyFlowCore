% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dplyr_verbs.R
\name{select.flowFrame}
\alias{select.flowFrame}
\title{Keep or drop columns using their names and types.}
\usage{
\method{select}{flowFrame}(.data, ...)
}
\arguments{
\item{.data}{A \code{\link[flowCore]{flowFrame}}}

\item{...}{One or more unquoted expressions separated by commas. Variables names
(as specified by \code{\link[flowCore]{featureNames}}) can be used as if they
were positions in the \code{\link[flowCore]{flowFrame}}). Supports tidyselection.}
}
\value{
A \code{\link[flowCore]{flowFrame}}. The output has the following properties:
* Rows are not affected.
* Output columns are a subset of input columns, potentially with a different order. Columns will be renamed if new_name = old_name form is used.
* The \code{\link[flowCore]{flowFrame}}'s \code{\link[flowCore]{identifier}} will be preserved.
}
\description{
Keep or drop columns using their names and types.
}
\examples{

my_flowframe <- simulate_cytometry_data()$flowframe

 my_flowframe |>
   dplyr::select(feature_1)

}
