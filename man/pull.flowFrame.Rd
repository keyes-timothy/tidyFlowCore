% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dplyr_verbs.R
\name{pull.flowFrame}
\alias{pull.flowFrame}
\title{Extract a single column.}
\usage{
\method{pull}{flowFrame}(.data, var = -1, name = NULL, ...)
}
\arguments{
\item{.data}{A \code{\link[flowCore]{flowFrame}}.}

\item{var}{A variable specified as:
* a literal variable name
* a positive integer, giving the position counting from the left
* a negative integer, giving the position counting from the right.}

\item{name}{An optional parameter that specifies the column to be used as
names for a named vector. Specified in a similar manner as var.}

\item{...}{For use by methods.}
}
\value{
A vector the same size as .data.
}
\description{
pull() is similar to $. It's mostly useful because it looks a little nicer
in pipes.
}
\examples{

my_flowframe <- simulate_cytometry_data()$flowframe

 my_flowframe |>
   dplyr::pull(feature_1)

}
