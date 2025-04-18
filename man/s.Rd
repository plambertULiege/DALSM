% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/s.R
\name{s}
\alias{s}
\title{Internal function to specify smooth terms in DALSM formulas}
\usage{
s(x, ref = NULL)
}
\arguments{
\item{x}{Name of the variable for which an additive term is requested.}

\item{ref}{(Optional) reference value for \code{x} where the additive term is zero.}
}
\value{
The submitted variable for which an additive term is required.
}
\description{
Internal function to specify smooth terms in DALSM formulas
}
\keyword{internal}
