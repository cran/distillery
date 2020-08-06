\name{distillery-package}
\alias{distillery-package}
\alias{distillery}
\docType{package}
\title{
distillery: Methods to Distill Information from R Objects
}
\description{
\pkg{distillery} contains primarily method functions to distill out pertinent information from R objects, as well as to compute confidence intervals.  It now also contains new fairly general bootstrap functions.
}

\details{

Primary functions include:

\code{distill}: Typically, to distill pertinent information from a complicated (usually a list) object and return a named vector.

\code{ci}: Calculate confidence intervals.  This is a method function for calculating confidence intervals.  Includes methods for numeric vectors and matrices, whereby the mean is taken (column-wise for matrices) and normal approximation confidence intervals for the mean are calculated and returned.

\code{booter}, \code{pbooter} and \code{tibber}: Functions to perform bootstrap resampling that work with \code{ci} (\code{booter} and \code{pbooter}).  Allows for m < n bootstrap resampling, circular block bootstrapping, parametric bootstrap resampling (\code{pbooter}), and the test-inversion bootstrap approach (\code{tibber}).

}

\author{
Eric Gilleland
}
\examples{
## See help files for above named functions and datasets
## for specific examples.
}
\keyword{package}
