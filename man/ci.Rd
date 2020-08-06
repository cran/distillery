\name{ci}
\alias{ci}
\alias{ci.matrix}
\alias{ci.numeric}
\alias{print.ci}
\title{
Find Confidence Intervals
}
\description{
Method function for finding confidence intervals.
}
\usage{
ci(x, alpha = 0.05, ...)

\method{ci}{matrix}(x, alpha = 0.05, ...)

\method{ci}{numeric}(x, alpha = 0.05, ...)

\method{print}{ci}(x, ...)
}
\arguments{
  \item{x}{
\code{ci}: an R object that has a \code{ci} method function for it.

\code{print}: output from \code{ci}.
}
  \item{alpha}{
number between zero and one giving the 1 - \code{alpha} confidence level.
}
  \item{\dots}{
Optional arguments depending on the specific method function.  In the case of those for \code{ci.matrix} and \code{ci.numeric}, these are any optional arguments to \code{mean} and \code{var}.

Not used by \code{print} method function.
}
}
\details{
\code{ci.numeric}: Calculates the mean and normal approximation CIs for the mean.

\code{ci.matrix}: Does the same as \code{ci.numeric}, but applies to each column of \code{x}.
}
\value{
\code{ci.numeric}: a numeric vector giving the CI bounds and mean value.

\code{ci.matrix}: a matrix giving the mean and CI bounds for each column of \code{x}.
}
\author{
Eric Gilleland
}

\examples{
ci(rnorm(100, mean=10, sd=2))

ci(matrix(rnorm(10000, mean=40, sd=10), 100, 100))
}
\keyword{ htest }
\keyword{ methods}
