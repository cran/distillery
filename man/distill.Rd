\name{distill}
\alias{distill}
\alias{distill.list}
\alias{distill.matrix}
\alias{distill.data.frame}
\title{
Distill An Object
}
\description{
Distill a complex object to something easier to manage, like a numeric vector.
}
\usage{
distill(x, ...)

\method{distill}{list}(x, ...)

\method{distill}{matrix}(x, ...)

\method{distill}{data.frame}(x, ...)

}
\arguments{
  \item{x}{
A list, vector, matrix or data frame, or other object that has a \code{distill} method, e.g., \code{fevd} objects.
}
  \item{\dots}{
Not used.
}
}
\details{
Perhaps a fine line exists between functions such as \code{c}, \code{print}, \code{str}, \code{summary}, etc.  The idea behind the \code{distill} method is to have a function that \dQuote{distills} out the most pertinent information from a more complex object.  For example, when fitting a model to a number of spatial locations, it can be useful to pull out only certain information into a vector for ease of analysis.  With many models, it might not be feasible to store (or analyze) large complicated data objects.  In such a case, it may be useful to keep only a vector with the most pertinent information (e.g., parameter estimates, their standard errors, the likelihood value, AIC, BIC, etc.).  For example, this is used within \pkg{extRemes} >= 2.0 on the \dQuote{fevd} class objects with the aim at fitting models to numerous locations within an \code{apply} call so that something easily handled is returned, but with enough information as to be useful.

The data frame and matrix methods attempt to name each component of the vector.  The list method simply does \code{c(unlist(x))}.
}
\value{
numeric vector, possibly named.
}
\author{
Eric Gilleland
}
\seealso{
\code{\link{c}}, \code{\link{unlist}}, \code{\link{print}}, \code{\link{summary}}, \code{\link{str}}, \code{\link{args}}
}
\examples{
x <- cbind(1:3, 4:6, 7:9)
distill(x)

x <- data.frame(x=1:3, y=4:6, z=7:9)
distill(x)
}
\keyword{ manip }
\keyword{ methods }
\keyword{ misc }
