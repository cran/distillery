\name{datagrabber}
\alias{datagrabber}
\title{
Get Original Data from an R Object
}
\description{
Get the original data set used to obtain the resulting R object for which a method function exists.
}
\usage{
datagrabber(x, ...)

}
\arguments{
  \item{x}{
An R object that has a method function for \code{datagrabber}.
}
  \item{\dots}{
Not used.
}
}
\details{
Often when applying functions to data, it is handy to be able to grab the original data for subsequent routines (e.g., plotting, etc.).  In some cases, information about where to obtain the original data might be available (more difficult) and in other cases, the data may simply be contained within a fitted object.  This method function is generic, but some packages (e.g., \pkg{extRemes} >= 2.0, \pkg{SpatialVx} >= 1.0) have datagrabber functions specific to particular object types.
}
\value{
The original pertinent data in whatever form it takes.
}
\author{
Eric Gilleland
}

\examples{
\dontrun{
## From the extRems (>= 2.0) package.
y <- rnorm(100, mean=40, sd=20)
y <- apply(cbind(y[1:99], y[2:100]), 1, max)
bl <- rep(1:3, each=33)

ydc <- decluster(y, quantile(y, probs=c(0.95)), r=1, blocks=bl)

yorig <- datagrabber(ydc)
all(y - yorig == 0)

}

}
\keyword{ data }
\keyword{ manip }
