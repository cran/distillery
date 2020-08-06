\name{is.formula}
\alias{is.formula}
\title{
Is the R Object a Formula
}
\description{
Tests to see if an object is a formula or not.
}
\usage{
is.formula(x)
}
\arguments{
  \item{x}{
An R object.
}
}
\details{
This function is a very simple one that simplifies checking whether or not the class of an object is a formula or not.
}
\value{
single logical
}
\author{
Eric Gilleland
}

\examples{
is.formula(~1)
is.formula(1:3)
}
\keyword{ classes }
