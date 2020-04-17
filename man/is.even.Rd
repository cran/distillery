\name{is.even}
\alias{is.even}
\alias{is.odd}
\alias{even}
\alias{odd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Identify Even or Odd Numbers
}
\description{
	Simple functions to test for or return the even or odd numbers.
}
\usage{
is.even(x)
is.odd(x)
even(x)
odd(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ any numeric, but maybe makes the most sense with integers.}
}
\details{
	Return a logical vector/matrix of the same dimension as the argument \code{x} telling whether each component is odd (\code{is.odd}) or even (\code{is.even}), or return just the even (\code{even}) or odd (\code{odd}) numbers from the vector/matrix.  Uses \code{\%\%}.
}
\value{
	Returns a logical vector/matrix/array of the same dimension as x in the case of is.even and is.odd, and returns a vector of length less than or equal to x in the case of even and odd; or if no even/odd values, returns integer(0).
}
\author{
	Eric Gilleland
}
\seealso{
\code{\%\%} 
}
\examples{
is.even( 1:7 )
is.odd( 1:7 )
even( 1:7 )
odd( 1:7 )
}
\keyword{ manip }% use one of  RShowDoc("KEYWORDS")
\keyword{ arith }% __ONLY ONE__ keyword per line
\keyword{ logic }
