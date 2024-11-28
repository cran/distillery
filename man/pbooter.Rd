\name{pbooter}
\alias{pbooter}
\title{
Parametric Bootstrap Resampling
}
\description{
Creates sample statistics for several replicated samples derived by sampling from a parametric distribution.
}
\usage{
pbooter(x, statistic, B, rmodel, rsize, v.terms, verbose = FALSE, ...)
}
\arguments{
  \item{x}{
Original data set.  If it is a vector, then it is assumed to be univariate.  If it is a matrix, it is assumed to be multivariate where each column is a variate.
}
  \item{statistic}{
Function that minimally takes arguments: \code{data} and \code{\dots}.  The argument \code{data} must be the input data for which resamples are taken.  Must return a vector of all desired statistics.
}
  \item{B}{
number of bootstrap resamples to make.
}
  \item{rmodel}{
Function that generates the data to be applied to statistic.  Must have arguments \code{size}, giving the size of the data to be returned, and \code{\dots}.
}
  \item{rsize}{
Number giving the resample size for each bootstrap sample.  If missing and \code{x} is a vector, it will be the length of \code{x}, and if it is a matrix, it will be the number of rows of \code{x}.
}
  \item{v.terms}{
If \code{statistic} returns variance estimates for other parameters, then use this argument to specify the indices returned that give the variance estimates.  There must be a component for every other parameter returned, and they must be in the same order as the other parameters (see examples below).  If an estimate does not exist, an NA should be returned for that spot.
}
  \item{verbose}{logical, should progress information be printed to the screen?}
  \item{\dots}{
Optional arguments to \code{statistic} or \code{rmodel}.
}
}
\details{
Similar functionality to \code{boot} from \pkg{boot} when \code{sim} = \dQuote{parametric}.  In this case, the function is a little simpler, and is intended for use with \code{ci.booted}, or just \code{ci}.  It is similar to \code{booter}, but uses parametric sampling instead of resampling from the original data.
}
\value{
A list object of class \dQuote{booted} is returned with components:
  \item{call }{the function call}
  \item{data }{original data series}
  \item{statistic }{statistic argument passed in}
  \item{statistic.args}{all other arguments passed by \dots}
  \item{B}{ Number of bootstrap replicate samples}
  \item{v.terms}{if variance terms are returned by statistic, the argument is repeated in the returned object.}
  \item{rsize}{the size of the bootstrap resamples.}
  \item{rdata}{rsize by B matrix giving the rmodel generated data.}
  \item{v}{B length vector or B column matrix (if statistic returns a vector) giving the estimated parameter variances for each bootstrap replicate.}
  \item{orig.v}{vector giving the parameter variances (i.e. se^2) of statistic when applied to the original data.}
  \item{original.est}{vector giving the estimated parameter values when statistic is applied to the original data.}
  \item{results}{B length vector or B column matrix giving the parameter estimates for each bootstrap resample.}
  \item{type}{character stating whether the resample method is iid or cbb.}
}

\references{
Efron, B. and Tibshirani, R. J. (1998) \emph{An Introduction to the Bootstrap}.  Chapman and Hall, Boca Raton, Florida, 436 pp.
}
\author{
Eric Gilleland
}

\seealso{
\code{\link{booter}}, \code{\link{ci.booted}}
\code{\link{tibber}}
}
\examples{
z <- rnorm( 100 )

zfun <- function( data, ... ) {

    return( c( mean( data ), var( data ), mean( data^2 ), var( data^2 ) ) )

} # end of 'zfun' function.

rfun <- function( size, ... ) rnorm( size, ... )

res <- pbooter( x = z, statistic = zfun, rmodel = rfun, B = 500,
    rsize = 100, v.terms = c(2, 4) )

print( res )

\dontrun{ci( res )}

}
\keyword{ distribution }% use one of  RShowDoc("KEYWORDS")
\keyword{ datagen }% __ONLY ONE__ keyword per line
