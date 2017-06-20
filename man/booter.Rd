\name{booter}
\alias{booter}
\title{
Bootstrap Resampling
}
\description{
Generate B bootstrap replicates of size rsize and apply a statistic to them.  Can do IID or Circular Block Bootstrap (CBB) methods.
}
\usage{
booter(x, statistic, B, rsize, block.length = 1, v.terms, shuffle = NULL,
    replace = TRUE, ...)
}
\arguments{
  \item{x}{
Original data series.  May be a vector, matrix or data frame object.
}
  \item{statistic}{
Function that minimally takes arguments: \code{data} and \code{\dots}.  The argument \code{data} must be the input data for which resamples are taken.  Must return a vector of all desired statistics.
}
  \item{B}{
number of bootstrap resamples to make.
}
  \item{rsize}{
Number giving the resample size for each bootstrap sample.  Must be between 1 and the length of \code{x}, if \code{x} is a vector, else the number of rows of \code{x}.  Default is to use the size of the original data.
}
  \item{block.length}{
Number giving the desired block lengths.  Default (\code{block.length} = 1) is to do IID resamples.  Should be longer than the length of dependence in the data, but much shorter than the size of the data.
}
  \item{replace}{
logical, should the resamples be taken with replacement?
}
  \item{v.terms}{
If \code{statistic} returns variance estimates for other parameters, then use this argument to specify the indices returned that give the variance estimates.  There must be a component for every other parameter returned, and they must be in the same order as the other parameters (see examples below).  If an estimate does not exist, an NA should be returned for that spot.
}
  \item{ shuffle}{ \code{rsize} by \code{B} matrix giving the indices for each bootstrap replication.  If provided, \code{B} may be missing.}
  \item{\dots}{
Optional arguments passed to \code{statistic}.
}
}
\details{
Similar functionality to \code{boot} from package \pkg{boot}, but allows for easier implementation of certain other approaches.  For example, m-out-of-n bootstrap resampling (appropriate for heavy-tail distributed data) can be performed via the \code{rsize} argument.  The \code{ci} function is used to obtain subsequent confidence limits.  For parameteric bootstrap resampling, see \code{pbooter}.

For more complicated bootstrap resampling, e.g., Bayesian bootstrap sampling, the \code{shuffle} argument may prove useful.  That is, no weighting is allowed with this function through the standard mechanism, but the same result may be obtained by supplying your own indices through the \code{shuffle} argument.  For parametric bootstrap resampling, see the \code{pbooter} function, but for certain types of parametric resampling, the \code{shuffle} argument could prove useful.

If the block length is > 1, then \code{rsize} overlapping blocks of this length are sampled from the data.  In order to minimize over or under sampling of the end points, the blocks are circular (cf. Lahiri 2003).

Many good books and other materials are available about bootstrap resampling.  One good text on IID bootstrap resampling is Efron and Tibshirani (1998) and for the block bootstrap, Lahiri (2003).
}
\value{
A list object of class \dQuote{booted} is returned with components:
  \item{call }{the function call}
  \item{data }{original data series}
  \item{statistic }{statistic argument passed in}
  \item{statistic.args}{all other arguments passed by \dots}
  \item{B}{ Number of bootstrap replicate samples}
  \item{block.length}{ The block length used}
  \item{replace}{logical stating whether the samples are taken with replacement or not. }
  \item{v.terms}{if variance terms are returned by statistic, the argument is repeated in the returned object.}
  \item{rsize}{the size of the bootstrap resamples.}
  \item{indices}{rsize by B matrix giving the resample indices used (rows) for each bootstrap resample (columns).}
  \item{v}{B length vector or B column matrix (if statistic returns a vector) giving the estimated parameter variances for each bootstrap replicate.}
  \item{orig.v}{vector giving the parameter variances (i.e. se^2) of statistic when applied to the original data.}
  \item{original.est}{vector giving the estimated parameter values when statistic is applied to the original data.}
  \item{results}{B length vector or B column matrix giving the parameter estimates for each bootstrap resample.}
  \item{type}{character stating whether the resample method is iid or cbb.}
}
\references{
Efron, B. and Tibshirani, R. J. (1998) \emph{An Introduction to the Bootstrap}.  Chapman \& Hall, Boca Raton, Florida, 436 pp.

Lahiri, S. N. (2003) \emph{Resampling Methods for Dependent Data}.  Springer-Verlag, New York, New York, 374 pp.
}
\author{
Eric Gilleland
}

\seealso{
\code{\link{pbooter}}, \code{\link{ci.booted}}
}
\examples{

z <- rnorm( 100 )

zfun <- function( data, ... ) {

    return( c( mean( data ), var( data ), mean( data^2 ), var( data^2 ) ) )

} # end of 'zfun' function.

res <- booter( x = z, statistic = zfun, B = 500, v.terms = c(2, 4) )

print( res )

ci( res )
}
\keyword{ datagen }% use one of  RShowDoc("KEYWORDS")
\keyword{ distribution }% __ONLY ONE__ keyword per line
\keyword{ nonparametric }
\keyword{ htest}
