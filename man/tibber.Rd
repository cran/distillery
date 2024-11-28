\name{tibber}
\alias{tibber}
\alias{tibberRM}
\title{
Test-Inversion Bootstrap
}
\description{
Calculate (1 - alpha) * 100 percent confidence intervals for an estimated parameter using the test-inversion bootstrap method.
}
\usage{
tibber(x, statistic, B, rmodel, test.pars, rsize, block.length = 1, v.terms,
    shuffle = NULL, replace = TRUE, alpha = 0.05, verbose = FALSE, ...)

tibberRM(x, statistic, B, rmodel, startval, rsize, block.length = 1, 
    v.terms, shuffle = NULL, replace = TRUE, alpha = 0.05, step.size, 
    tol = 1e-04, max.iter = 1000, keep.iters = TRUE, verbose = FALSE, 
    ...)
}
\arguments{
  \item{x}{
numeric vector or data frame giving the original data series.
}
  \item{statistic}{
function giving the estimated parameter value.  Must minimally contain arguments \code{data} and \code{\ldots}.
}
  \item{B}{
number of replicated bootstrap samples to use.
}
  \item{rmodel}{
function that simulates data based on the nuisance parameter provided by \code{test.pars}.  Must minimally take arguments: \code{data}, \code{par}, \code{n}, and \code{\ldots}.  The first, \code{data}, is the data series (it need not be used by the function, but it must have this argument, and the original data are passed to it via this argument), \code{par} is the nuisance parameter, \code{n} is the sample size, and \code{\ldots} are any additional arguments that might be needed.
}
  \item{test.pars}{
single number or vector giving the nuisance parameter value.  If a vector of length greater than one, then the interpolation method will be applied to estimate the confidence bounds.
}
  \item{startval}{one or two numbers giving the starting value for the nuisance parameter in the Robbins-Monro algorithm.  If two numbers are given, the first is used as the starting value for the lower bound, and the second for the upper.}
  \item{rsize}{
(optional) numeric less than the length of the series given by \code{x}, used if an m-out-of-n bootstrap sampling procedure should be used.
}
  \item{block.length}{
(optional) length of blocks to use if the circular block bootstrap resampling scheme is to be used (default is iid sampling).
}
  \item{v.terms}{
(optional) gives the positions of the variance estimate in the output from \code{statistic}.  If supplied, then Studentized intervals are returned instead of (\code{tibberRM}) of in addition to (\code{tibber}) the regular intervals.  Generally, such intervals are not ideal for the test-inversion method.
}
  \item{shuffle}{
\code{n} (or \code{rsize}) by \code{B} matrix giving the indices for the resampling procedure (obviates arguments \code{block.length} and \code{B}).
}
  \item{replace}{
logical stating whether or not to sample with replacement.
}
  \item{alpha}{
significance level for the test.
}
  \item{step.size}{ Step size for the Robbins-Monro algorithm.}
  \item{tol}{ tolerance giving the value for how close the estimated p-value needs to be to \code{alpha} before stopping the Robbins-Monro algorithm.}
  \item{max.iter}{Maximum number of iterations to perform before stopping the Robbins-Monro algorithm.}
  \item{keep.iters}{logical, should information from each iteration of the Robbins-Monro algorithm be saved?}
  \item{verbose}{
logical should progress information be printed to the screen.
}
  \item{\dots}{
Optional arguments to \code{booter}, \code{statistic} and \code{rmodel}.
}
}
\details{
The test-inversion bootstrap (Carpenter 1999; Carpenter and Bithell 2000; Kabaila 1993) is a parametric bootstrap procedure that attempts to take advantage of the duality between confidence intervals and hypothesis tests in order to create bootstrap confidence intervals.  Let X = X_1,\dots,X_n be a series of random variables, T, is a parameter of interest, and R(X) is an estimator for T.  Further, let x = x_1,\ldots,x_n be an observed realization of X, and r(x) an estimate for R(X), and let x* be a bootstrap resample of x, etc.  Suppose that X is distributed according to a distribution, F, with parameter T and nuisance parameter V.

The procedure is carried out by estimating the p-value, say p*, from r*_1, ..., r*_B estimated from a simulated sample from \code{rmodel} assuming a specific value of V by way of finding the sum of r*_i < r(x) (with an additional correction for the ties r*_i = r(x)).  The procedure is repeated for each of k values of V to form a sample of p-values, p*_1, ..., p*_k.  Finally, some form of root-finding algorithm must be employed to find the values r*_L and r*_U that estimate the lower and upper values, resp., for R(X) associated with (1 - alpha) * 100 percent confidence limits.  For \code{tibber}, the routine can be executed one time if \code{test.pars} is of length one, which will enable a user to employ their own root-finding algorithm.  If \code{test.pars} is a vector, then an interpolation estimate is found for the confidence end points.  \code{tibberRM} makes successive calls to \code{tibber} and uses the Robbins-Monro algorithm (Robbins and Monro 1951) to try to find the appropriate bounds, as suggested by Garthwaite and Buckland (1992).
}
\value{
For tibber, if test.pars is of length one, then a 3 by 1 matrix is returned (or, if \code{v.terms} is supplied, then a 4 by 1 matrix) where the first two rows give estimates for R(X) based on the original simulated series and the median from the bootstrap samples, respectively.  the last row gives the estimated p-value.  If \code{v.terms} is supplied, then the fourth row gives the p-value associated with the Studentized p-value.

If test.pars is a vector with length k > 1, then a list object of class \dQuote{tibbed} is returned, which has components:
  \item{results}{3 by k matrix (or 4 by k, if \code{v.terms} is not missing) giving two estimates for R(X) (one from the simulated series and one of the median of the bootstrap resamples, resp.) and the third row giving the estimated p-value for each value of V. }
  \item{TIB.interpolated, STIB.interpolated}{numeric vector of length 3 giving the lower bound estimate, the estimate from the original data (i.e., r(x)), and the estimated upper bound as obtained from interpolating over the vector of possible values for V given by test.pars.  The Studentized TIB interval, \code{STIB.interpolated}, is only returned if \code{v.terms} is provided.}
  \item{Plow, Pup, PstudLow, PstudUp}{Estimated p-values used for interpolation of p-value. }
  \item{call}{the original function call.}
  \item{data}{the original data passed by the x argument.}
  \item{statistic, B, rmodel, test.pars, rsize, block.length, alpha, replace}{arguments passed into the orignal function call.}
  \item{n}{original sample size.}
  \item{total.time}{Total time it took for the function to run.}

For tibberRM, a list of class \dQuote{tibRMed} is returned with components:

  \item{call}{the original function call.}
  \item{x, statistic, B, rmodel, rsize, block.length, alpha, replace}{arguments passed into the orignal function call.}
  \item{result}{vector of length 3 giving the estimated confidence interval with the original parameter estimate in the second component.}
  \item{ lower.p.value, upper.p.value }{Estimated achieved p-values for the lower and upper bounds.}
  \item{ lower.nuisance.par, upper.nuisance.par }{ nuisance parameter values associated with the lower and upper bounds. }
  \item{ lower.iterations, upper.iterations}{ number of iterations of the Robbins-Monro algorithm it took to find the lower and upper bounds.}
  \item{total.time}{Total time it took for the function to run.}

}
\references{
Carpenter, James (1999) Test inversion bootstrap confidence intervals.  \emph{J. R. Statist. Soc. B}, \bold{61} (1), 159--172.

Carpenter, James and Bithell, John (2000) Bootstrap confidence intervals: when, which, what?  A practical guide for medical statisticians.  \emph{Statist. Med.}, \bold{19}, 1141--1164.

Garthwaite, P. H. and Buckland, S. T. (1992) Generating Monte Carlo confidence intervals by the Robbins-Monro process.  \emph{Appl. Statist.}, \bold{41}, 159--171.

Kabaila, Paul (1993) Some properties of profile bootstrap confidence intervals.  \emph{Austral. J. Statist.}, \bold{35} (2), 205--214.

Robbins, Herbert and Monro, Sutton (1951) A stochastic approximation method.  \emph{Ann. Math Statist.}, \bold{22} (3), 400--407.
}
\author{
Eric Gilleland
}

\seealso{
\code{\link{booter}}, \code{\link{pbooter}}
}
\examples{
# The following example follows the example provided at:
#
# http://influentialpoints.com/Training/bootstrap_confidence_intervals.htm
#
# which is provided with a creative commons license:
#
# https://creativecommons.org/licenses/by/3.0/ 
#
y <- c( 7, 7, 6, 9, 8, 7, 8, 7, 7, 7, 6, 6, 6, 8, 7, 7, 7, 7, 6, 7,
        8, 7, 7, 6, 8, 7, 8, 7, 8, 7, 7, 7, 5, 7, 7, 7, 6, 7, 8, 7, 7,
        8, 6, 9, 7, 14, 12, 10, 13, 15 )

trm <- function( data, ... ) {

    res <- try( mean( data, trim = 0.1, ... ) )
    if( class( res ) == "try-error" ) return( NA )
    else return( res )

} # end of 'trm' function.

genf <- function( data, par, n, ... ) {

    y <- data * par
    h <- 1.06 * sd( y ) / ( n^( 1 / 5 ) )
    y <- y + rnorm( rnorm( n, 0, h ) )
    y <- round( y * ( y > 0 ) )

    return( y )

} # end of 'genf' function.

look <- tibber( x = y, statistic = trm, B = 500, rmodel = genf,
    test.pars = seq( 0.85, 1.15, length.out = 100 ) )

look

plot( look )
# outer vertical blue lines should cross horizontal blue lines
# near where an estimated p-value is located.

\dontrun{
tibber( x = y, statistic = trm, B = 500, rmodel = genf, test.pars = 1 )


look2 <- tibberRM(x = y, statistic = trm, B = 500, rmodel = genf, startval = 1,
    step.size = 0.03, verbose = TRUE )

look2
# lower achieved est. p-value should be close to 0.025
# upper should be close to 0.975.

plot( look2 )

trm2 <- function( data, par, n, ... ) {

    a <- list( ... )
    res <- try( mean( data, trim = a$trim ) )
    if( class( res ) == "try-error" ) return( NA )
    else return( res )

} # end of 'trm2' function.

tibber( x = y, statistic = trm2, B = 500, rmodel = genf,
    test.pars = seq( 0.85, 1.15, length.out = 100 ), trim = 0.1 )

# Try getting the STIB interval.  v.terms = 2 below because mfun
# returns the variance of the estimated parameter in the 2nd position.
#
# Note: the STIB interval can be a bit unstable.

mfun <- function( data, ... ) return( c( mean( data ), var( data ) ) )

gennorm <- function( data, par, n, ... ) {

    return( rnorm( n = n, mean = mean( data ), sd = sqrt( par ) ) )

} # end of 'gennorm' function.

set.seed( 1544 )
z <- rnorm( 50 )
mean( z )
var( z )

# Trial-and-error is necessary to get a good result with interpolation method.
res <- tibber( x = z, statistic = mfun, B = 500, rmodel = gennorm,
    test.pars = seq( 0.95, 1.10, length.out = 100 ), v.terms = 2 )

res

plot( res )

# Much trial-and-error is necessary to get a good result with RM method.
# If it fails to converge, try increasing the tolerance.
res2 <- tibberRM( x = z, statistic = mfun, B = 500, rmodel = gennorm,
    startval = c( 0.95, 1.1 ), step.size = 0.003, tol = 0.001, v.terms = 2,
    verbose = TRUE )
# Note that it only gives the STIB interval.

res2

plot( res2 )

}
}
\keyword{ htest }
