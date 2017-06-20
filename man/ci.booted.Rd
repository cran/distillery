\name{ci.booted}
\alias{ci.booted}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bootstrap Confidence Intervals
}
\description{
Calculate confidence intervals for objects output from the booter function.
}
\usage{
\method{ci}{booted}(x, alpha = 0.05, ..., type = c("perc", "basic", "stud", "bca", "norm"))
}
\arguments{
  \item{x}{
object of class \dQuote{booted} as returned by the \code{booter} function.
}
  \item{alpha}{
Significance level for which the (1 - alpha) * 100 percent confidence intervals are determined.
}
  \item{\dots}{
Not used.
}
  \item{type}{
character stating which intervals are to be reutrned.  Default will do them all.
}
}
\details{
Many methods exist for sampling parameters associated with a data set, and many methods for calculating confidence intervals from those resamples are also available.  Some points to consider when using these methods are the accuracy of the intervals, and whether or not they are range-preserving and/or transformation-respecting.  An interval that is range-preserving means that if a parameter can only take on values within a specified range, then the end points of the interval will also fall within this range.  Transformation-respecting means that if a parameter, say phi, is transformed by a monotone function, say m(phi), then the (1 - alpha) * 100 percent confidence interval for m(phi) can be derived by applying m() to the limits of the (1 - alpha) * 100 percent interval for phi.  That is [L(phi), U(phi)] = [m(L(phi)), m(U(phi))].

For accuracy, a (1 - 2 * alpha) * 100 percent confidence interval, (L, U), is presumed to have probability alpha of not covering the true value of the parameter from above or below.  That is, if theta is the true value of the parameter, then Pr( theta < L ) = alpha, and Pr( theta > U ) = alpha.  A second-order accurate interval means that the error in these probabilities tends to zero at a rate that is inversely proportional to the sample size.  On the other hand, first-order accuracy means that the error tends to zero more slowly, at a rate inversely proportional to the square root of the sample size.

the types of intervals available, here, are described below along with some considerations for their use.

Percentile intervals (\code{type} = \dQuote{perc}) are 1st order accurate, range-preserving, and transformation-respecting.  However, they may have poor coverage in some situations.  They are given by (L, U) where L and U are the 1 - alpha / 2 and alpha / 2 quantiles of the non-parametric distribution obtained through bootstrap resampling.

The basic interval (\code{type} = \dQuote{basic}) is the originally proposed interval and is given by (2 * theta - U, 2 * theta - L ), where U and L are as for the percentile interval.  This interval is 1st order accurate, but is not range-preserving or transformation-respecting.

Studentized (or Bootstrap-t) intervals (\code{type} = \dQuote{stud}) are 2nd order accurate, but not range-preserving or transformation-respecting, and they can be erratic for small samples, as well as sensitive to outliers.  They are obtained by the basic bootstrap, but where U and L are taken from the studentized version of the resampled parameter estimates.  That is, T' is taken for each bootstrap replicate, b, to be:

T'(b) = (theta'(b) - theta) / (se'(b)), where theta'(b) and se'(b) are the estimated value of the parameter and its estimated standard error, resp., for bootstrap replicate b, and theta is the estimated parameter value using the original data.

The bias-corrected and accelerated (BCa, \code{type} = \dQuote{bca}) method applies a bias correction and adjustment to the percentile intervals.  The intervals are 2nd order accurate, range-preserving and transformation-respecting.  However, the estimation performed, here (Eq 14.15 in Efron and Tibshirani 1998), requires a further jacknife resampling estimation, so the computational burden can be more expensive.  The estimates for the bias-correction and acceleration adjustment can be found in Efron and Tibshirani (1998) p. 178 to 201.  The bias-correction factor includes an adjustment for ties.  

Finally, the normal approximation interval (\code{type} = \dQuote{norm}) uses the average of the estimated parameters from the bootstrap replicates, call it m, and their standard deviation, call is s, to make the usual normal approximation interval.  An assumption of normality for the parameter estimates is assumed, which means that they will be symmetric.  This method yields 1st order accurate intervals that are not range-preserving or transformation-respecting.
}
\value{
A list object of class \dQuote{ci.booted} is returned with components depending on which types of intervals are calculated.
  \item{booted.object }{The object passed through the x argument.}
  \item{perc, basic, stud, bca, norm }{vectors of length 3 or 3-column matrices giving the intervals and original parameter estimates for each CI method.}
  \item{bias.correction, accelerated}{If type includes \dQuote{bca}, then the estiamted bias correction factor and acceleration are given in these components. }
}
\references{
Efron, B. and Tibshirani, R. J. (1998) \emph{An Introduction to the Bootstrap}.  Chapman \& Hall, Boca Raton, Florida, 436 pp.
}
\author{
Eric Gilleland
}

\seealso{
\code{\link{booter}}, \code{\link{pbooter}}
}
\examples{
##
## See the help file for booter and/or pbooter for examples.
##
}
\keyword{ htest }% use one of  RShowDoc("KEYWORDS")
