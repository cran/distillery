\name{MatrixSqrt}
\alias{MatrixSqrt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Square-Root of a Square Matrix
}
\description{
	Find the (approximate ) square-root of a square matrix that is possibly not positive definite using the singular-value decomposition. 
}
\usage{
MatrixSqrt( Sigma, verbose = getOption("verbose") )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Sigma}{
matrix for which the square root is to be taken.
}
  \item{verbose}{
	  logical, should progress information be printed to the screen.
}
}
\details{
	The \code{eigen} function is first called in order to obtain the eigen values and vectors.  If any are complex then a symmetry transformation is applied (i.e., Sigma = 0.5 * ( Sigma + t( Sigma ) ) ) and then the \code{eigen} function is called again.  Eigen values that are less than zero, but close to zero, are set to zero.  If the matrix is positive definite, then the \code{chol} function is called in order to return the Cholesky decomposition.  Otherwise, U sqrt( D ) U' is returned, where U is the matrix of eigen vectors and D a diagonal matrix whose diagonal contains the eigen values.  The function will try to find the square root even if it is not positive definite, but it may fail.
}
\value{
	A matrix is returned.
}
\references{
	Hocking, R. R. (1996) \emph{Methods and Applications of Linear Models}.  Wiley Series in Probability and Statistics, New York, NY, 731 pp.
}
\author{
	Eric Gilleland
}

\seealso{
\code{\link{eigen}}, \code{\link{chol}}
}
\examples{
# Simulate 3 random variables, Y, X1 and X2, such that
# Y is correlated with both X1 and X2, but X1 and X2
# are uncorrelated.

set.seed( 2421 );

Z <- matrix( rnorm( 300 ), 100, 3 );
R1 <- cbind( c( 1, 0.8, 0.6 ), c( 0.8, 1, 0 ), c( 0.6, 0, 1 ) );
R2 <- MatrixSqrt( R1 );

# R1;
# R2 %*% t( R2 );
# zapsmall( R2 %*% t( R2 ) );

Z <- Z %*% R2;
Y <- Z[,1];
X1 <- Z[,2];
X2 <- Z[,3];
cor( Y, X1 );
cor( Y, X2 );
cor( X1, X2 );
plot( Y, X1, pch = 20, col = "darkblue",
     bg = "darkblue", cex = 1.5 );
points( Y, X2, col = "darkgray", pch = "+", cex = 1.5 );
plot( X1, X2 );

\dontrun{
# The following line will give an error message.
# chol( R1 );
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ array }% use one of  RShowDoc("KEYWORDS")
\keyword{ algebra }% __ONLY ONE__ keyword per line
