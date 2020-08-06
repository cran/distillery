MatrixSqrt <-
function(Sigma, verbose=getOption("verbose")) {
# Find square root of 'Sigma'.  If 'Sigma' is positive definite, then
# use the Cholesky square root, otherwise use svd.

UDU <- eigen( Sigma, symmetric = TRUE )
if( is.complex( UDU$values) || is.complex( UDU$vector ) ) {

	Sigma <- 0.5*Sigma + 0.5*t( Sigma)
	UDU <- eigen( Sigma, symmetric=T)
	if( is.complex( UDU$values ) || is.complex( UDU$vector ) ) {

		warning.msg <- paste( "UDU still complex after symmetry transformation.",
			"Using only the real part for square root.",
			"Range of imaginary part is: ", " ", sep="\n")
		# warning( warning.msg)
		if( verbose ) {

			print( warning.msg )
		       	print( range( Im( UDU), na.rm=T))
			cat("\n")

		} # end of if 'verbose' stmt.

		UDU <- Re( UDU)

	} # end of inner if complex stmt
} # end of if complex stmt

if( any( UDU$values <= 0 ) || any( UDU$vectors <= 0 ) ) {

        if( verbose ) cat("Sigma is not positive-definite!")
        if( verbose) cat("\n", "Using SVD method. \n")
# A = U D U', where U are orthonormal matrices and D is diagonal matrix of
# eigenvalues.
# UDU <- eigen( Sigma, symmetric=T)
        if( is.complex( UDU$vector) || is.complex( UDU$values) && verbose )
		cat( "\nComplex Eigen Values!  Only using real part.\n" ) 
	U <- Re( UDU$vector )
        D <- Re( UDU$values )
# Make sure any eigenvalues near zero are set to zero.
# Note, if problems remain with eigenvalues less than zero, reset digits
# options.  Use ' getOption("digits") ' to see current digits used.
        D <- zapsmall( D)
# Check that eigenvalues are >= 0 and if not, stop the program.
if( any( D < 0) )
        stop("Still have eigenvalues less than zero.  Maybe reset digits?")

# Find Square Root of Sigma (call it Sigma) here.  Note that
# Sigma = U sqrt(D) U' since U are orthonormal.  Thus, we have
# (U sqrt(D) U')(U sqrt(D) U')' = U sqrt(D) U'U sqrt(D) U' = U D U' = Sigma.
# Since D is diagonal matrix, use D * t( U) for computing efficiency.
        Sigma <- U %*% ( sqrt(D)* t( U))
        } else {
# if( verbose) cat("\n", "Using Cholesky Decomposition. \n")
# Find Cholesky square root of Sigma.
        Sigma <- t( chol( Sigma))
        } # end of if else positive definiteness stmt.
return( Sigma)
}
