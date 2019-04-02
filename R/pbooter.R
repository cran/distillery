pbooter <- function( x, statistic, B, rmodel, rsize, v.terms, verbose = FALSE, ... ) {

    begin.tiid <- Sys.time()

    theCall <- match.call()

    if( missing( statistic ) ) stop( "booter: must specify statistic argument." )
    if( missing( B ) ) stop( "booter: must specify B argument." )
    if( B < 1 ) stop( "booter: invalid B argument." )

    # if( missing( rsize ) ) stop( "pbooter: must specify rsize argument." )
    # if( rsize < 1 ) stop( "pbooter: invalid rsize argument." )

    # Determine if 'x' is a vector or a matrix.  Also, specify 'rsize' if it is not
    # specified by the user.
    xdim <- dim( x )
    if( is.null( xdim ) ) {

	isv <- TRUE
	if( missing( rsize ) ) rsize <- length( x )

    } else {

	isv <- FALSE
	if( missing( rsize ) ) rsize <- xdim[ 1 ]

    } # end of if else 'x' is a vector a matrix stmt.

    if( rsize < 1 ) stop( "pbooter: invalid rsize argument." )

    bfun <- function( x, statistic, verbose, ... ) {

        return( do.call( statistic, c( list( data = x ), list( ... ) ) ) )

    } # end of 'bfun' function.

    if( verbose ) cat( "Simulating data from ", deparse( substitute( rmodel ) ), "\n" )
    xdat <- do.call( rmodel, c( list( size = rsize * B ), list( ... ) ) )

    if( verbose ) {

	cat( "Simulations found.  Calculating statistic for each simulation.\n" )

    } # end of if 'verbose' stmt.

    if( isv ) {

        xdat <- array( xdat, dim = c( rsize, xdim[ 2 ], B ) )
        res <- apply( xdat, 2, bfun, statistic = statistic, verbose = verbose, ... )

    } else {

        hold <- array( dim = c( rsize, xdim[ 2 ], B ) )
	bb <- rep( 1:B, each = rsize )
        for( b in 1:B ) hold[,, b ] <- xdat[ bb == b, ]
        xdat <- hold
        res <- apply( xdat, 3, bfun, statistic = statistic, ... )

    } # end of does 'rmodel' return a vector or a matrix stmts.

    if( verbose ) cat( "Bootstrap samples obtained.  Calculating original estimates.\n" )

    original.est <- do.call( statistic, c( list( data = x ), list( ... ) ) )

    out <- list()
    out$call <- theCall
    out$data <- x
    out$statistic <- statistic
    out$statistic.args <- list( ... )
    out$B <- B
    if( !missing( v.terms ) ) out$v.terms <- v.terms
    out$rsize <- rsize
    out$rdata <- xdat

    if( !missing( v.terms ) ) {

	out$v <- res[ v.terms, ]
        res <- res[ -v.terms, ]
        out$orig.v <- original.est[ v.terms ]
        original.est <- original.est[ -v.terms ]

    } # end of if variance terms calculated or not stmt.

    out$original.est <- original.est
    out$results <- res

    out$process.time <- Sys.time() - begin.tiid
    if( verbose ) print( out$process.time )

    out$type <- "parametric"
    class( out ) <- "booted"

    return( out )

} # end of 'pbooter' function.
