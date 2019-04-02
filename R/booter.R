booter <- function( x, statistic, B, rsize, block.length = 1, v.terms, shuffle = NULL, replace = TRUE, ... ) {

    theCall <- match.call()

    if( block.length < 1 ) stop( "booter: invalid block.length argument." )
    if( missing( statistic ) ) stop( "booter: must specify statistic argument." )
    if( missing( B ) && is.null( shuffle ) ) stop( "booter: must specify B or shuffle arguments." )
    if( !missing( B ) ) if( B < 1 ) stop( "booter: invalid B argument." )

    if( is.null( dim( x ) ) ) N <- length( x )
    else N <- dim( x )[ 1 ]

    if( is.null( shuffle ) ) {

        if( missing( rsize ) ) rsize <- N

        if( rsize > N ) stop( "booter: invalid rsize argument." )

        sz <- rsize / block.length

        id <- sample( 1:N, size = sz * B, replace = replace )

        if( block.length > 1 ) {

	    id <- apply( matrix( id, ncol = 1 ), 1, function( x ) return( x:( x + block.length - 1 ) ) )
	    if( any( id > N ) ) id[ id > N ] <- id[ id > N ] - N
	    id <- c( id )

        } # end of if do cbb bootstrap stmt.

        id <- matrix( id, nrow = sz * block.length, ncol = B )

    } else {

	sdim <- dim( shuffle )
	if( missing( B ) ) B <- sdim[ 2 ]

	if( missing( rsize ) ) rsize <- sdim[ 1 ]

	if( !missing( block.length ) ) warning( "booter: shuffle is provided so ignoring block.length argument." )

	id <- shuffle

    }

    bfun <- function( ind, x, statistic, ... ) {

	if( is.null( dim( x ) ) ) z <- x[ ind ]
	else z <- x[ ind, ]

	out <- do.call( statistic, c( list( data = z ), list( ... ) ) )

	return( out )

    } # end of internal 'bfun' function.

    res <- apply( id, 2, bfun, x = x, statistic = statistic, ... )
    original.est <- do.call( statistic, c( list( data = x ), list( ... ) ) )

    out <- list()
    out$call <- theCall
    out$data <- x
    out$statistic <- statistic
    out$statistic.args <- list( ... )
    out$B <- B
    out$block.length <- block.length
    out$replace <- replace
    if( !missing( v.terms ) ) out$v.terms <- v.terms
    out$rsize <- rsize
    out$indices <- id

    if( !missing( v.terms ) ) {

        out$v <- res[ v.terms, ]
        res <- res[ -v.terms, ]
	out$orig.v <- original.est[ v.terms ]
	original.est <- original.est[ -v.terms ]

    }

    out$original.est <- original.est
    out$results <- res

    if( block.length == 1 ) out$type <- "iid"
    else out$type <- "cbb"

    class( out ) <- "booted"

    return( out )

} # end of 'booter.vector' function.

ci.booted <- function( x, alpha = 0.05, ..., type = c( "perc", "basic", "stud", "bca", "norm" ) ) {

    if( missing( type ) && is.null( x$v ) ) type <- type[ -3 ]
    # else if( is.element( "stud", type ) && is.null( x$v ) ) stop( "ci: no variance to do Studentized intervals." )

    out <- list()

    conf.level <- ( 1 - alpha ) * 100

    nomen <- c( paste( conf.level, "% lower CI", sep = "" ), "Estimate",
                                paste( conf.level, "% upper CI", sep = "" ) )

    est <- x$original.est
    res <- x$results
    isv <- !is.matrix( res )

    npar <- length( est )

    if( is.element( "perc", type ) || is.element( "basic", type ) ) {

	if( isv ) {

	    perc <- c( quantile( res, probs = alpha / 2 ), est, quantile( res, probs = 1 - alpha / 2 ) )
	    basic <- c( 2 * est - perc[ 3 ], est, 2 * est - perc[ 1 ] )

	    names( perc ) <- names( basic ) <- nomen

	} else {

	    perc <- cbind( c( apply( res, 1, quantile, probs = alpha / 2 ) ), est,
				apply( res, 1, quantile, probs = 1 - alpha / 2 ) )

	    basic <- cbind( 2 * est - perc[, 3 ], est, 2 * est - perc[, 1 ] )

	    colnames( perc ) <- colnames( basic ) <- nomen

	}

	if( is.element( "perc", type ) ) {

	    attr( perc, "data.name" ) <- ""
	    attr( perc, "method" ) <- "Percentile Method"
	    attr( perc, "conf.level" ) <- conf.level
	    class( perc ) <- "ci"
	    out$perc <- perc

	}

	if( is.element( "basic", type ) ) {

	    attr( basic, "data.name" ) <- ""
            attr( basic, "method" ) <- "Basic Method"
            attr( basic, "conf.level" ) <- conf.level
            class( basic ) <- "ci"
            out$basic <- basic

        }

    } # end of if do percentile intervals stmt.

    if( is.element( "stud", type ) ) {

	v <- sqrt( x$v )
	ov <- sqrt( x$orig.v )

	if( isv ) {

	    Tstar <- ( res - est ) / v
	    stud <- c( est - ov * quantile( Tstar, probs = 1 - alpha / 2 ), est,
			est - ov * quantile( Tstar, probs = alpha / 2 ) )

	    names( stud ) <- nomen

	} else {

	    Tstar <- ( res - matrix( est, nrow = npar, ncol = x$B ) ) / v

	    stud <- cbind( est - ov * apply( Tstar, 1, quantile, probs = 1 - alpha / 2 ), est,
			    est - ov * apply( Tstar, 1, quantile, probs = alpha / 2 ) )

	    colnames( stud ) <- nomen

	}

	attr( stud, "data.name" ) <- ""
	attr( stud, "method" ) <- "Boostrap-t (Studentized) Method"
        attr( stud, "conf.level" ) <- conf.level
        class( stud ) <- "ci"

	out$stud <- stud

    } # end of if do studentized intervals stmt.


    if( is.element( "bca", type ) ) {

	if( !is.matrix( x$data ) ) N <- length( x$data )
	else N <- dim( x$data )[ 1 ]

	id <- matrix( rep( 1:N, N ), N, N )
	diag( id ) <- NA
	id <- matrix( id[ !is.na( id ) ], N - 1, N - 1 )

	zalpha <- qnorm( alpha / 2 )


	jfun <- function( ind, x, statistic, args ) {

	    if( is.null( dim( x ) ) ) z <- x[ ind ]
	    else z <- x[ ind, ]

	    res <- do.call( statistic, c( list( data = z ), args ) )

	    return( res )

	} # end of internal 'jfun' function.

	jack <- apply( id, 2, jfun, x = x$data, statistic = x$statistic, args = x$statistic.args )
	if( !is.null( x$v.terms ) ) jack <- jack[ -x$v.terms, ]

	if( isv ) {

	    z0hat <- ( sum( res < est, na.rm = TRUE ) + sum( res == est ) / 2 ) / x$B
	    z0hat <- qnorm( z0hat )
	    mth <- mean( jack, na.rm = TRUE )
	    th <- mth - jack
	    ahat <- sum( th^3, na.rm = TRUE ) / ( 6 * sqrt( sum( th^2, na.rm = TRUE )^3 ) )

	    zL <- z0hat + ( z0hat - zalpha ) / ( 1 - ahat * ( z0hat - zalpha ) )
	    zU <- z0hat + ( z0hat + zalpha ) / ( 1 - ahat * ( z0hat + zalpha ) )

	    a1 <- pnorm( zL )
	    a2 <- pnorm( zU )

	    bca <- c( quantile( res, probs = a1, na.rm = TRUE ), est, quantile( res, probs = a2, na.rm = TRUE ) )
	    names( bca ) <- nomen

	    if( !any( is.na( bca ) ) ) {

	        if( bca[ 1 ] > bca[ 3 ] || bca[ 1 ] == bca[ 3 ] ) bca[ c( 1, 3 ) ] <- NA
	        if( !is.na( bca[ 2 ] ) ) if( bca[ 1 ] > bca[ 2 ] || bca[ 2 ] > bca[ 3 ] ) bca[ c( 1, 3 ) ] <- NA

	    } else bca[ c( 1, 3 ) ] <- NA

	} else {

	    est2 <- matrix( est, nrow = npar, ncol = x$B )
	    z0hat <- ( rowSums( res < est2, na.rm = TRUE ) + rowSums( res == est2 ) / 2 ) / x$B
	    z0hat <- qnorm( z0hat )
	    mth <- apply( jack, 1, mean, na.rm = TRUE )
	    th <- mth - jack
	    ahat <- rowSums( th^3, na.rm = TRUE ) / ( 6 * sqrt( rowSums( th^2, na.rm = TRUE )^3 ) )

	    zL <- z0hat + ( z0hat + zalpha ) / ( 1 - ahat * ( z0hat + zalpha ) )
            zU <- z0hat + ( z0hat - zalpha ) / ( 1 - ahat * ( z0hat - zalpha ) )

	    a1 <- pnorm( zL )
            a2 <- pnorm( zU )

	    thL <- thU <- numeric( 0 )
	    for( i in 1:npar ) {

		thL <- c( thL, quantile( res[ i, ], probs = a1[ i ], na.rm = TRUE ) )
		thU <- c( thU, quantile( res[ i, ], probs = a2[ i ], na.rm = TRUE ) )

	    } # end of for 'i' loop.

	    bca <- cbind( thL, est, thU )
	    rownames( bca ) <- rownames( res )
	    colnames( bca ) <- nomen

	    bca.bad <- bca[, 1 ] > bca[, 3 ] | bca[, 1 ] == bca[, 3 ]  | bca[, 1 ] > bca[, 2 ] | bca[, 2 ] > bca[, 3 ]
	    if( any( bca.bad ) ) bca[ bca.bad, c( 1, 3 ) ] <- NA

	}

	attr( bca, "data.name" ) <- ""
	attr( bca, "method" ) <- "Bias-Corrected and accelerated (BCa) Percentile Method"
        attr( bca, "conf.level" ) <- conf.level
        class( bca ) <- "ci"

        out$bca <- bca
	out$bias.correction <- z0hat
	out$acceleration <- ahat

    } # end of if do basic intervals stmt.


    if( is.element( "norm", type ) ) {

	z <- qnorm( alpha / 2 )

	if( isv ) {

	    v <- sd( res, na.rm = TRUE )
	    m <- mean( res, na.rm = TRUE )
	    nrm <- c( 2 * est - m + v * z, est, 2 * est - m - v * z )
	    names( nrm ) <- c( paste( conf.level, "% lower CI", sep = "" ), "Estimate",
                                paste( conf.level, "% upper CI", sep = "" ) )

	} else {

	    v <- apply( res, 1, sd, na.rm = TRUE )
	    m <- apply( res, 1, mean, na.rm = TRUE )
	    nrm <- cbind( 2 * est - m + z * v, est, 2 * est - m - z * v )
	    colnames( nrm ) <- c( paste( conf.level, "% lower CI", sep = "" ), "Estimate",
                                paste( conf.level, "% upper CI", sep = "" ) )

	}

	attr( nrm, "data.name" ) <- ""
	attr( nrm, "method" ) <- "Normal Approximation Method"
        attr( nrm, "conf.level" ) <- conf.level
        class( nrm ) <- "ci"

        out$norm <- nrm

    } # end of if do normal approx intervals stmt.

    out$booted.object <- x

    class( out ) <- "ci.booted"

    return( out )

} # end of 'ci.booted' function.

print.booted <- function( x, ... ) { 

    print( x$call )

    if( x$type == "iid" ) {

	if( x$replace ) cat( "IID Resample Bootstrap with replacement with ", x$B, " resamples of size ", x$rsize, ".\n" )
	else cat( "IID Resample Bootstrap without replacement with ", x$B, " resamples of size ", x$rsize, ".\n" )

    } else if( x$type == "cbb" ) {

	if( x$replace ) {

	    cat( "Circular Block Bootstrap Resampling with blocks of length ", x$block.length,
	        ",\n", x$B, " resamples with replacement of size ", x$rsize, ".\n" )

	} else {

	    cat( "Circular Block Bootstrap Resampling with blocks of length ", x$block.length,
                ",\n", x$B, " resamples without replacement of size ", x$rsize, ".\n" )

	}

    } else if( x$type == "parametric" ) {

	cat( "Parametric Bootstrap with ", x$B, " resamples of size ", x$rsize, ".\n" )

    }

    invisible()

} # end of 'print.booted' function.

summary.booted <- function( object, ... ) {

    x <- object

    args <- list( ... )

    if( !is.null( args$silent ) ) silent <- args$silent
    else silent <- FALSE

    if( !is.logical( silent ) ) silent <- FALSE

    if( !silent ) print( x )

    if( is.matrix( x$results ) ) { 

	    m <- apply( x$results, 1, mean, na.rm = TRUE )
	    v <- apply( x$results, 1, var, na.rm = TRUE )

    } else {

	    m <- mean( x$results )
	    v <- var( x$results )
    }

    out <- cbind( x$original.est, x$original.est - m, v )

    if( !is.null( x$v ) ) {

	out <- cbind( out, x$orig.v )
    	colnames( out ) <- c( "Estimate", "Bias", "Variance*", "Variance.orig" )

    } else colnames( out ) <- c( "Estimate", "Bias", "Variance" )
    # else colnames( out ) <- c( "Estimate", "Bias" )

    if( !silent ) print( out )

    invisible( out )

} # end of 'summary.booted' function.

print.ci.booted <- function( x, ... ) {

    print( x$booted.object )

    if( !is.null( x$perc ) ) print( x$perc )
    if( !is.null( x$basic ) ) print( x$basic )
    if( !is.null( x$stud ) ) print( x$stud )
    if( !is.null( x$bca ) ) print( x$bca )
    if( !is.null( x$norm ) ) print( x$norm )

} # end of 'print.ci.booted' function.

