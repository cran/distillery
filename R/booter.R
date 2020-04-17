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

} # end of 'booter' function.

plot.booted <- function( x, ..., which = 1, type = c( "ecdf", "qq", "box", "hist", "qq2" ), col = 1, distn = "norm" , distn.args = NULL, distnLongName = "Standard Normal", add.legend = TRUE ) {

	if( missing( type ) ) type <- c( "ecdf", "qq", "box", "hist" )
	else type <- match.arg( type )

	if( is.element( "qq2", type ) && length( which ) != 2 ) stop( "plot.booted: invalid which argument for qq2 option." )

	if( length( type ) > 1 ) {

		opar <- par( "mfrow" )

		if( length( type ) == 2 ) par( mfrow = c(1, 2) )
		else if( length( type ) %in% c( 3, 4 ) ) par( mfrow = c(2, 2) )
		else par( mfrow = c(2, 3) )

	} # end of set up plot display.

	if( missing( col ) ) {

		if( length( type ) == 1 ) col <- "black"
		else if( length( type ) == 2 ) col <- c( "black", "gray" )
		else col <- 1:length( col )

	} # end of if missing 'col' argument stmts.

	if( !is.null( x$parnames ) ) pnames <- x$parnames[ which ]
	else {

		pnames <- paste( "theta.hat", 1:length( which ), "*", sep = "" )
		# pnames <- character(0)
		# for( i in 1:length( which ) ) pnames <- c( pnames, expression( hat( theta )[i] ) )

	} # end of if else 'parnames' not null stmt.

	## Make empirical CDF plot if type is "ecdf".
	if( "ecdf" %in% type ) {
		
		plotBootECDF( x = x, ..., which = which, col = col,
			     distn = distn, distn.args = distn.args,
			     distnLongName = distnLongName, pnames = pnames,
			     add.legend = add.legend )

	} # end of if make ecdf plot stmt.

	if( "qq" %in% type ) {

		plotBootQQ( x = x, ..., which = which, col = col,
			     distn = distn, distn.args = distn.args,
			     distnLongName = distnLongName, pnames = pnames,
			     add.legend = add.legend )

	} # end of if make qq plot stmt.

	if( "box" %in% type ) plotBootBox( x = x, ..., which = which, col = col, pnames = pnames )
	# end of if make box plot stmt.

	if( "hist" %in% type ) {
	       
		plotBootHist( x = x, ..., which = which, col = col,
				distn = distn, distn.args = distn.args,
				distnLongName = distnLongName, pnames = pnames,
				add.legend = add.legend )

	} # end of if make histogram stmt.

	if( "qq2" %in% type ) plotBootQQ2( x = x, ..., which = which, col = col, pnames = pnames )

	if( length( type ) > 1 ) par( mfrow = opar )

	invisible()

} # end of 'plot.booted' function.

plotBootECDF <- function(  x, ..., which = 1, col, distn, distn.args, distnLongName, pnames, add.legend ) {

	if( is.matrix( x$data ) ) N <- nrow( x$data )
	else N <- length( x$data )

	z <- x$results[ which[ 1 ], ]
	Fn <- ecdf( z )

        plot( Fn, col = col[ 1 ], main = "Empirical CDF" )

	# Add the parametric distribution if desired.
	if( !is.na( distn ) && distn != "" ) {

		pdistn <- paste( "p", distn, sep = "" )
		qdistn <- paste( "q", distn, sep = "" )
		zq <- do.call( qdistn, args = c( list( p = ( 1 - 0.5 ) / N ), distn.args ) )
		# zd <- approx( seq( 0, 1,, 100 ), zq, n = N )$y
		zp <- do.call( pdistn, args = c( list( q = zq ), distn.args ) )
		lines( zq, zp, lty = 2, ... )
		theory.message <- paste( distnLongName, "CDF" )

	} else theory.message <- NULL
	# end of add theoretical distribution stmt.

        if( length( which ) > 1 ) {

		Fn2 <- list()
		if( length( col ) == 1 ) col <- rep( col, length( which ) )

	        for( i in 2:length( which ) ) {

			Fn2[[ i - 1 ]] <- ecdf( x$results[ which[ i ], ] )
			plot( Fn2[[ i - 1 ]], add = TRUE, col = col[ i ] )

		} # end of for 'i' loop.

	} # end of if 'which > 1' stmt

	pnames <- c( pnames, theory.message )
	if( !is.null( theory.message ) ) {

	       	col <- c( col, "black" )
		ltype <- c( rep( 1, length( which ) ), 2 )

	} else ltype <- 1
	if( add.legend ) legend( "topleft", legend = pnames, col = col, lty = ltype, bty = "n" )

	invisible()

} # end of 'plotBootECDF' function. 

plotBootQQ <- function( x, ..., which, col, distn, distn.args, distnLongName, pnames, add.legend ) {

	z <- x$results[ which[ 1 ], ]
	N <- length( z )
	qdistn <- paste( "q", distn, sep = "" )
	zq <- do.call( qdistn, args = c( list( p = ( 1:N - 0.5 ) / N ), distn.args ) )

	qqplot( zq, z, xlab = paste( distnLongName, "Quantiles" ), ylab = "Bootstrap Sample Quantiles", col = col[ 1 ] )

	if( length( which ) > 1 ) {

		if( length( col ) == 1 ) col <- c( col, 2:length( which ) )

		for( i in 2:length( which ) ) {

			newq <- qqplot( x = zq, y = x$results[ which[ i ], ], plot.it = FALSE )
			points( newq$x, newq$y, col = col[ i ], pch = "+" )

		} # end of for 'i' loop.

		if( add.legend ) legend( "topleft", legend = pnames, col = col, bty = "n" )

	} # end of if add more quantiles to plot statement.

	invisible()

} # end of 'plotBootQQ' function.

plotBootBox <- function( x, ..., which, col, pnames ) {

	y <- as.data.frame( t( x$results[ which, , drop = FALSE ] ) )
	colnames( y ) <- pnames

	boxplot( y, col = col, ... )

	invisible()

} # end of 'plotBootBox' function.

plotBootHist <- function( x, ..., which, col, distn, distn.args, distnLongName, pnames, add.legend ) {

	if( length( which ) > 1 )
		warning( "plotBootHist: which argument has length > 1. Only first is used in the histogram." )

	y <- x$results[ which[ 1 ], ]
	hist( y, main = paste( "Histogram of ", pnames[ 1 ], sep = "" ), xlab = pnames[ 1 ], freq = FALSE )

	if( !is.na( distn ) && distn != "" ) {

		r <- range( y, finite = TRUE )
		z <- seq( r[ 1 ], r[ 2 ], , 200 )
		ddistn <- paste( "d", distn, sep = "" )
		zd <- do.call( ddistn, args = c( list( x = z ), distn.args ) )
		lines( z, zd, ... )
		if( add.legend ) legend( "topleft", legend = paste( distnLongName, "PDF" ), lty = 1, bty = "n" )

	} # end of if add theoretical density stmt.

	invisible()

} # end of 'plotBootHist' function.

plotBootQQ2 <- function( x, ..., which, col, pnames ) {

	if( length( which ) != 2 ) stop( "plotBootQQ2: invalid which argument.  Must have length 2." )
	if( length( col ) > 1 ) col <- col[ 1 ]

	y <- x$results[ which[ 1 ], ]
	z <- x$results[ which[ 2 ], ]

	qqplot( y, z, xlab = pnames[ 1 ], ylab = pnames[ 2 ] )

	invisible()

} # end of 'plotBootQQ2' function.

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

