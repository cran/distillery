tibber <- function( x, statistic, B, rmodel, test.pars, rsize, block.length = 1, v.terms, 
    shuffle = NULL, replace = TRUE, alpha = 0.05, verbose = FALSE, ... ) {

    begin.tiid <- Sys.time()

    theCall <- match.call()

    out <- list()

    #
    # 'rmodel' is a random generator function that minimally takes arguments:
    #	 'data', 'par', 'n', and '...', where 'par' is a single parameter value
    #	 'n' is the sample size to draw.
    #

    if( missing( x ) ) stop( "tibber: Must supply x argument." )
    if( missing( statistic ) ) stop( "tibber: Must supply statistic argument." )
    if( missing( B ) ) stop( "tibber: Must supply B argument." )
    if( missing( rmodel ) ) stop( "tibber: Must supply rmodel argument." )
    if( missing( test.pars ) ) stop( "tibber: test.pars argument not supplied." )

    nT <- length( test.pars )

    if( missing( rsize ) ) {

	if( is.null( dim( x ) ) ) n <- length( x )
	else n <- dim( x )[ 1 ]

	rsize <- n

    } else n <- rsize

    obs <- do.call( statistic, c( list( data = x ), list( ... ) ) )

    if( !missing( v.terms ) ) {

	obs.v <- obs[ v.terms ]
	obs <- obs[ -v.terms ]

    } # end of 'v.terms' missing or not stmts.

    tfun <- function( tp, x, st, B, rmodel, N, rsize, bl, v.terms, 
		shuffle, replace, o, osd, ... ) {

	xstar <- do.call( rmodel, c( list( data = x, par = tp, n = N ), list( ... ) ) )

	if( verbose ) cat( "Nuisance parameter = ", tp, "\n" )

	if( missing( rsize ) && missing( v.terms ) ) {
	
	    bootobj <- booter( x = xstar, statistic = st, B = B, block.length = bl,
			    shuffle = shuffle, replace = replace, ... )

	} else if( missing( rsize ) ) {

	    bootobj <- booter( x = xstar, statistic = st, B = B, block.length = bl,
                            v.terms = v.terms, shuffle = shuffle, replace = replace, ... )

	} else if( missing( v.terms ) ) {

	    bootobj <- booter( x = xstar, statistic = st, B = B, rsize = rsize, block.length = bl,
                            shuffle = shuffle, replace = replace, ... )

	} else {

	    bootobj <- booter( x = xstar, statistic = st, B = B, rsize = rsize, block.length = bl,
                            v.terms = v.terms, shuffle = shuffle, replace = replace, ... )

	}

	# Qk <- do.call( st, c( list( data = xstar[ c( bootobj$indices ) ], list( ... ) ) ) )
	Qk <- bootobj$original.est

	if( !missing( v.terms ) ) Qk.v <- bootobj$orig.v

	# TIB is not defined for more than one parameter at a time, but for now allowing
	# for the possibility.
	if( is.null( dim( bootobj$results ) ) ) {

	    B2 <- sum( !is.na( bootobj$results ), na.rm = TRUE )

	    Qmk <- quantile( bootobj$results, probs = 0.5 )
	    Pk  <- ( sum( bootobj$results < o, na.rm = TRUE ) + sum( bootobj$results == o, na.rm = TRUE ) / 2 ) / B2

	    res <- c( Qk, Qmk, Pk )

	    if( !missing( v.terms ) )  {

		Qmk.sd <- sqrt( bootobj$v )
		Tstar <- ( bootobj$results - Qk ) / Qmk.sd
		ot <- ( o - Qk ) / osd
		
		Pk.stud <- ( sum( Tstar < ot, na.rm = TRUE ) + sum( Tstar == ot, na.rm = TRUE ) / 2 ) / B2

		# Qk.stud <- quantile( Tstar, probs = 0.5 )

		res <- c( res, Pk.stud )

	    } # end of if do studentized intervals too stmt.

	} else {

	    B2 <- colSums( !is.na( bootobj$results ), na.rm = TRUE )

	    Qmk <- apply( bootobj$results, 2, quantile, probs = 0.5 ) 
	    bigO <- matrix( o, nrow = nrow( bootobj$results ), ncol = B )
	    o1  <- bootobj$results < bigO
	    o2  <- bootobj$results == bigO
	    Pk  <- rowSums( o1 + o2 / 2, na.rm = TRUE ) / B2

	    res <- c( Qk, Qmk, Pk )

	    if( !missing( v.terms ) ) {

		Qmk.sd <- matrix( sqrt( bootobj$v ), nrow = nrow( bootobj$results ), ncol = B )
		Qkmat <- matrix( Qk, nrow = nrow( bootobj$results ), ncol = B )
		Tstar <- ( bootobj$results - Qkmat ) / Qmk.sd
		ot <- matrix( ( o - Qk ) / osd, nrow = nrow( bootobj$results ), ncol = B )

		o1 <- Tstar < ot
		o2 <- Tstar == ot

		Pk.stud <- ( rowSums( o1 + o2 / 2, na.rm = TRUE ) ) / B2

		# Qk.stud <- apply( Tstar, 1, quantile, probs = 0.5 )

		res <- c( res, Pk.stud )

	    } # end of if do studentized intervals too stmt.

	}

	return( res )

    } # end of internal 'tfun' function.

    if( missing( rsize ) && missing( v.terms ) ) {

        res <- apply( cbind( test.pars ), 1, tfun, x = x, st = statistic, B = B,
		rmodel = rmodel, N = n, bl = block.length, shuffle = shuffle,
	 	replace = replace, o = obs, ... )

    } else if( missing( rsize ) ) {

	res <- apply( cbind( test.pars ), 1, tfun, x = x, st = statistic, B = B,
                v.terms = v.terms, rmodel = rmodel, N = n, bl = block.length, shuffle = shuffle,
                replace = replace, o = obs, osd = sqrt( obs.v ), ... )

    } else if( missing( v.terms ) ) {

	res <- apply( cbind( test.pars ), 1, tfun, x = x, st = statistic, B = B, rsize = rsize,
                rmodel = rmodel, N = n, bl = block.length, shuffle = shuffle,
                replace = replace, o = obs, ... )

    } else {

	res <- apply( cbind( test.pars ), 1, tfun, x = x, st = statistic, B = B, rsize = rsize,
                v.terms = v.terms, rmodel = rmodel, N = n, bl = block.length, shuffle = shuffle,
                replace = replace, o = obs, osd = sqrt( obs.v ), ... )

    }

    if( missing( v.terms ) ) rownames( res ) <- c( "Est. Parameter", "Bootstrap Median Est. Parameter", "Est. P-value" )
    else rownames( res ) <- c( "Est. Parameter", "Bootstrap Median Est. Parameter", "Est. P-value", "Est. P-value (Stud.)" )

    colnames( res ) <- test.pars

    out$results <- res

    if( nT > 1 ) { # Do interpolation method.

        # 'res' is a 3 (or 4 if v.terms present) by k matrix.
        P <- res[ 3, ]
        Q <- res[ 1, ]
    
        a2 <- alpha / 2
    
        q1 <- min( which( P < 1 - a2 ), na.rm = TRUE )
        q2 <- max( which( P >= 1 - a2 ), na.rm = TRUE )
    
        P1 <- P[ q1 ]
        P2 <- P[ q2 ]
    
        Q1 <- Q[ q1 ]
        Q2 <- Q[ q2 ]
    
        L <- ( Q1 - Q2 ) * ( P2 - ( 1 - a2 ) ) / ( P2 - P1 ) + Q2
    
        q1 <- min( which( P <= a2 ), na.rm = TRUE )
        q2 <- max( which( P > a2 ), na.rm = TRUE )
    
        P1 <- P[ q1 ]
        P2 <- P[ q2 ]
    
        Q1 <- Q[ q1 ]
        Q2 <- Q[ q2 ]
    
        U <- ( Q1 - Q2 ) * ( P2 - a2 ) / ( P2 - P1 ) + Q2
    
        # Note: if it is possible to have multiple statistics, then the following will need modification.
        tib <- c( L, obs, U )
        names( tib ) <- c( "Lower", "Estimate", "Upper" )
    
        out$TIB.interpolated <- tib
    
        if( !missing( v.terms ) ) {
    
    	    Pstud <- res[ 4, ]
    	    # Qstud <- res[ 5, ]
    
    	    q1stud <- min( which( Pstud < 1 - a2 ), na.rm = TRUE )
    	    q2stud <- max( which( Pstud >= 1 - a2 ), na.rm = TRUE )
    
    	    P1stud <- Pstud[ q1stud ]
    	    P2stud <- Pstud[ q2stud ]
    
    	    Q1 <- Q[ q1stud ]
    	    Q2 <- Q[ q2stud ]
    
    	    Lstud <- ( Q1 - Q2 ) * ( P2stud - ( 1 - a2 ) ) / ( P2stud - P1stud ) + Q2
    
    	    q1stud <- min( which( Pstud <= a2 ), na.rm = TRUE )
            q2stud <- max( which( Pstud > a2 ), na.rm = TRUE )
    
    	    P1stud <- Pstud[ q1stud ]
            P2stud <- Pstud[ q2stud ]
    
    	    Q1 <- Q[ q1stud ]
            Q2 <- Q[ q2stud ]
    
    	    Ustud <- ( Q1 - Q2 ) * ( P2stud - a2 ) / ( P2stud - P1stud ) + Q2
   
	    stib <- c( Lstud, obs, Ustud ) 
	    names( stib ) <- c( "Lower", "Estimate", "Upper" )
    	    out$STIB.interpolated <- stib
    
        } # end of if do STIB interval stmt.
    
        out$call <- theCall
    
        out$data <- x
        out$statistic <- statistic
        out$B <- B
        out$rmodel <- rmodel
        out$n <- n
        out$test.pars <- test.pars
        out$rsize <- rsize
        out$block.length <- block.length
        if( !missing( v.terms ) ) out$v.terms <- v.terms
        out$shuffle <- shuffle
        out$replace <- replace
        out$alpha <- alpha
    
        tot.tiid <- Sys.time() - begin.tiid
        out$total.time <- tot.tiid
        if( verbose ) print( tot.tiid )
    
        class( out ) <- "tibbed"

	return( out )

    } else return( res )

} # end of 'tibber' function.

tibberRM <- function( x, statistic, B, rmodel, startval, rsize, block.length = 1, v.terms,
    shuffle = NULL, replace = TRUE, alpha = 0.05, step.size, tol = 1e-4, max.iter = 1000,
    keep.iters = TRUE, verbose = FALSE, ... ) {

    if( missing( startval ) ) stop( "tibberRM: must specify startval." )

    begin.tiid <- Sys.time()

    theCall <- match.call()

    out <- list()

    out$call <- theCall

    out$x <- x
    out$statistic <- statistic
    out$B <- B
    out$rmodel <- rmodel
    if( !missing( v.terms ) ) out$v.terms <- v.terms
    out$shuffle <- shuffle
    out$alpha <- alpha
    out$step.size <- step.size
    out$tolerance <- tol
    out$extra.args <- list( ... )
    out$block.length <- block.length
    out$replace <- replace

    if( missing( rsize ) ) {

        if( is.null( dim( x ) ) ) n <- length( x )
        else n <- dim( x )[ 1 ]

        rsize <- n

    } else n <- rsize

    out$rsize <- rsize

    # Specify starting values for the nuisance parameter.
    if( length( startval ) == 1 ) Lstart <- Ustart <- startval
    else if( length( startval ) == 2 ) {

	Lstart <- startval[ 1 ]
	Ustart <- startval[ 2 ]

    } else stop( "tibberRM: invalid startval argument." )

    a2 <- alpha / 2

    if( verbose ) {

	cat( "\nBeginning to find lower bound.\n" )
	cat( "Obtaining bootstrap estimated p-value for nuisance parameter = ", Lstart, "\n" )

    } # end of if 'verbose' stmt.

    # Do lower end-point.
    if( missing( v.terms ) ) {

        theta.star <- tibber( x = x, statistic = statistic, B = B, rmodel = rmodel,
			test.pars = Lstart, rsize = rsize, block.length = block.length,
			shuffle = shuffle, replace = replace, alpha = alpha,
			verbose = verbose, ... )

    } else {

	theta.star <- tibber( x = x, statistic = statistic, B = B, rmodel = rmodel,
                        test.pars = Lstart, rsize = rsize, block.length = block.length, v.terms = v.terms,
			shuffle = shuffle, replace = replace, alpha = alpha,
			verbose = verbose, ... )

    } # end of find starting value for parameter of lower end-point.

    par.star  <- theta.star[ 1 ]
    if( missing( v.terms ) ) pval.star <- theta.star[ 3 ]
    else pval.star <- theta.star[ 4 ]

    dP <- abs( pval.star - a2 )

    th.star <- Lstart

    if( keep.iters ) lower.iters <- c( par.star, th.star, pval.star, dP )

    i <- 0

    if( verbose ) cat( "Initial estimate (p-value) = ", par.star, "(", pval.star, ")\n" )

    while( ( dP > tol ) && ( i < max.iter ) ) {

	i <- i + 1
	th.n <- th.star - step.size * ( ( 1 - a2 ) - pval.star ) / i

	# if( verbose ) cat( "Obtaining bootstrap estimated p-value for nuisance parameter = ", th.n, "\n" )

	if( missing( v.terms ) ) {

            theta.star <- tibber( x = x, statistic = statistic, B = B, rmodel = rmodel,
                        test.pars = th.n, rsize = rsize, block.length = block.length, shuffle = shuffle,
                        replace = replace, alpha = alpha, verbose = verbose, ... )

        } else {

            theta.star <- tibber( x = x, statistic = statistic, B = B, rmodel = rmodel,
                        test.pars = th.n, rsize = rsize, block.length = block.length, v.terms = v.terms,
                        shuffle = shuffle, replace = replace, alpha = alpha,
                        verbose = verbose, ... )

        } 

	par.star  <- theta.star[ 1 ]
        if( missing( v.terms ) ) pval.star <- theta.star[ 3 ]
	else pval.star <- theta.star[ 4 ]

	th.star <- th.n

	if( keep.iters ) lower.iters <- rbind( lower.iters, c( par.star, th.star, pval.star, abs( pval.star - ( 1 - a2 ) ) ) )

        dP <- abs( pval.star - ( 1 - a2 ) )

    } # end of 'while' loop.

    if( dP > abs( pval.star - a2 ) ) warning( "tibberRM: did not converge for lower end-point." )
    if( i >= max.iter - 1 ) warning( "tibberRM: iterations reached max.iter for lower end-point." )

    L <- par.star

    if( keep.iters ) {

	if( is.matrix( lower.iters ) ) colnames( lower.iters ) <- c( "parameter", "nuisance parameter", "p-value", "|p-value - (1 - alpha/2)|" )
	else names( lower.iters ) <- c( "parameter", "nuisance parameter", "p-value", "|p-value - (1 - alpha/2)|" )
	out$lower.iters <- lower.iters

    }
    out$lower.p.value <- 1 - pval.star
    out$lower.nuisance.par <- th.star
    out$lower.iterations <- i

    if( verbose ) cat( "\nLower bound estimate = ", L, "\n\n" )

    # Do upper end-point.
    if( verbose ) {

        cat( "\nNow find upper bound.\n" )
        cat( "Obtaining bootstrap estimated p-value for nuisance parameter = ", Ustart, "\n" )

    } # end of if 'verbose' stmt.

    if( missing( v.terms ) ) {

        theta.star <- tibber( x = x, statistic = statistic, B = B, rmodel = rmodel,
                        test.pars = Ustart, rsize = rsize, block.length = block.length, shuffle = shuffle,
                        replace = replace, alpha = alpha, verbose = verbose, ... )

    } else {

        theta.star <- tibber( x = x, statistic = statistic, B = B, rmodel = rmodel,
                        test.pars = Ustart, rsize = rsize, block.length = block.length, v.terms = v.terms,
                        shuffle = shuffle, replace = replace, alpha = alpha,
                        verbose = verbose, ... )

    } # end of find starting value for parameter of upper end-point.

    par.star  <- theta.star[ 1 ]
    if( missing( v.terms ) ) pval.star <- theta.star[ 3 ]
    else pval.star <- theta.star[ 4 ]

    dP <- abs( pval.star - ( 1 - a2 ) )
    th.star <- Ustart

    if( keep.iters ) upper.iters <- c( par.star, th.star, pval.star, dP )

    i <- 0

    while( ( dP > tol ) && ( i < max.iter ) ) {

	i <- i + 1
	th.n <- th.star - step.size * ( a2 - pval.star ) / i

	# if( verbose ) cat( "Obtaining bootstrap estimated p-value for nuisance parameter = ", th.n, "\n" )

        if( missing( v.terms ) ) {

            theta.star <- tibber( x = x, statistic = statistic, B = B, rmodel = rmodel,
                        test.pars = th.n, rsize = rsize, block.length = block.length, shuffle = shuffle,
                        replace = replace, alpha = alpha, verbose = verbose, ... )

        } else {

            theta.star <- tibber( x = x, statistic = statistic, B = B, rmodel = rmodel,
                        test.pars = th.n, rsize = rsize, block.length = block.length, v.terms = v.terms,
                        shuffle = shuffle, replace = replace, alpha = alpha,
                        verbose = verbose, ... )
        
        }

	par.star  <- theta.star[ 1 ]
        if( missing( v.terms ) ) pval.star <- theta.star[ 3 ]
        else pval.star <- theta.star[ 4 ]

        th.star <- th.n

	if( keep.iters ) upper.iters <- rbind( upper.iters, c( par.star, th.star, pval.star, abs( pval.star - a2 ) ) )

        dP <- abs( pval.star - a2 )

    } # end of (second) 'while' loop.

    if( dP > abs( pval.star - ( 1 - a2 ) ) ) warning( "tibberRM: did not converge for upper end-point." )
    if( i >= max.iter - 1 ) warning( "tibberRM: iterations reached max.iter for upper end-point." )

    U <- par.star

    if( keep.iters ) {

        if( is.matrix( upper.iters ) ) colnames( upper.iters ) <- c( "parameter", "nuisance parameter", "p-value", "|p-value - alpha/2|" )
	else names( upper.iters ) <- c( "parameter", "nuisance parameter", "p-value", "|p-value - alpha/2|" )

	out$upper.iters <- upper.iters

    }

    out$upper.p.value <- 1 - pval.star
    out$upper.nuisance.par <- th.star
    out$upper.iterations <- i

    if( verbose ) cat( "\nUpper bound estimate = ", U, "\n\n" )

    obs <- do.call( statistic, c( list( data = x ), list( ... ) ) )

    if( !missing( v.terms ) ) {

	obs.v <- obs[ v.terms ]
	obs <- obs[ -v.terms ]

	out$type <- "STIB" 
	# warning( "tibberRM: v.terms specified, but STIB not yet supported with Robbins-Monro algorithm." )

    } else out$type <- "TIB"

    res <- c( L, obs, U )
    names( res ) <- c( "Lower", "Estimate", "Upper" )

    out$result <- res

    class( out ) <- "tibRMed"
    tot.tiid <- Sys.time() - begin.tiid
    out$total.time <- tot.tiid

    return( out )

} # end of 'tibberRM' function.

print.tibbed <- function( x, ... ) {

    cat( "\nCall:\n" )
    print( x$call )

    cat( "\nTest-inversion Bootstrap ", ( 1 - x$alpha ) * 100, "% CI (interpolated)\n" )
    print( x$TIB )

    if( !is.null( x$STIB ) ) {

	cat( "Studentized version:\n" )
	print( x$STIB )

    }

    invisible()

} # end of 'print.tibbed' function.

plot.tibbed <- function( x, ..., type = c("pvalue", "estimates") ) {

    type <- match.arg( type )

    theta1 <- x$results[ 1, ]
    if( type == "pvalue" ) {

	y  <- x$results[ 3, ]
	xl <- "Bootstrap estimated p-value"

    } else {

	y <- x$results[ 2, ]
	xl <- "Median Parameter Estimate from Bootstrap Resamples"

    }

    if( type == "pvalue" && dim( x$results )[ 1 ] == 4 ) {

	y2 <- x$results[ 4, ]
	leg <- TRUE

    } else leg <- FALSE

    a <- list( ... )

    if( is.null( a$xlab ) && is.null( a$ylab ) ) plot( theta1, y, xlab = xl, ylab = "Est. Parameter", ... )
    else if( is.null( a$xlab ) ) plot( theta1, y, xlab = xl, ... )
    else if( is.null( a$ylab ) ) plot( theta1, y, ylab = "Est. Parameter", ... )
    else plot( theta1, y, ... )

    if( type == "pvalue" ) abline( v = x$TIB.interpolated, h = c( x$alpha / 2, 1 - x$alpha / 2 ), col = "blue" )

    if( !is.null( a$pch ) ) a.pch <- a$pch
    else a.pch <- "o"

    if( !is.null( a$col ) ) a.col <- a$col
    else a.col <- 1

    if( !is.null( a$cex ) ) a.cex <- a$cex
    else a.cex <- 1

    if( leg && type == "pvalue" ) {

	if( a.col != "gray" ) b.col <- "gray"
	else b.col <- "black"

	points( theta1, y2, pch = a.pch, col = b.col, cex = a.cex )

	abline( v = x$STIB.interpolated, h = c( x$alpha / 2, 1 - x$alpha / 2 ),
		col = "blue", lty = 2 )

	legend( "topright", legend = c( "TIB p-values (CI in vertical blue solid line)",
		"STIB p-values (CI in vertical blue dashed line)" ), col = c( a.col, b.col ),
		cex = a.cex, pch = a.pch, bty = "n" )

    } # end of if add STIB p-values and legend to plot stmts.

} # end of 'plot.tibbed' function.

print.tibRMed <- function( x, ... ) {

    cat( "\nCall:\n" )
    print( x$call )

    cat( "\nAchieved estimated p-value for lower = ", x$lower.p.value, "\n" )
    cat( "Associated nusiance parameter = ", x$lower.nuisance.par, "\n\n" )

    cat( "Achieved estimated p-value for upper = ", x$upper.p.value, "\n" )
    cat( "Associated nusiance parameter = ", x$upper.nuisance.par, "\n\n" )

    if( x$type == "TIB" ) cat( "Test-inversion Bootstrap ", ( 1 - x$alpha ) * 100, "% CI (Robbins-Monro Algorithm)\n" )
    else if( x$type == "STIB" ) cat( "(Studentized) Test-inversion Bootstrap ", ( 1 - x$alpha ) * 100, "% CI (Robbins-Monro Algorithm)\n" )
    print( x$result )

    invisible()

} # end of 'print.tibRMed' function.

plot.tibRMed <- function( x, ..., xlab = "Estimate", ylab = "p-value" ) {

    if( is.null( x$lower.iters ) ) stop( "plot: Must use keep.iters = TRUE to make plot." )

    lx <- x$lower.iters[, 1 ]
    ly <- x$lower.iters[, 3 ]

    ux <- x$upper.iters[, 1 ]
    uy <- x$upper.iters[, 3 ]

    x1 <- c( lx, ux )
    x2 <- c( ly, uy )

    plot( x1, x2, xlab = xlab, ylab = ylab, ..., type = "n" )
    points( lx, ly, pch = "l", col = "darkblue" )
    points( ux, uy, pch = "u", col = "lightblue" )

    abline( v = x$result, h = c( x$alpha / 2, 1 - x$alpha / 2 ), col = "blue" )

    invisible()

} # end of 'plot.tibRMed' function.



