is.even <- function( x ) return( ifelse(  x %% 2 == 0, TRUE, FALSE ) )
is.odd <- function( x ) return( ifelse( x %% 2 == 1, TRUE, FALSE ) )

even <- function( x ) return( x[ is.even( x ) ] )
odd <- function( x ) return( x[ is.odd( x ) ] )
