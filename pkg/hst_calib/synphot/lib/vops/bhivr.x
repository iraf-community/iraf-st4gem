# BHIVR -- Compute the high value (maximum) of a vector, avoiding INDEFRs

real procedure bhivr (a, npix)

real	a[ARB]
int	npix
real	high, pixval
int	i

begin
	high = -INDEFR

	do i = 1, npix {
	    pixval = a[i]
	    if (pixval > high && !IS_INDEFR (pixval))
	        high = pixval
	}

	if ( high == -INDEFR )
	    high = INDEFR

	return (high)
end
