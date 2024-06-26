# Copyright(c) 1986 Association of Universities for Research in Astronomy Inc.

# BLOV -- Compute the low value (minimum) of a vector.

real procedure blovr (a, npix)

real	a[ARB]
int	npix
real	low, pixval
int	i

begin
	low = a[1]

	do i = 1, npix {
	    pixval = a[i]
	    if (pixval < low && !IS_INDEFR (pixval))
		low = pixval
	}

	return (low)
end
