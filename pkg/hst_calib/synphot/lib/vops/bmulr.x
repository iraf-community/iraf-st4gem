# BMULR -- Multiply two vectors, excluding INDEFRs

procedure bmulr (a, b, c, npix)

real	a[ARB], b[ARB], c[ARB]
int	npix, i

begin
	do i = 1, npix
	    if ( !IS_INDEFR (a[i]) && !IS_INDEFR (b[i]) )
	        c[i] = a[i] * b[i]
	    else
	        c[i] = INDEFR
end
