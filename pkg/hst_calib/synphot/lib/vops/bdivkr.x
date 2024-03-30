# BDIVK -- Divide a vector by a constant (generic).  No divide by zero checking
# is performed.  Check for INDEFs

procedure bdivkr (a, b, c, npix)

real	a[ARB]
real	b
real	c[ARB]
int	npix, i

begin
	do i = 1, npix
	    if ( !IS_INDEFR (a[i]) )
	        c[i] = a[i] / b
	    else
	        c[i] = INDEFR
end
