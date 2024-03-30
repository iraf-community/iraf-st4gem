# BMULKR -- Multiply a constant times a vector avoiding INDEFRs

procedure bmulkr (a, b, c, npix)

real	a[ARB]
real	b
real	c[ARB]
int	npix, i

begin
	do i = 1, npix {
	    if ( !IS_INDEFR (a[i]) )
	       c[i] = a[i] * b
	    else
	       c[i] = a[i]
	}
end
