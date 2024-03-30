# BSUB -- Subtract two real vectors

procedure bsubr (a, b, c, npix)

real	a[ARB], b[ARB], c[ARB]
int	npix, i

begin
	do i = 1, npix {
	    if ( !IS_INDEFR (a[i]) && !IS_INDEFR (b[i]) )
	       c[i] = a[i] - b[i]
	    else
	       c[i] = INDEFR
	}
end
