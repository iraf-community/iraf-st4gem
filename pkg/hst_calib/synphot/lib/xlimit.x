define	MARGIN		0.05

# XLIMIT -- Derive x plot limits for function

# B.Simon	24-Feb-93	original

procedure xlimit (xplot, nplot, xmin, xmax)

real	xplot[ARB]	# i: Independent variable, sorted
int	nplot		# i: Number of data points
real	xmin		# u: Minimum x value
real	xmax		# u: Maximum x value
#--
real	xlo, xhi, xdif

begin
	if (nplot <= 0)
	    return

	# Find the smallest and largest x values
	# Array is assumed to be sorted

	xlo = min (xplot[1], xplot[nplot])
	xhi = max (xplot[1], xplot[nplot])
	xdif = xhi - xlo

	if (IS_INDEFR (xmin)) {
	    xmin = xlo - MARGIN * xdif
	} else {
	    xmin = min (xmin, xlo - MARGIN * xdif)
	}

	if (IS_INDEFR (xmax)) {
	    xmax = xhi + MARGIN * xdif
	} else {
	    xmax = max (xmax, xhi + MARGIN * xdif)
	}

end
