define	MARGIN		0.05

# XYLIMIT -- Derive y plot limits for function dependent on x limits

# B.Simon	22-Feb-93	original

procedure xylimit (xplot, yplot, nplot, xmin, xmax, ymin, ymax)

real	xplot[ARB]	# i: Independent variable, sorted
real	yplot[ARB]	# i: Dependent variable
int	nplot		# i: Number of data points
real	xmin		# u: Minimum x value
real	xmax		# u: Maximum x value
real	ymin		# u: Minimum y value
real	ymax		# u: Maximum y value
#--
bool	first
int	i
real	ylo, yhi, ydif

begin
	if (nplot <= 0)
	    return

	# Find the smallest and largest y values that 
	# lie in the domain of the x values

	first = true
	do i = 1, nplot {
	    if (IS_INDEFR (yplot[i])) 
		next

	    if (xplot[i] < xmin || xplot[i] > xmax)
		next

	    if (first) {
		first = false
		ylo = yplot[i]
		yhi = yplot[i]
	    } else {
		ylo = min (ylo, yplot[i])
		yhi = max (yhi, yplot[i])
	    }
	}

	# Reset the minimum and maximum y values as necessary

	ydif = yhi - ylo

	if (IS_INDEFR (ymin)) {
	    ymin = ylo - MARGIN * ydif
	} else {
	    ymin = min (ymin, ylo - MARGIN * ydif)
	}

	if (IS_INDEFR (ymax)) {
	    ymax = yhi + MARGIN * ydif
	} else {
	    ymax = max (ymax, yhi + MARGIN * ydif)
	}

end
