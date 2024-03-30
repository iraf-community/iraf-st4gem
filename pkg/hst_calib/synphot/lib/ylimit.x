define	MARGIN		0.05

# YLIMIT -- Determine plot y limits from array of ordinates

# B.Simon	23-Feb-93	original

procedure ylimit (yplot, nplot, ymin, ymax)

real	yplot[ARB]	# i: Dependent variable
int	nplot		# i: Number of data points
real	ymin		# u: Minimum y value
real	ymax		# u: Maximum y value
#--
int	i
real	ylo, yhi, ydif

begin
	if (nplot <= 0)
	    return

	do i = 1, nplot {
	    if (IS_INDEFR (yplot[i])) 
		next

	    ylo = min (ylo, yplot[i])
	    yhi = max (yhi, yplot[i])
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
