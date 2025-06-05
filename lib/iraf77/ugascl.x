include <iraf77.h>

# UGASCL -- Scale the data in X, in Y or in both coordinates.

procedure ugascl (gp, x, y, npts, axis, istat)

pointer	gp			# graphics descriptor
real	x[ARB], y[ARB]		# points defining the data
int	npts			# number of data points
int	axis			# 1 = X axis ; 2 = Y axis ; 3 = both axes
int	istat

begin
	istat = ER_OK

	# We need at least two points.
	if (npts < 2)  {
	   istat = ER_GRAPHBADNPTS
	   return
	}

	if (axis < 1 || axis > 3)  {
	   istat = ER_GRAPHBADAXES
	   return
	}

	# What do we have to scale ?
	if (axis == 1 || axis == 3)  
	    iferr (call gascale (gp, x, npts, 1)) {	# Scale in X
		istat = ER_GRAPHSCALE
		return
	    }

	if (axis == 2 || axis == 3)
	    iferr (call gascale (gp, y, npts, 2)) {	# Scale in Y
		istat = ER_GRAPHSCALE
		return
	    }
end
