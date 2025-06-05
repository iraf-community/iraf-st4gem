include <iraf77.h>

# UGMARK -- Polymarker.  Plot a sequence of npts markers at the positions
# given by successive WCS coordinates pairs (X[i],Y[i]).
# All the markers will be of the same type and size.

procedure ugmark (gp, x, y, npts, marktype, xsize, ysize, istat)

pointer	gp			# graphics descriptor
real	x[ARB], y[ARB]		# points defining the polyline
int	npts
int	marktype		# Type of marker (like in GIO)
real	xsize			# size of marker (like in GIO )
real	ysize			# size of marker (like in GIO )
int	istat

begin
	istat = ER_OK

	# We need at least one point.
	if (npts < 1)  {
	   istat = ER_GRAPHBADNPTS
	   return
	}

	# Mark the plot
	iferr (call gpmark (gp, x, y, npts, marktype, xsize, ysize))
	   istat = ER_GRAPHMARKER
end
