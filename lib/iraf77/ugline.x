include <iraf77.h>

# UGLINE -- Polyline.  Draw a line connecting the points (X[i],Y[i]), i.e.,
# move to the first point and draw a straight line from there to the second
# point, from the second to the third, and so on.

procedure ugline (grdscr, x, y, npts, istat)

pointer	grdscr			# graphics descriptor
real	x[ARB], y[ARB]		# points defining the polyline
int	npts
int	istat

begin
	istat = ER_OK

	# We need at least two points.
	if (npts < 2)  {
	   istat = ER_GRAPHBADNPTS
	   return
	}

	# draw the line
	iferr (call gpline (grdscr, x, y, npts))
	   istat = ER_GRAPHLINE
end
