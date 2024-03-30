# HGLINE -- Draw a stepped curve of the histogram data.

procedure hgline (gp, xdata, ydata, npts )

pointer	gp		# Graphics descriptor
real	xdata[ARB]	# X coordinates data
real	ydata[ARB]	# Y coordinates of data
int	npts		# Number of datapoints

int	pixel
real	x, y

begin

	# Do the first horizontal line
	x = xdata[1]
	y = ydata[1]
	call gamove (gp, x, y)
	x = xdata[2]
	call gadraw( gp, x, y )

	do pixel = 2, npts-1 {
	    y = ydata[pixel]
	    # vertical connection
	    call gadraw (gp, x, y)
	    # horizontal line
	    x = xdata[pixel+1]
	    call gadraw (gp, x, y)
	}
	# Draw last segments
	y = ydata[npts]
	call gadraw (gp, x, y)
end
