include	<gset.h>

# PLTBAR -- Plot an array of horizontal or vertical bars

procedure pltbar (gp, x, y, npts, marktype, size)

pointer	gp		# i: graphics descriptor
real	x[ARB]		# i: x coordinate of bar center
real	y[ARB]		# i: y coordinate of bar center
int	npts		# i: number of bars to plot
int	marktype	# i: bar type
real	size[ARB]	# i: bar size
#--
int	i

begin
	switch (marktype) {
	case GM_HLINE:
	    do i = 1, npts
		call gmark (gp, x[i], y[i], marktype, size[i], 1.0)
	case GM_VLINE:
	    do i = 1, npts
		call gmark (gp, x[i], y[i], marktype, 1.0, size[i])
	default:
	    do i = 1, npts
		call gmark (gp, x[i], y[i], marktype, size[i], size[i])
	}
end


