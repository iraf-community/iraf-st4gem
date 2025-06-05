include <iraf77.h>
include	<gset.h>

# UGERBS -- Plot symmetrical error bars in X or Y with different
# magnitude for each data point. 

procedure ugerbs (gp, x, y, e, npts, axis, istat)

pointer	gp			# Graphics descriptor
real	x[ARB], y[ARB]		# Data values
real	e[ARB]			# Errors (half amplitude)
int	npts			# Number of values
int	axis			# X or Y axis
int	istat

int	pix
pointer	sp, size

begin
	istat = ER_OK
	# Allocate a buffer for the array of error bar sizes
	call smark  (sp)
	call salloc (size, npts, TY_REAL)

	# Find the size of the markers
	call amulkr (e, 2.0, Memr[size], npts)

	# Specify the marker size in WCS
	call anegr  (Memr[size], Memr[size], npts)

	if (axis == 1)  {		# X
	    do pix = 1, npts
		# Draw the horizontal bar
		iferr {
		    call gmark (gp, x[pix], y[pix], 
			GM_HEBAR, Memr[size+pix-1], 1.0)
		} then
		    istat = ER_GRAPHERRBAR

	} else if (axis == 2)  {	# Y
	    do pix = 1, npts 
		# Draw the vertical bar
		iferr {
		    call gmark (gp, x[pix], y[pix], 
			GM_VEBAR, 1.0, Memr[size+pix-1])
		} then
		    istat = ER_GRAPHERRBAR
	}

	call sfree (sp)
end
