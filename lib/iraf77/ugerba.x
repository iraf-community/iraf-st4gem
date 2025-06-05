include <iraf77.h>
include	<gset.h>

# UGERBA -- Plot asymmetrical error bars in X or Y with different
# magnitude for each data point.

procedure ugerba (gp, x, y, ne, pe, npts, axis, istat)

pointer	gp			# Graphics descriptor
real	x[ARB], y[ARB]		# Data values
real	ne[ARB], pe[ARB]	# Low and high errors 
int	npts			# Number of values
int	axis			# X or Y axis?
int	istat

int	pix
pointer	sp, loerr, hierr, size, cent

begin
	istat = ER_OK
	# Allocate buffers for the array of error bar sizes and 
	# corrected positions
	call smark  (sp)
	call salloc (loerr, npts, TY_REAL)
	call salloc (hierr, npts, TY_REAL)
	call salloc (size, npts, TY_REAL)
	call salloc (cent, npts, TY_REAL)

	# Make sure the lower and upper errors are what we think they are
	call aabsr (ne, loerr, npts)
	call aabsr (pe, hierr, npts)

	# Find the size of an equivalent symmetrical error bar
	call aaddr (loerr, hierr, Memr[size], npts)

	# Specify the marker size in WCS
	call anegr (Memr[size], Memr[size], npts)

	# Find the equivalent location of a symmetrical error bar
	call asubr (Memr[hierr], Memr[loerr], Memr[cent], npts)
	call adivkr (Memr[cent], 2.0, Memr[cent], npts)

	if (axis == 1) {		# X 
	    call asubr (x, Memr[cent], Memr[cent], npts)

	    do pix = 1, npts
		# Draw the horizontal bar
		iferr (call gmark (gp, Memr[cent+pix-1], y[pix], 
		    GM_HEBAR, Memr[size+pix-1], 1.0))
	 	   istat = ER_GRAPHERRBAR
	} else if (axis == 2) {		# Y 
	    call asubr (y, Memr[cent], Memr[cent], npts)

	    do pix = 1, npts
		# Draw the vertical bar
		iferr (call gmark (gp, x[pix], Memr[cent+pix-1], 
		    GM_VEBAR, 1.0, Memr[size+pix-1]))
	 	   istat = ER_GRAPHERRBAR
	}

	call sfree (sp)
end
