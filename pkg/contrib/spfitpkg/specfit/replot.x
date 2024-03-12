###########################################################################
#                    Center for Astrophysical Sciences
#                        Johns Hopkins University

#  Synopsis:	procedure replot (gfd, gt, x, pix, npts)
#		pointer	gfd
#		pointer	gt
#		real	x[ARB]
#		real	pix[ARB]
#		int	npts

#  Synopsis:	procedure overplot (gfd, gt, x, pix, npts)
#		pointer	gfd
#		pointer	gt
#		real	x[ARB]
#		real	pix[ARB]
#		int	npts

#  Description:	REPLOT -- Replot the current array
#  Description:	OVERPLOT -- Overplot the specified array on the cuurent plot

#  Arguments:	pointer	gfd, gt		Pointers for graphics control
#		real	x[ARB]		List of x positions
#		real	pix[ARB]	The data array to plot
#		int	npts		Number of points in pix

#  Returns:	None

#  Notes:	Adapted from splot in the onedspec IRAF package

#  History:	June	1987	Gerard Kriss

###########################################################################

procedure replot (gfd, gt, x, pix, npts)
pointer	gfd
pointer	gt
real	x[ARB]
real	pix[ARB]
int	npts

begin
	call gclear (gfd)
	call gswind (gfd, x[1], x[npts], INDEF, INDEF)
	call gascale (gfd, pix, npts, 2)
	call gt_swind (gfd, gt)
	call gt_labax (gfd, gt)
	call gpline (gfd, x, pix, npts)
end

procedure overplot (gfd, gt, x, pix, npts)
pointer	gfd
pointer	gt
real	x[ARB]
real	pix[ARB]
int	npts

begin
	call gpline (gfd, x, pix, npts)
end







