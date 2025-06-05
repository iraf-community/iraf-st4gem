include <iraf77.h>

#  UGHIST -- Draw a histogram  (stepped curve or bar graph) through the
#  points. 

procedure ughist (gp, xdata, ydata, npts, hstfmt, istat)

pointer	gp		# Graphics descriptor
real	xdata[ARB]	# X coordinates of the line endpoints
real	ydata[ARB]	# Y coordinates of the line endpoints
int	npts		# Number of line endpoints
int	hstfmt		# Histogram type 
int	istat		# Status return

begin
	istat = ER_OK

	switch (hstfmt)  {
	case GR_STEP:
	iferr (call hgline (gp, xdata, ydata, npts))
	   istat = ER_GRAPHISTGM
	case GR_BAR:
	iferr (call fhglin (gp, xdata, ydata, npts))
	   istat = ER_GRAPHISTGM
	default:
	   istat = ER_GRAPHISTGM
	}
end
