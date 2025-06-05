include <iraf77.h>

# UGPVPT  -- Set viewport

procedure ugpvpt (gp, x1, x2, y1, y2, istat)

pointer	gp
real	x1, x2		# range of world coordinates in X
real	y1, y2		# range of world coordinates in Y
int	istat

begin
	istat = ER_OK

	# Set viewport
	iferr (call gsview (gp, x1, x2, y1, y2))
	   istat = ER_GRAPHSETVP
end
