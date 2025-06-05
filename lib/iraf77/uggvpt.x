include <iraf77.h>

# UGGVPT  -- Get viewport

procedure uggvpt (gp, x1, x2, y1, y2, istat)

pointer	gp
real	x1, x2		# range of world coordinates in X
real	y1, y2		# range of world coordinates in Y
int	istat

begin
	istat = ER_OK

	# Get viewport
	iferr (call ggview (gp, x1, x2, y1, y2))
	   istat = ER_GRAPHGETVP
end
