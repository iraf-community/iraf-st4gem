include <iraf77.h>

# UGGWND  -- Get window

procedure uggwnd (gp, x1, x2, y1, y2, istat)

pointer	gp
real	x1, x2		# range of world coordinates in X
real	y1, y2		# range of world coordinates in Y
int	istat

begin
	istat = ER_OK

	# Get window
	iferr ( call ggwind (gp, x1, x2, y1, y2))
	   istat = ER_GRAPHGETWND
end
