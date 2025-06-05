include <iraf77.h>

# UGPWND  -- Set data window

procedure ugpwnd (gp, x1, x2, y1, y2, istat)

pointer	gp
real	x1, x2		# range of world coordinates in X
real	y1, y2		# range of world coordinates in Y
int	istat

begin
	istat = ER_OK

	# Set data window
	iferr (call gswind (gp, x1, x2, y1, y2))
	   istat = ER_GRAPHSETWND
end
