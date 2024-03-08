# xyz_limits -- find the limits in X and Y
# The values x1, x2, y1, y2 were gotten from the cl.  If x1 was not
# specified (i.e. is INDEF), the minimum X value in the data array will be
# assigned to x1.  If x2 is INDEF it will be assigned the maximum value in
# the X array.  Similarly for y1 and y2.
# Y values are not changed if ndim is one.
#
# Phil Hodge, 15-Dec-1993  Subroutine created.

procedure xyz_limits (x, y, nrows, ndim, x1, x2, y1, y2)

double	x[ARB]		# i: X values
double	y[ARB]		# i: Y values
int	nrows		# i: size of x and y arrays
int	ndim		# i: dimension of image (1 or 2)
double	x1, x2		# io: limits on X values
double	y1, y2		# io: limits on Y values
#--
double	xmin, xmax	# mininum and maximum X values from table
double	ymin, ymax	# mininum and maximum Y values from table
int	i
bool	findx, findy	# true if we need to find either x1 or x2 (y1 or y2)

begin
	# Find out whether we need to examine the data for either
	# the minimum or maximum values.  If we need to find either,
	# we'll find both but only update the unspecified one(s).
	findx = (IS_INDEFD(x1) || IS_INDEFD(x2))
	findy = (IS_INDEFD(y1) || IS_INDEFD(y2))

	if (ndim == 1)
	    findy = false

	if (findx) {

	    xmin = x[1]			# initial values
	    xmax = x[1]

	    do i = 1, nrows {
		xmin = min (xmin, x[i])
		xmax = max (xmax, x[i])
	    }
	    if (IS_INDEFD(x1))
		x1 = xmin
	    if (IS_INDEFD(x2))
		x2 = xmax
	}

	if (findy) {

	    ymin = y[1]			# initial values
	    ymax = y[1]

	    do i = 1, nrows {
		ymin = min (ymin, y[i])
		ymax = max (ymax, y[i])
	    }
	    if (IS_INDEFD(y1))
		y1 = ymin
	    if (IS_INDEFD(y2))
		y2 = ymax
	}
end
