# FNEXP.COM -- Common block to hold exponential object variables

common	/ fnexp / radius, ar, pa, scale, nx, ny, xgrid, ygrid

real	radius		# object radius
real	ar		# object aspect ratio
real	pa		# object position angle
real	scale		# computed scaling factor
int	nx		# number of pixels along x dimension
int	ny		# number of pixels along y dimension
pointer	xgrid		# x pixel grid points
pointer	ygrid		# y pixel grid points

