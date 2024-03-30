# FNGAUSS.COM -- Common block to hold gaussian object variables

common	/ fngauss / radius, ar, pa, beta, scale, nx, ny, xgrid, ygrid

real	radius		# object radius
real	ar		# object aspect ratio
real	pa		# object position angle
real	beta		# exponent of ,offat profile
real	scale		# computed scaling factor
int	nx		# number of pixels along x dimension
int	ny		# number of pixels along y dimension
pointer	xgrid		# x pixel grid points
pointer	ygrid		# y pixel grid points

