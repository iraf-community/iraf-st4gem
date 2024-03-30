# FNLIST.COM -- Common block to hold list profile object variables

common	/ fnlist / radius, ar, pa, coef, scale, prof, nx, ny, xgrid, ygrid

real	radius		# object radius
real	ar		# object aspect ratio
real	pa		# object position angle
real	coef		# profile size in pixels
real	scale		# computed scaling factor
pointer	prof		# interpolated profile
int	nx		# number of pixels along x dimension
int	ny		# number of pixels along y dimension
pointer	xgrid		# x pixel grid points
pointer	ygrid		# y pixel grid points

