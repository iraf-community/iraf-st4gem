# FNIMAGE.COM -- Common block to hold image object variables

common	/ fnimage / im, radius, npix, grid, nx, ny, buf, xcen, ycen 

pointer	im		# image descriptor
real	radius		# object radius
int	npix		# number of pixels in grid
pointer	grid		# pixel grid points
int	nx		# x length of image buffer
int	ny		# y length of image buffer
pointer	buf		# image buffer
real	xcen		# x coord of central pixel of image
real	ycen		# y coord of central pixel of image
