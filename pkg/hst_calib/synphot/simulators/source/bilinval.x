#* HISTORY *
#* B.Simon	17-Feb-95	original

# BILINVAL -- Compute a function value through bilinear interpolation

real procedure bilinval (x, y)

real	x		# i: x coordinate
real	y		# i: y coordinate
#--
include	"bilin.com"

int	ix, iy, p1, p2, p3, p4
real	ax, bx, ay, by, z

begin
	# Locate grid points which bound the interpolated point

	ix = min (nx-1, max (1, int(x)))
	iy = min (ny-1, max (1, int(y)))

	p1 = ny * (iy - 1) + (ix - 1)
	p3 = p1 + ny
	p2 = p1 + 1
	p4 = p3 + 1

	# Compute interpolated value

	bx = x - ix
	ax = 1.0 - bx

	by = y - iy
	ay = 1.0 - by

	z = ax * ay * Memr[buf+p1] + bx * ay * Memr[buf+p2] +
	    ax * by * Memr[buf+p3] + bx * by * Memr[buf+p4]

	return (max (z, 0.0))
end
