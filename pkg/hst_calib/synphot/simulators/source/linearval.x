#* HISTORY *
#* B.Simon	11-Jul-95	derived from bilinval

# LINEARVAL -- Compute a function value through linear interpolation

real procedure linearval (x)

real	x		# i: x coordinate
#--
include	"linear.com"

int	ix
real	ax, bx, z

begin
	# Locate grid points which bound the interpolated point

	ix = min (nx-1, max (1, int(x)))

	# Compute interpolated value

	bx = x - ix
	ax = 1.0 - bx

	z = ax * Memr[buf+ix-1] + bx * Memr[buf+ix]
	return (max (z, 0.0))

end
