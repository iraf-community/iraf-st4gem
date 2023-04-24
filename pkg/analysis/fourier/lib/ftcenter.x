# ft_center -- shift data in arrays
# The values in arrays REALPT and IMAGPT are shifted by SHIFT elements
# to the right if SHIFT > 0 and to the left if SHIFT < 0.  Elements
# shifted off one end fill in the spaces vacated at the other end
# (i.e. the shift is cyclic).
# FT_CRPIX is not updated.
#
# Phil Hodge, 20-Jul-1988  Subroutine created.

procedure ft_center (npts, shift, realpt, imagpt, dx)

int	npts		# i: length of arrays
int	shift		# i: amount of shift (+ or -)
real	realpt[npts]	# io: real part
real	imagpt[npts]	# io: imaginary part
real	dx[ARB]		# io: scratch space
#--
int	sh		# abs (shift)
int	n		# npts - abs(shift)

begin
	if (shift > 0) {

	    n = npts - shift

	    call amovr (realpt, dx, n)
	    call amovr (realpt[n+1], realpt[1], shift)
	    call amovr (dx, realpt[shift+1], n)

	    call amovr (imagpt, dx, n)
	    call amovr (imagpt[n+1], imagpt[1], shift)
	    call amovr (dx, imagpt[shift+1], n)

	} else {

	    sh = abs (shift)
	    n = npts - sh

	    call amovr (realpt[sh+1], dx, n)
	    call amovr (realpt[1], realpt[n+1], sh)
	    call amovr (dx, realpt[1], n)

	    call amovr (imagpt[sh+1], dx, n)
	    call amovr (imagpt[1], imagpt[n+1], sh)
	    call amovr (dx, imagpt[1], n)
	}
end
