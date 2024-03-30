include <imhdr.h>

#* HISTORY *
#* B.Simon	17-Feb-95	original

# LSFINTERP -- Interpolate lsf image on a pixel grid

procedure lsfinterp (im, pcenter, grid, npix, array)

pointer	im		# i: psf image descriptor
real	pcenter		# i: psf center
real	grid[ARB]	# i: interpolant grid
int	npix		# i: size of interpolant array
real	array[ARB]	# o: interpolated values
#--
int	ixlo, ixhi, nx, ix
pointer	buf
real	xlo, xhi, xcen, term, total

real    linearval() # Compute a function value through linear interpolation
extern	linearval
pointer	imgs1r()

begin
	# Read buffer containing lsf image

	ixlo = max (int (pcenter + grid[1]), 1)
	ixhi = min (int (pcenter + grid[npix+1] + 1.0), IM_LEN(im,1))
	nx = (ixhi - ixlo) + 1

	buf = imgs1r (im, ixlo, ixhi)

	# Initialize bilinear interpolator

	call linearset (buf, nx)

	# Integrate interpolant over pixel boundary for interpolated value
	# (Grid is defined to have npix+1 values)

	total = 0.0
	term = INDEFR
	xcen = pcenter - (ixlo - 1)

	do ix = 1, npix {
	    xlo = grid[ix] + xcen
	    xhi = grid[ix+1] + xcen

	    call simpson1 (xlo, xhi, linearval, term, array[ix])
	    total = total + array[ix]
	}

	# Noramlize interpolated array to one

	if (total > 0.0) {
	    total = 1.0 / total
	    call amulkr (array, total, array, npix)
	}

end
