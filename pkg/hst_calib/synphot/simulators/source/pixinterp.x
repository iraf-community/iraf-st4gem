include <imhdr.h>

#* HISTORY *
#* B.Simon	17-Feb-95	original

# PIXINTERP -- Interpolate image on a pixel grid

procedure pixinterp (im, xpcen, ypcen, grid, npix, array)

pointer	im		# i: psf image descriptor
real	xpcen		# i: psf center
real	ypcen		# i: psf center
real	grid[ARB]	# i: interpolant grid
int	npix		# i: size of interpolant array
real	array[npix,npix]# o: interpolated values
#--
int	ixlo, ixhi, iylo, iyhi, nx, ny, ix, iy
pointer	buf
real	xlo, xhi, ylo, yhi, xcen, ycen, sum, total

real    bilinval() # Compute a function value through bilinear interpolation
extern	bilinval
pointer	imgs2r()

begin
	# Read buffer containing image object

	ixlo = max (int (xpcen + grid[1]), 1)
	ixhi = min (int (xpcen + grid[npix+1] + 1.0), IM_LEN(im,1))
	iylo = max (int (ypcen + grid[1]), 1)
	iyhi = min (int (ypcen + grid[npix+1] + 1.0), IM_LEN(im,2))
	nx = (ixhi - ixlo) + 1
	ny = (iyhi - iylo) + 1

	buf = imgs2r (im, ixlo, ixhi, iylo, iyhi)

	# Initialize bilinear interpolator

	call bilinset (buf, nx, ny)

	# Integrate interpolant over pixel boundary for interpolated value
	# (Grid is defined to have npix+1 values)

	total = 0.0
	xcen = xpcen - (ixlo - 1)
	ycen = ypcen - (iylo - 1)

	do iy = 1, npix {
	    ylo = grid[iy] + ycen
	    yhi = grid[iy+1] + ycen
	    sum = INDEFR

	    do ix = 1, npix {
		xlo = grid[ix] + xcen
		xhi = grid[ix+1] + xcen

		call simpson2 (xlo, xhi, ylo, yhi, bilinval, sum, array[ix,iy])
		total = total + array[ix,iy]
	    }
	}

	# Noramlize interpolated array to one

	if (total > 0.0) {
	    total = 1.0 / total
	    call amulkr (array, total, array, npix*npix)
	}

end
