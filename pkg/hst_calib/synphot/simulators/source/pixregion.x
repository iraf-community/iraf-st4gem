include <imio.h>
include <imhdr.h>

define	NC		5	# Box size of region to search for center

#* HISTORY *
#* B.Simon	17-Feb-95	original
#* B.Simon	28-Sep-95	compute center to fractional pixel
#* B.Simon	10-Oct-95	subtract background

# PIXREGION -- Find region of image containing psf

procedure pixregion (im, dynrange, xpcen, ypcen, npix)

pointer	im		# i: image descriptor
real	dynrange	# i: dynamic range
real	xpcen		# o: x center of psf
real	ypcen		# o: y center of psf
int	npix		# o: length of a side of the psf
#--
int	linenum, ixc, iyc, xstep, ystep, nstep, step
int	xlo, xhi, ylo, yhi, nold, nbuf
pointer	buf
real	xstart, ystart, mean, sigma, max, oldtot, total, avgtot, range

string	baddimen   "PSF image is not two dimensional"

int	aravr(), checkdim()
pointer	imgl2r(), imgs2r()
real	imgetr(), asumr()

begin
	# Check dimensions of psf

	if (checkdim (im) != 2)
	    call printerr_str (baddimen, IM_NAME(im))

	# Read psf center location from image header

	iferr {
	    xpcen = imgetr (im, "XCENTER")
	    ypcen = imgetr (im, "YCENTER")

	} then {
	    xpcen = 0.0
	    ypcen = 0.0

	} else {
	    if (xpcen < 1.0 || xpcen > real(IM_LEN(im,1)))
		xpcen = 0.0

	    if (ypcen < 1 || ypcen > real(IM_LEN(im,2)))
		ypcen = 0.0
	}

	# If center not found or invalid, compute by centroiding

	if (xpcen == 0.0 || ypcen == 0.0) {
	    xstart = 0.5 * (IM_LEN(im,1) + 1)
	    ystart = 0.5 * (IM_LEN(im,2) + 1)

	    call centroid2 (im, xstart, ystart, NC, xpcen, ypcen)
	}

	# Estimate background from a line in the image

	linenum = 0.25 * IM_LEN(im,2)
	buf = imgl2r (im, linenum)
	npix = aravr (Memr[buf], IM_LEN(im,1), mean, sigma, 0.0)

	# Find rectangle containing most of flux

	ixc = xpcen
	iyc = ypcen
	xstep = min (ixc - 1, IM_LEN(im,1) - ixc)
	ystep = min (iyc - 1, IM_LEN(im,2) - iyc)
	nstep = (min (xstep, ystep))
	range = 1.0 / dynrange

	do step = 1, nstep {
	    xlo = ixc - step
	    xhi = ixc + step
	    ylo = iyc - step
	    yhi = iyc + step

	    if (xlo < 1 || xhi > IM_LEN(im,1) || 
		ylo < 1 || yhi > IM_LEN(im,2)   ) {
		npix = npix - 2
		break
	    }

	    buf = imgs2r (im, xlo, xhi, ylo, yhi)
	    if (step == 1) {
		max = Memr[buf+4] - mean
		oldtot = max
		nold = 1
	    }

	    npix = 2 * step + 1
	    nbuf = npix ** 2

	    total = asumr (Memr[buf], nbuf) - nbuf * mean
	    avgtot = (total - oldtot) / (nbuf - nold)

	    if ((avgtot / max) < range) {
		npix = npix - 2
		break
	    }

	    oldtot = total
	    nold = nbuf
	}

end
