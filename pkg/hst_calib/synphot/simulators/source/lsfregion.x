include <imio.h>
include <imhdr.h>
define 	NC		5

#* HISTORY *
#* B.Simon	10-Jul-95	derived from pixregion and pixsize
#* B.Simon	29-Sep-95	compute center to fractional pixel
#* B.Simon	10-Oct-95	subtract background

# LSFREGION -- Find extent and region of line spread function

procedure lsfregion (im, dynrange, npix, pcenter, lsfscale)

pointer	im		# i: image descriptor
real	dynrange	# i: dynamic range
int	npix		# o: length of a side of the psf
real	pcenter		# o: xcenter of psf
real	lsfscale	# o: size of lsf pixel
#--
int	nstep, step
pointer	sp, type, buf, ptr, mw
real	start, mean, sigma, max, total, range, r[1], w[1], cd[1,1]

string	badaxis    "|pixels|"
string	baddimen   "LSF image is not one dimensional"

int	aravr(), checkdim(), strdic()
pointer	imgl1r(), mw_openim()
real	imgetr()

begin
	# Allocate memory for temporary strings

	call smark (sp)
	call salloc (type, SZ_FNAME, TY_CHAR)

	# Check dimensions of psf

	if (checkdim (im) != 1)
	    call printerr_str (baddimen, IM_NAME(im))

	# Read lsf center location from image header

	iferr {
	    pcenter = imgetr (im, "CENTER")

	} then {
	    pcenter = 0.0

	} else {
	    if (pcenter < 1.0 || pcenter > real(IM_LEN(im,1)))
		pcenter = 0.0
	}

	# If center not found or invalid, compute by centroiding

	if (pcenter == 0.0) {
	    start = 0.5 * (IM_LEN(im,1) + 1)
	    call centroid1 (im, start, NC, pcenter)
	}

	# Find pixels containing most of flux

	nstep = IM_LEN(im,1) / 2
	range = 1.0 / dynrange

	buf = imgl1r (im)
	npix = aravr (Memr[buf], IM_LEN(im,1), mean, sigma, 0.0)

	ptr = int (pcenter) + buf - 1
	max = Memr[ptr] - mean

	do step = 1, nstep {
	    total = 0.5 * (Memr[ptr+step] + Memr[ptr-step]) - mean

	    if (total / max < range) {
		npix = 2 * step + 1
		break
	    }
	}

	# Calculate pixel size from cdmatrix

	mw = mw_openim (im)
	call mw_gwtermr (mw, r, w, cd, 1)

	lsfscale = cd[1,1]

	call mw_close (mw)

	# Check for bad axis types

	call imgstr (im, "CTYPE1", Memc[type], SZ_FNAME)
	call strfix (Memc[type])

	if (strdic (Memc[type], Memc[type], SZ_FNAME, badaxis) != 0)
	    lsfscale = 0.0


end

