include "otf.h"
define	none_		90

#* HISTORY *
#* B.Simon	11-Aug-95	original

# DSFIMAGE -- Create a structure  holding the point spread functions

pointer procedure dsfimage (dsfcat, obsmode, apscale, dynrange)

char	dsfcat[ARB]	# i: catalog of psf file names
char	obsmode[ARB]	# i: observation mode
double	apscale		# i: aperture scale (degrees)
real	dynrange	# i: dynamic range of psf
#--
int	npix, jrow
pointer	sp, dsf, im, tp, cp, filename, grid, buffer
real	xpcen, ypcen, range, psfscale, scale

string	obscol   "OBSMODE"
string	filecol  "FILENAME"
string	badsize  "DSF has invalid scale; detector scale assumed"
string	toobig   "DSF size too large"

int	access(), imaccess()
pointer	immap(), tbtopn()

begin
	# Allocate memory for temporary strings

	call smark (sp)
	call salloc (filename, SZ_FNAME, TY_CHAR)

	# Open dsf catalog

	call lastfile (dsfcat, Memc[filename], SZ_FNAME)
	if (access (Memc[filename], 0, 0) == NO)
	    goto none_

	# If the "catalog file" is a single image, open that
	# image as the sole dsf. Otherwise, open the catalog and
	# search for the matching dsf

	if (imaccess (Memc[filename], READ_ONLY) == YES &&
	    imaccess (Memc[filename], NEW_FILE) == YES    ) {

	    im = immap (Memc[filename], READ_ONLY, 0)

	} else {
	    tp = tbtopn (Memc[filename], READ_ONLY, 0)
	    call findmode (tp, obscol, obsmode, jrow)

	    if (jrow == 0)
		goto none_
		
	    call syncolptr (tp, filecol, 2, cp)
	    call tbegtt (tp, cp, jrow, Memc[filename], SZ_FNAME)
	    im = immap (Memc[filename], READ_ONLY, 0)

	    call tbtclo (tp)
	}


	# Compute center of dsf and pixel scale

	call pixregion (im, dynrange, xpcen, ypcen, npix) 
	call pixsize (im, psfscale)
	if (psfscale <= 0.0) {
	    call synphotwarn (badsize, Memc[filename])
	    psfscale = apscale
	}

	range = npix * psfscale

	# Compute dsf array size and allocate arrays based on this size

	npix = int (range / apscale) + 1
	if (npix > MAXOTF)	# sanity check on dsf size
	    call printerr_int (toobig, npix)

	call salloc (grid, npix+1, TY_REAL)
	call malloc (buffer, npix*npix, TY_REAL)

	# Interpolate dsf on detector pixel grid

	scale = apscale / psfscale
	call setgrid (scale, Memr[grid], npix+1)

	call pixinterp (im, xpcen, ypcen, Memr[grid], npix, Memr[buffer])
	call imunmap (im)

	# Allocate and fill psf structure

	call malloc (dsf, LEN_OTFSTRUCT, TY_STRUCT)

	OTF_NUMBER(dsf) = 1
	OTF_NXPIX(dsf) = npix
	OTF_NYPIX(dsf) = npix
	OTF_WAVPTR(dsf) = 5500.
	OTF_BUFFER(dsf) = buffer

	# Return point spread function descriptor

	call sfree (sp)
	return (dsf)

	# Alternate return if point spread function not found

none_	call sfree (sp)
	return (NULL)

end
