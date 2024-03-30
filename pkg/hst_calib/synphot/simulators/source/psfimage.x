include	<imhdr.h>
include "otf.h"

#* HISTORY *
#* B.Simon	17-Feb-95	original

# PSFIMAGE -- Create a structure  holding the point spread functions

pointer procedure psfimage (psfcat, obsmode, apscale, dynrange, nsub)

char	psfcat[ARB]	# i: catalog of psf file names
char	obsmode[ARB]	# i: observation mode
double	apscale		# i: aperture scale (degrees)
real	dynrange	# i: dynamic range of psf
int	nsub		# i: pixel subsampling
#--
int	numpsf, ipsf, npix
pointer	sp, xpcen, ypcen, psfscale, image
pointer	psf, im, wavepsf, grid, buffer, array
real	range, scale

string	badsize  "PSF has invalid scale; detector scale assumed"
string	toobig   "PSF size too large"

begin
	# Retrieve psf images that match this obsmode

	call otfopen (psfcat, obsmode, im, wavepsf, numpsf)

	# Allocate arrays based on number of psfs

	call smark (sp)
	call salloc (xpcen, numpsf, TY_REAL)
	call salloc (ypcen, numpsf, TY_REAL)
	call salloc (psfscale, numpsf, TY_REAL)
	call salloc (image, SZ_FNAME, TY_CHAR)

	# Compute largest extent of any psf

	range = 0.0
	do ipsf = 0, numpsf-1 {
	    call pixregion (Memi[im+ipsf], dynrange, Memr[xpcen+ipsf], 
			    Memr[ypcen+ipsf], npix)

	    call pixsize (Memi[im+ipsf], Memr[psfscale+ipsf])
	    if (Memr[psfscale+ipsf] <= 0.0) {
		call strcpy (IM_HDRFILE(Memi[im+ipsf]), Memc[image], SZ_FNAME)
		call synphotwarn (badsize, Memc[image])
		Memr[psfscale+ipsf] = apscale
	    }

	    range = max (npix * Memr[psfscale+ipsf], range)
	}

	# Compute psf array size and allocate arrays based on this size

	npix = int (range / apscale) + 1
	if (npix > MAXOTF)	# sanity check on psf size
	    call printerr_int (toobig, npix)

	npix = npix * nsub
	call salloc (grid, npix+1, TY_REAL)
	call malloc (buffer, numpsf*npix*npix, TY_REAL)

	# Interpolate psf on common pixel grid

	array = buffer
	do ipsf = 0, numpsf-1 {
	    scale = apscale / (nsub * Memr[psfscale+ipsf])
	    call setgrid (scale, Memr[grid], npix+1)

	    call pixinterp (Memi[im+ipsf], Memr[xpcen+ipsf], Memr[ypcen+ipsf],
			    Memr[grid], npix, Memr[array])

	    call imunmap (Memi[im+ipsf])
	    array = array + npix * npix
	}

	# Allocate and fill psf structure

	call malloc (psf, LEN_OTFSTRUCT, TY_STRUCT)

	OTF_NUMBER(psf) = numpsf
	OTF_NXPIX(psf) = npix
	OTF_NYPIX(psf) = npix
	OTF_WAVPTR(psf) = wavepsf
	OTF_BUFFER(psf) = buffer

	# Free temporary arrays

	call mfree (im, TY_INT)
	call sfree (sp)

	return (psf)
end
