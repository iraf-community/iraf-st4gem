include	<imhdr.h>
include "otf.h"

#* HISTORY *
#* B.Simon	10-Jul-95	derived from psfimage

# LSFIMAGE -- Create a structure holding the line spread functions

pointer procedure lsfimage (lsfcat, obsmode, apscale, dynrange, nsub)

char	lsfcat[ARB]	# i: catalog of lsf file names
char	obsmode[ARB]	# i: observation mode
double	apscale		# i: aperture scale (degrees)
real	dynrange	# i: dynamic range of lsf
int	nsub		# i: pixel subsampling
#--
int	numlsf, idx, npix
pointer	sp, pcenter, lsfscale, image
pointer	lsf, im, wavelsf, grid, buffer, array
real	range, scale

string	badsize  "LSF has invalid scale; detector scale assumed"
string	toobig   "LSF size too large"

begin
	# Retrieve lsf images that match this obsmode

	call otfopen (lsfcat, obsmode, im, wavelsf, numlsf)

	# Allocate arrays based on number of lsfs

	call smark (sp)
	call salloc (pcenter, numlsf, TY_REAL)
	call salloc (lsfscale, numlsf, TY_REAL)
	call salloc (image, SZ_FNAME, TY_CHAR)

	# Compute largest extent of any lsf

	range = 0.0
	do idx = 0, numlsf-1 {
	    call lsfregion (Memi[im+idx], dynrange, npix, Memr[pcenter+idx], 
			    Memi[lsfscale+idx])

	    if (Memr[lsfscale+idx] <= 0.0) {
		call strcpy (IM_HDRFILE(Memi[im+idx]), Memc[image], SZ_FNAME)
		call synphotwarn (badsize, Memc[image])
		Memr[lsfscale+idx] = apscale
	    }

	    range = max (npix * Memr[lsfscale+idx], range)
	}

	# Compute lsf array size and allocate arrays based on this size

	npix = int (range / apscale) + 1
	if (npix > MAXOTF)	# sanity check on lsf size
	    call printerr_int (toobig, npix)

	npix = npix * nsub
	call salloc (grid, npix+1, TY_REAL)
	call malloc (buffer, numlsf*npix, TY_REAL)

	# Interpolate lsf on common pixel grid

	array = buffer
	do idx = 0, numlsf-1 {
	    scale = apscale / (nsub * Memr[lsfscale+idx])
	    call setgrid (scale, Memr[grid], npix+1)

	    call lsfinterp (Memi[im+idx], Memr[pcenter+idx],
			    Memr[grid], npix, Memr[array])

	    call imunmap (Memi[im+idx])
	    array = array + npix
	}

	# Allocate and fill lsf structure

	call malloc (lsf, LEN_OTFSTRUCT, TY_STRUCT)

	OTF_NUMBER(lsf) = numlsf
	OTF_NXPIX(lsf) = npix
	OTF_NYPIX(lsf) = 1
	OTF_WAVPTR(lsf) = wavelsf
	OTF_BUFFER(lsf) = buffer

	# Free temporary arrays

	call mfree (im, TY_INT)
	call sfree (sp)

	return (lsf)
end
