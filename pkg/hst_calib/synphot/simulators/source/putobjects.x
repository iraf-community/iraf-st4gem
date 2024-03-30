include "otf.h"
include "object.h"

#* HISTORY *
#* B.Simon	17-Feb-95	original
#* B.Simon	04-May-95	Changed function descriptor
#* B.Simon	07-Jun-95	Increased size of output object
#* B.Simon	29-Jun-95	Revised so npx != npy
#* B.Simon	11-Nov-95	Added call to psfshift

# PUTOBJECTS -- Put objects in output image, convolved with psf

procedure putobjects (out, nix, niy, obj, psf, apscale, dynrange, nsub)

real	out[nix,niy]	# u: output image buffer
int	nix		# i: first dimension of buffer
int	niy		# i: second dimension of buffer 
pointer	obj		# i: object list descriptor
pointer	psf		# i: point spread function
double	apscale		# i: aperture scale
real	dynrange	# i: dynamic range of object
int	nsub		# i: number of subpixels
#--
int	nweight, idx, npx, npy, ntx, nty, nox, noy
pointer	shp, apsf, spsf, anobj, newfunc, delfunc, valfunc, sizfunc
real	ox, oy, flux

begin
	# Allocate psf array

	npx = OTF_NXPIX(psf)
	npy = OTF_NYPIX(psf)
	nweight = OTF_NUMBER(psf)
	call malloc (apsf, npx*npy, TY_REAL)
	call malloc (spsf, npx*npy, TY_REAL)

	do idx = 1, OBJ_NUMBER(obj) {
	    ox = OBJ_XPOS(obj,idx)
	    oy = OBJ_YPOS(obj,idx)
	    flux = OBJ_FLUX(obj,idx)
	    shp = OBJ_SHAPE(obj,idx)

	    # Print diagnostic message

	    call obj_message (ox,oy, flux, shp)

	    # Calculate psf for this object

	    call otfcalc (psf, ox, oy, nsub, nweight, OBJ_PWEIGHT(obj,idx,1), 
			  npx, npy, Memr[apsf])

	    # Shift the psf to account for the offset of the object position
	    # from the center of the pixel

	    call psfshift (ox, oy, nsub, npx, npy, Memr[apsf], Memr[spsf])

	    if (shp == NULL) {
		# Convolve star  with psf

		ntx = npx / nsub + 2
		nty = npy / nsub + 2

		call malloc (anobj, ntx*nty, TY_REAL)

		call starconv (ox, oy, flux, nsub, npx, npx, ntx, nty, 
			       Memr[spsf], Memr[anobj])
	    } else {
		# Convolve extended object with psf

		call shapefuncs (shp, newfunc, delfunc, valfunc, sizfunc)

		call zcall4 (newfunc, shp, apscale, dynrange, nsub)
		call zcall2 (sizfunc, nox, noy)

		ntx = (nox + npx) / nsub + 2
		nty = (noy + npy) / nsub + 2

		call malloc (anobj, ntx*nty, TY_REAL)

		call objconv (ox, oy, flux, nsub, nox, noy, npx, npy, ntx, nty,
			      valfunc, Memr[spsf], Memr[anobj])

		call zcall0 (delfunc)
	    }

	    # Add convolved object to output buffer

	    call addobject (ox, oy, ntx, nty, nix, niy, Memr[anobj], out)
	    call mfree (anobj, TY_REAL)
	}

	call mfree (apsf, TY_REAL)
	call mfree (spsf, TY_REAL)
end
