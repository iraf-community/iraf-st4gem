include "source.h"

#* HISTORY *
#* B.Simon	12-Jul-95	original
#* B.Simon	11-Nov-95	Added call to psfshift

# PUTSOURCE -- Add source spectrum to output image

procedure putsource (src, apshape, apscale, grating, 
		     dynrange, nsub, nix, niy, out)

pointer	src		# i: source descriptor
pointer	apshape		# i: aperture shape descriptor
double	apscale		# i: detector pixel size
real	grating[ARB]	# i: grating dispersion parameters
real	dynrange	# i: dynamic range of object
int	nsub		# i: number of subpixels
int	nix		# i: first dimension of output buffer
int	niy		# i: second dimension of output buffer 
real	out[nix,niy]	# o: output buffer
#--
pointer	sp, mask, anobj, asrc, apsf, spsf
pointer	newfunc, delfunc, valfunc, sizfunc
int	npx, npy, ndx, ndy, ntx, nty, nox, noy
real	ox, oy, flux

int	makemask()

begin
	# Get position and flux of source

	ox = SRC_XPOS(src)
	oy = SRC_YPOS(src)
	flux = SRC_FLUX(src)

	# Write diagnostic message

	call obj_message (ox, oy, flux, SRC_SHAPE(src))
	# Get source psf

	npx = SRC_NXPSF(src)
	npy = SRC_NYPSF(src)
	apsf = SRC_PSFPTR(src)

	# Calculate size of convolved source

	if (SRC_SHAPE(src) == NULL) {
	    ntx = npx / nsub + 2
	    nty = npy / nsub + 2

	} else {
	    call shapefuncs (SRC_SHAPE(src), newfunc, 
			     delfunc, valfunc, sizfunc)

	    call zcall4 (newfunc, SRC_SHAPE(src), apscale, dynrange, nsub)
	    call zcall2 (sizfunc, nox, noy)

	    ntx = (nox + npx) / nsub + 2
	    nty = (noy + npy) / nsub + 2
	}

	ndx = ntx
	ndy = nty + SRC_NXLSF(src) / nsub

	# Allocate memory for convoved source and aperture mask

	call smark (sp)
	call salloc (mask, ntx*nty, TY_REAL)
	call salloc (anobj, ntx*nty, TY_REAL)
	call salloc (asrc, ndx*ndy, TY_REAL)
	call salloc (spsf, npx*npy, TY_REAL)

	# Shift the psf to account for the offset of the object position
	# from the center of the pixel

	call psfshift (ox, oy, nsub, npx, npy, Memr[apsf], Memr[spsf])

	# Create aperture mask. If mask does not overlap source, skip source

	if (makemask (apshape, apscale, ox, oy, nsub, 
		      nix, niy, ntx, nty, Memr[mask]) == YES) {

	    # Convolve source with point spread function

	    if (SRC_SHAPE(src) == NULL) {
		# Point source convolution routine

		call starconv (ox, oy, flux, nsub, npx, 
			       npy, ntx, nty, Memr[spsf], 
			       Memr[anobj])

	    } else {
		# Extended source convolution routine

		call objconv (ox, oy, flux, nsub, nox, noy, npx, 
			      npy, ntx, nty, valfunc, Memr[spsf], 
			      Memr[anobj])

		call zcall0 (delfunc)
	    }

	    # Mask the convolved source

    	    call addmask (ntx, nty, Memr[mask], Memr[anobj], flux)
	    SRC_FLUX(src) = flux

	    # Convolve masked source with line spread function

	    if (SRC_LSFPTR(src) == NULL) {
		call amovr (Memr[anobj], Memr[asrc], ntx*nty)

	    } else {
		call lsfconv (flux, nsub, SRC_NXLSF(src), SRC_LSF(src), 
			      ntx, nty, Memr[anobj], ndx, ndy, Memr[asrc])
	    }

	    # Compute dispersion solution and add source to output image

	    call addsource (src, grating, ndx, ndy, 
			    nix, niy, Memr[asrc], out)
	}

	call sfree (sp)

end
