include "source.h"
include "grating.h"

#* HISTORY *
#* B.Simon	18-Jul-95	original

# ADDSOURCE -- Compute dispersion solution and add source to output image

procedure addsource (src, grating, ndx, ndy, nix, niy, asrc, out)

pointer	src		# i: source descriptor
real	grating[LEN_G]	# i: grating dispersion parameters
int	ndx		# i: x dimension of masked source
int	ndy		# i: y dimension of masked source
int	nix		# i: x dimension of output image
int	niy		# i: y dimension of output image
real	asrc[ndx,ndy]	# i: masked source
real	out[nix,niy]	# u: output image
#--
int	nwave, m1, m2, order, idx
pointer	sp, ypos, wave, spec, blaze, line, unity
real	ox, oy, flux, frac, xl, yl

int	specconv()

begin
	# Get position and flux of source

	ox = SRC_XPOS(src)
	oy = SRC_YPOS(src)
	flux = SRC_FLUX(src)
	if (flux <= 0.0)
	    return

	# Allocate memory for temporary arrays

	nwave = niy

	call smark (sp)
	call salloc (ypos, nwave, TY_INT)
	call salloc (wave, nwave, TY_REAL)
	call salloc (spec, nwave, TY_REAL)
	call salloc (blaze, nwave, TY_REAL)
	call salloc (line, nwave, TY_REAL)
	call salloc (unity, ndx*ndy, TY_REAL)

	# Normalize source flux to unity

	frac = 1.0 / flux
	call amulkr (asrc, frac, Memr[unity], ndx*ndy)

	m1 = grating[M1]
	m2 = grating[M2]
	do order = m1, m2 {
	    # Convert spectrum from wavelength to pixel space

	    call disperse (ox, oy, order, grating, SRC_NWAVE(src), 
			   SRC_WAVE(src), SRC_SPEC(src), nwave, 
			   Memi[ypos], Memr[wave], Memr[spec])

	    # Print diagnostic message

	    call order_message (order, nwave, Memr[wave])

	    # Compute efficiency of grating

	    call blazefunc (ox, oy, order, grating, 
			    nwave, Memr[wave], Memr[blaze])

	    call amulr (Memr[spec], memr[blaze], Memr[spec], nwave)

	    # Copy spectrum to output array

	    xl = ox - 0.5 * (ndx - 1)
	    yl = oy

	    do idx = 1, ndx {
		# Convolve spectrum with source line
		# Skip if source line has no flux

		if (specconv (nwave, Memr[wave], Memr[spec], order, grating, 
			      idx, ndx, ndy, Memr[unity], Memr[line]) == YES) {

		    # Add convolved spectrum to output

		    call addline (xl, yl, grating, nwave, Memi[ypos],
				  Memr[wave], Memr[line], nix, niy, out)

		}

		xl = xl + 1.0
	    }
	}

	call sfree (sp)
end
