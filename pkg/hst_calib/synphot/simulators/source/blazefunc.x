include	<math.h>
include "grating.h"

#* HISTORY *
#* B.Simon	14-Jul-95	original

# BLAZEFUNC -- Compute grating blaze correction to spectrum

procedure blazefunc (ox, oy, order, grating, nwave, wave, blaze)

real	ox		# i: source x position
real	oy		# i: source y position
int	order		# i: spectral order
real	grating[ARB]	# i: grating dispersion coefficients
int	nwave		# i: length of wavelength array
real	wave[ARB]	# i: wavelength array
real	blaze[ARB]	# o: blaze function
#--
int	iwave
pointer	sp, xpos, ypos
real	pixscale, offset, outoff, inangle, outangle, small, large, temp

begin
	# Allocate memory for temporary array

	call smark (sp)
	call salloc (xpos, nwave, TY_REAL)
	call salloc (ypos, nwave, TY_REAL)

	# Compute grid of pixel positions

	pixscale = atan2 (grating[SIZE] / nwave, grating[F])
	call setgrid (pixscale, Memr[ypos], nwave)
	offset = pixscale * (0.5 * (nwave + 1) - oy)

	# Compute incident (input) and reflection (output) angles

	inangle = (offset + grating[BETA_Y] + grating[DELTA_Y] + 
		   grating[THETA_Y] ) / RADIAN

	outoff = (offset + grating[BETA_Y] - grating[DELTA_Y] + 
		  grating[THETA_Y] ) / RADIAN

	# Compute echelle efficiency

	do iwave = 1, nwave {
	    outangle = outoff + Memr[ypos+iwave-1] 
	    small = min (inangle, outangle)	# bb
	    large = max (inangle, outangle)	# aa

	    temp = PI * order * cos (large) * sin (0.5 * (small + large) -
		   DEGTORAD(grating[BETA_Y])) / sin (0.5 * (small + large))

	    blaze[iwave] = cos (large) / cos (small)
	    if (abs(temp) > 1.0e-4)
		blaze[iwave] = blaze[iwave] * (sin (temp) / temp) ** 2
	}

	# Correct for cross-disperser

	if (grating[GX] > 0.0) {
	    # Compute cross disperser offsets from wavelengths
	    # offset assumes square detector

	    call crossdisp (ox, oy, grating, nwave, wave, Memr[xpos])
	    offset = pixscale * (0.5 * (nwave + 1) - ox)	

	    inangle = (offset + grating[BETA_X] + grating[DELTA_X] + 
		       grating[THETA_X] ) / RADIAN

	    outoff = (offset + grating[BETA_X] - grating[DELTA_X] + 
		      grating[THETA_X] ) / RADIAN
			   
	    # Compute incident (input) and reflection (output) angles

	    do iwave = 1, nwave {
		outangle = outoff + pixscale * Memr[xpos+iwave-1] / RADIAN

		small = min (inangle, outangle)	# bb
		large = max (inangle, outangle)	# aa

		temp = PI * cos (large) * sin (0.5 * (small + large) -
		       DEGTORAD(grating[BETA_X])) / sin (0.5 * (small + large))

		blaze[iwave] = blaze[iwave] * cos (large) / cos (small)
		if (abs(temp) > 1.0e-4)
		    blaze[iwave] = blaze[iwave] * (sin (temp) / temp) ** 2
	    }
	}

	call sfree (sp)
end
