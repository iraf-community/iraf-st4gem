include	<math.h>
include "grating.h"

#* HISTORY *
#* B.Simon	17-Jul-95	original

# CROSSDISP -- Compute cross dispersion solution from wavelength set

procedure crossdisp (ox, oy, grating, nwave, wave, xpos)

real	ox		# i: source x position
real	oy		# i: source y position
real	grating[ARB]	# i: grating dispersion coefficients
int	nwave		# i: length of wavelength set
real	wave[ARB]	# i: wavelength set
real	xpos[ARB]	# o: x position of dispersed wavelength
#--
pointer	iwave
real	pixscale, offset, outoff, factor, inangle, outangle

begin
	# Check for no cross-disperser

	if (grating[GX] <= 0.0) {
	    call amovkr (ox, xpos, nwave)
	    return
	}

	# Compute incident (input) and reflection (output) angles

	pixscale = atan2 (grating[SIZE] / nwave, grating[F])
	offset = pixscale * (0.5 * (nwave + 1) - ox)
	inangle = (offset + grating[BETA_X] + grating[DELTA_X] + 
		   grating[THETA_X] ) / RADIAN

	outoff = (offset + grating[BETA_X] - grating[DELTA_X] + 
		  grating[THETA_X] ) / RADIAN

	factor = cos(DEGTORAD(grating[SIGMA_X])) * grating[GX]

	# Compute positions corresponding to wavelengths

	do iwave = 1, nwave {
	    outangle = asin (wave[iwave] / factor - sin(inangle))
	    xpos[iwave] = ox + (outangle - outoff) / pixscale
	}

end
