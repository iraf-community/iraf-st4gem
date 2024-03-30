include	<math.h>
include	"simtwo.h"
include "grating.h"

#* HISTORY *
#* B.Simon	14-Jul-95	original

# DISPERSE -- Compute dispersion solution for spectrum

procedure disperse (ox, oy, order, grating, nowave, owave, ospec, 
		    nwave, ypos, wave, spec)

real	ox		# i: source x position
real	oy		# i: source y position
int	order		# i: spectral order
real	grating[LEN_G]	# i: grating dispersion coefficients
int	nowave		# i: length of input wavelength set
real	owave[ARB]	# i: input wavelength set
real	ospec[ARB]	# i: input spectrum
int	nwave		# i: length of resampled wavelength set
int	ypos[ARB]	# o: y position of resampled spectrum
real	wave[ARB]	# o: resampled wavelength set
real	spec[ARB]	# o: resampled spectrum
#--
int	iwave, done
pointer	sp, pixpos
real	pixscale, offset, outoff, factor, inangle, outangle

string	objform   OBJ_UNITS

int	phottoany()

begin
	# Allocate memory for temporary array

	call smark (sp)
	call salloc (pixpos, nwave, TY_REAL)

	# Compute grid of pixel positions

	pixscale = atan2 (grating[SIZE] / nwave, grating[F])
	call setgrid (pixscale, Memr[pixpos], nwave)
	offset = pixscale * (0.5 * (nwave + 1) - oy)

	# Compute incident (input) and reflection (output) angles

	inangle = (offset + grating[BETA_Y] + grating[DELTA_Y] + 
		   grating[THETA_Y] ) / RADIAN

	outoff = (offset + grating[BETA_Y] - grating[DELTA_Y] + 
		  grating[THETA_Y] ) / RADIAN

	factor = cos(DEGTORAD(grating[SIGMA_Y])) * grating[GY] / order

	# Compute wavelengths which fall at pixel centers

	do iwave = 1, nwave {
	    outangle = outoff + Memr[pixpos+iwave-1] 

	    # Wavelength increases with decreasing y pixel number

	    ypos[iwave] = (nwave - iwave) + 1
	    wave[iwave] = factor * (sin(inangle) + sin(outangle))
	}

	# Interpolate flux on new wavelength set

	call syninterp (nowave, owave, ospec, nwave, wave, spec)

	# Convert spectrum to units of counts

	done = phottoany (objform, nwave, wave, spec)

	call sfree (sp)
end
