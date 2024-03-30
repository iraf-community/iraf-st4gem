include	<math.h>
include "grating.h"

#* HISTORY *
#* B.Simon	01-Dec-95	original

# MAGDISP -- Compute the magnification factor in the dispersion direction
#
# Use the dispersion equation to solve for how a change in position of 
# an object leads to a change in the position of the spectrum at a given
# wavelength. The ratio of the size of the input profile to the output 
# profile detrmines the magnification factor.

real procedure magdisp (nprof, npix, cenwave, order, grating)

int	nprof		# i: lnegth of object profile in pixels
int	npix		# i: number of pixels on a detector side
real	cenwave		# i: central wavelength of spectrum
int	order		# i: spectral order
real	grating[ARB]	# i: grating coefficients
#--
int	ipos
real	pixscale, factor, inoff, outoff, inangle, outangle, mag
real	inextent[2], outextent[2]

begin
	# Compute the angular extent of the input line profile

	pixscale = atan2 (grating[SIZE] / npix, grating[F])

	inextent[1] = 0.5 * pixscale * (nprof - 1)
	inextent[2] = - inextent[1]

	# Compute the corresponding change in the position
	# of a fixed wavelength in the dispersed spectrum

	factor = cos(DEGTORAD(grating[SIGMA_Y])) * grating[GY] / order

	inoff =  (grating[BETA_Y] + grating[DELTA_Y] + 
		  grating[THETA_Y]) / RADIAN

	outoff = (grating[BETA_Y] - grating[DELTA_Y] + 
		  grating[THETA_Y]) / RADIAN

	do ipos = 1, 2 {
	    inangle =  inoff + (inextent[ipos] / RADIAN)
	    outangle = asin (cenwave / factor - sin (inangle))
	    outextent[ipos] = outangle - outoff
	}

	# The ratio of the output to the input extent is the magnification

	mag = abs ((outextent[2] - outextent[1]) / (inextent[2] - inextent[1]))
	return (mag)
end
