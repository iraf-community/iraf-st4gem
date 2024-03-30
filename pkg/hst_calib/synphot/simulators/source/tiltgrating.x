include	<math.h>
include "grating.h"

define	MAXITER		10
define	DELTA_THETA	1.0
define	EPS_WAVE	0.001
define	EPS_THETA	1.0e-5

#* HISTORY *
#* B.Simon	06-Dec-95	original

# TILTGRATING -- Tilt a grating to get a specified central wavelength
#
# This procedure iteratively solves for the values of theta_y and theta_x 
# that would put the specified order and wavelength at the center of the
# detector.

procedure tiltgrating (cenorder, cenwave, grating)

int	cenorder	# i: order to place at center of detector
real	cenwave		# i: wavelength to place at center of detector
real	grating[LEN_G]	# u: grating parameters
#--
int	order, idir, iter
real	theta1, theta2, wave1, wave2, temp, dtheta

string	no_wave  "Central wavelength is undefined"
string	no_order "Central order is undefined"
string	badorder "Illegal central order for this grating"
string	badwave  "Cannot set central wavelength to this value"

real	wave_error()

begin
	# Check central order and wavelength

	if (IS_INDEFI(cenorder) && IS_INDEFR(cenwave))
	    return

	if (IS_INDEFR(cenwave))
	    call printerr_real (no_wave, cenwave)

	order = cenorder
	if (IS_INDEFI(order)) {
	    if (grating[M1] == grating[M2]) {
		order = grating[M1]
	    } else {
		call printerr_int (no_order, order)
	    }
	}

	if (order < grating[M1] || order > grating[M2])
	    call printerr_int (badorder, order)

	# Use dispersion equation to solve for theta given the values of
	# the order and wavelngth. New value of theta overwrites input value
	# Solve first in the y direction, then in the x direction.

	do idir = 1, 2 {
	    # Skip if no grating in this direction

	    if (GCON(grating,idir) <= 0.0)
		next

	    # Solution by secant method, which requires two initial guesses

	    theta1 = THETA(grating,idir)
	    wave1 = wave_error (order, cenwave, idir, grating)

	    theta2 = theta1 + DELTA_THETA
	    THETA(grating,idir) = theta2
	    wave2 = wave_error (order, cenwave, idir, grating)

	    # Swap values if first guess is closer than second

	    if (abs(wave1) < abs(wave2)) {
		temp = theta1
		theta1 = theta2
		theta2 = temp

		temp = wave1
		wave1 = wave2
		wave2 = temp

		THETA(grating,idir) = theta2
	    }

	    for (iter = 1; iter <= MAXITER; iter = iter + 1) {
		# Estimate theta where wavelength error is zero
		# by assuming function can be approximated by a
		# straight line between the two prevous estimates

		dtheta = (theta1 - theta2) * wave2 / (wave2 - wave1)
		if (abs (dtheta) > 90.0)
		    call printerr_real (badwave, cenwave)

		# Discard older of two estimates and compute new estimate

		theta1 = theta2
		wave1 = wave2
		theta2 = theta2 + dtheta

		# Compute wavelength error of new estimate

		THETA(grating,idir) = theta2
		wave2 = wave_error (order, cenwave, idir, grating)

		# Exit when wavelength error is less than tolerance
		# or change in theta is less than tolerance

		if (abs(wave2) < EPS_WAVE || abs (dtheta) < EPS_THETA)
		    break
	    }

	    if (iter > MAXITER)
		call printerr_real (badwave, cenwave)

	    # Set order to one for x dispersion direction
	    order = 1
	}
end

# WAVE_ERROR -- Compute the error btw. the desired and actual wavelength

real procedure wave_error (cenorder, cenwave, idir, grating)

int	cenorder	# i: order to place at center of detector
real	cenwave		# i: wavelength to place at center of detector
int	idir		# i: x or y direction (y = 1, x = 2)
real	grating[LEN_G]	# i: grating parameters
#--
real	inangle, outangle, factor, wave

begin

	inangle = (BETA(grating,idir) + DELTA(grating,idir) +
		   THETA(grating,idir)) / RADIAN

	outangle = (BETA(grating,idir) - DELTA(grating,idir) +
		   THETA(grating,idir)) / RADIAN

	factor = cos (DEGTORAD(SIGMA(grating,idir))) * 
		 GCON(grating,idir) / cenorder

	wave = factor * (sin(inangle) + sin (outangle))

	return (wave - cenwave)
end
