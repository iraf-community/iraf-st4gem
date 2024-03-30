include	<pkg/dttext.h>
include	<math/curfit.h>

# EVCURFIT -- Evaluate a function that has been by the curfit routine, given
# the coefficients, x points, and order of fit.  Will evaluate Chebyshev,
# Legendre, spline and cubic spline.

procedure evcurfit( cv, xmin, xmax, xfit, yfit, nfit )

pointer	cv		# i: pointer to curfit data structure
real	xmin		# i: min x value to evaluate
real	xmax		# i: max x value to evaluate
real	xfit[ARB]	# o: points were function should be evaluated
real	yfit[ARB]	# o: y values of curve
int	nfit		# i: number of points to fit

real	dx
int	ic

begin

	dx = (xmax - xmin)/ max( 1, nfit-1)
	do ic = 1, nfit
	   xfit[ic] = xmin + (ic-1) * dx

	call cvvector( cv, xfit, yfit, nfit )

end
