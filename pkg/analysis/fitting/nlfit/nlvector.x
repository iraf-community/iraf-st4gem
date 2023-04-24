include "nlfit.h"

# NL_VECTOR -- Procedure to compute the vector of function values, given
# the independent (unscaled) variable. It assumes that unscaled coefficients
# were already stored in NL_PARAMS(nl).

procedure nl_vector (nl, x, y, zfit, npts)

pointer	nl			# i: Curve descriptor.
real	x[ARB]			# i: Independent x variable array.
real	y[ARB]			# i: Independent y variable array.
real	zfit[ARB]		# o: Evaluated dependent variable array.
int	npts			# i: Number of data points

include "nluserf.com"

errchk	nl_eval, nl_xscale, nl_scale

begin
	# Scale independent variables.
	call nl_xscale (nl, x, y, npts)

	# Z's not scaled
	NL_ZSCALE(nl) = 1.

	call realloc (NL_SX(nl), npts, TY_REAL)
	call realloc (NL_SY(nl), npts, TY_REAL)
	call amulkr (x, NL_XSCALE(nl), Memr[NL_SX(nl)], npts)
	call amulkr (y, NL_YSCALE(nl), Memr[NL_SY(nl)], npts)

	# If user function, realloc common.
	if (NL_FITFUNC(nl) == USER)
	    call realloc (lx, npts, TY_REAL)

	# Scale input parameters.
	call nl_scale (nl)

	# Evaluate function
	call nl_eval (nl, Memr[NL_SX(nl)], Memr[NL_SY(nl)], zfit,
	              Memr[NL_SPARAMS(nl)], npts)
end
