include "nlfit.h"

# NL_XSCALE --  Compute scaling factors for independent variables. 
# The same scaling is applied to both X and Y variables, to simplify 
# computations with 2-d Gaussians.

procedure nl_xscale (nl, x, y, npts)

pointer nl			# i: Curve descriptor.
real	x[ARB]			# i: Independant x variable array.
real	y[ARB]			# i: Independant y variable array.
int	npts			# i: Size of arrays

int	i, nl_units()
real	aux

errchk	nl_units, nl_scale

begin
	# Check independent variable units.
	if (nl_units (nl, x, npts) == ERR)
	    call error (0, "Impossible to fit.")

	# Do not scale if user function or Levenberg-Marquardt method.
	if (NL_FITFUNC(nl) == USER ||
            NL_METHOD(nl)  == MARQUARDT) {

	    NL_XSCALE(nl) = 1.
	    NL_YSCALE(nl) = 1.

	} else {

	    # Compute X scaling factor.
	    NL_XSCALE(nl) = 0.5 * (abs(x[1]) + abs(x[npts]))
	    if (NL_XSCALE(nl) == 0)
	        i = 0
	    else
	        i = -log10(NL_XSCALE(nl))
	    NL_XSCALE(nl) = 10.0d0 ** i

	    # Compute Y scaling factor.
	    if (NL_DIM(nl) == 2) {
	        NL_YSCALE(nl) = 0.5 * (abs(y[1]) + abs(y[npts]))
	        if (NL_YSCALE(nl) == 0)
	            i = 0
	        else
	            i = -log10(NL_YSCALE(nl))
	        NL_YSCALE(nl) = 10.0d0 ** i
	    } else 
	        NL_YSCALE(nl) = 0.

	    # Set same scaling factor for both independent variables.
	    if (NL_DIM(nl) == 2) {
	        aux = (NL_XSCALE(nl) + NL_YSCALE(nl)) / 2.
	        NL_XSCALE(nl) = aux
	        NL_YSCALE(nl) = aux
	    }   
	}
end
