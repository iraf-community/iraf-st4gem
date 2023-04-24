include "nlfit.h"

# NL_FDEV --  Compute function and its derivatives.
#
# This procedure is called by the Levenberg-Marquardt minimization method.

procedure nl_fdev (nl, ex, ey, coeff, ez, deriv)

pointer	nl		# i: Curve descriptor.
real	ex		# i: Independent x variable.
real	ey		# i: Independent y variable.
real	coeff[ARB]	# i: Coefficients.
real	ez		# o: Evaluated function.
real	deriv[ARB]	# o: Derivatives.
#--
int	i, npar
real	x[1],y[1],z[1]
real	ac

int	nl_stati()
real	nl_defc()

errchk	nl_eval

begin
	# Compute function at given point.
	x[1] = ex
	y[1] = ey
	call nl_eval (nl, x, y, z, coeff, 1)
	ez = z[1]

	# Compute derivatives by stepping 10% away from locus in 
        # coefficient space. If coefficient is exactly zero,
        # use function-dependent value.
	npar = nl_stati (nl, "npar")
	do i = 1, npar {
	    ac = coeff[i]
	    if (coeff[i] == 0.0)
	        coeff[i] = nl_defc (nl, i, ex, ez)
	    if (coeff[i] != 0.) {
	        coeff[i] = 1.1 * coeff[i]
	        call nl_eval (nl, x, y, z, coeff, 1)
	        deriv[i] = (z[1] - ez) / (coeff[i] - ac)
	    } else
	        deriv[i] = 0. 
	    coeff[i] = ac
	}
end
        
                                                               
                     
