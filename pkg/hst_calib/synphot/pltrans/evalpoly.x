# EVALPOLY -- Evaluate a polynomial at an evenly space x grid given 
# xmin, xmax and coeffs.  Returns an array of x and y values.

procedure evalpoly(xmin, xmax, coeff, nterm, npt, xfit, yfit )

real	xmin		# i: minimum x value
real	xmax		# i: maximum x value
real	coeff[ARB]	# i: array of coefficient values
int	nterm		# i: number of terms ( ndegree + 1)
int	npt		# i: number of points to evaluate polynomial at
real	xfit[ARB]	# o: array of xvalues
real	yfit[ARB]	# o: array of yvalues

int	ic, jc
real	dx

begin

	# Determine the x-axis spacing

	dx = (xmax - xmin)/(npt - 1.)

	# Evaluate the x-coordinate and the polynomial at each x
	do ic = 1, npt {
	   xfit[ic] = xmin + (ic - 1) * dx

	   yfit[ic] = coeff[nterm]
	   do jc = nterm - 1, 1, -1
	      yfit[ic] = yfit[ic] * xfit[ic] + coeff[jc]

	}
	      
end
