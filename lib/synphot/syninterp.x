include "libsynphot.h"

# SYNINTERP -- Interpolate a tabulated function on a new grid

procedure syninterp (numold, xold, yold, numnew, xnew, ynew)

int	numold		# i: number of values in tabulated function
real	xold[ARB]	# i: independent variable in tabulated function
real	yold[ARB]	# i: dependent variable in tabulated function
int	numnew		# i: number of values in new grid
real	xnew[ARB]	# i: independent variable in new grid
real	ynew[ARB]	# o: dependent variable in new grid
#--

begin
	# Interpolate, setting extrapolated values to zero

	call synextrap (0.0, numold, xold, yold, numnew, xnew, ynew)
end

# SYNEXTRAP -- Interpolate, setting extapolated values to a constant

procedure synextrap (extrap, numold, xold, yold, numnew, xnew, ynew)

real	extrap		# i: extrapolated value
int	numold		# i: number of values in tabulated function
real	xold[ARB]	# i: independent variable in tabulated function
real	yold[ARB]	# i: dependent variable in tabulated function
int	numnew		# i: number of values in new grid
real	xnew[ARB]	# i: independent variable in new grid
real	ynew[ARB]	# o: dependent variable in new grid
#--
int	iold, inew
real	dx, a, b

begin
	# Find a bracket such that xold[iold-1] < xnew[inew] < xold[iold]
	# Then interpolate to get ynew[inew]

	iold = 1
	inew = 1
	while (iold <= numold && inew <= numnew) {
	    # Change index for new bracket

	    if (xnew[inew] > xold[iold]) {
		iold = iold + 1

	    # Interpolate for ynew

	    } else if (xnew[inew] < xold[iold]) {
		# No lower end for bracket, so set to extrapolated value

		if (iold == 1) {
		    ynew[inew] = extrap

		# Both ends of bracket found, so do interpolation

		} else {
		    dx = xold[iold] - xold[iold-1]
		    a = (xold[iold] - xnew[inew]) / dx
		    b = 1.0 - a

		    ynew[inew] = a * yold[iold-1] + b * yold[iold]
		}
		inew = inew + 1

	    # Equality handled as a special case for speed

	    } else {
		ynew[inew] = yold[iold]
		iold = iold + 1
		inew = inew + 1
	    }
	}

	# Set ynew without upper bracket to extrapolated value

	if (inew <= numnew)
	    call amovkr (extrap, ynew[inew], numnew-inew+1)

end
