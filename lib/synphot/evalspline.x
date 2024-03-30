# EVALSPLINE -- Evaluate a cubic spline on a grid

procedure evalspline (numold, xold, yold, y2, numnew, xnew, ynew)

int	numold		# i: number of values in tabulated function
real	xold[ARB]	# i: independent variable in tabulated function
real	yold[ARB]	# i: dependent variable in tabulated function
real	y2[ARB]		# i: spline coefficients of tabulated function
int	numnew		# i: number of values in new grid
real	xnew[ARB]	# i: independent variable in new grid
real	ynew[ARB]	# o: dependent variable in new grid
#--
int	iold, inew
real	dx, a, b

string	badseq  "Spline interpolation error. Data out of sequence."

begin
	if (numold == 1) {
	    call amovkr (yold[1], ynew, numnew)
	    return
	}

	iold = 2
	do inew = 1, numnew {
	    # Find bracket for interpolated value

	    while (xold[iold] < xnew[inew]) {
		if (iold == numold)
		    break

		iold = iold + 1
	    }

	    # Formula from Numerical Recipies, Section 3.3
	    # Modified by factor of six for Sedgewick algorithm

	    dx = xold[iold] - xold[iold-1]
	    if (dx <= 0.0)
		call printerr_int (badseq, iold)

	    a = (xold[iold] - xnew[inew]) / dx
	    b = (xnew[inew] - xold[iold-1]) / dx

	    ynew[inew] = a * yold[iold-1] + b * yold[iold] +  (dx * dx) *
			 ((a ** 3 - a) * y2[iold-1] + (b ** 3 - b) * y2[iold])
	}

end
