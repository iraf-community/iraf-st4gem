# SYNSPLINE -- Compute coefficients for a natural cubic spline

procedure synspline (numold, xold, yold, y2)

int	numold		# i: number of tabulated points
real	xold[ARB]	# i: independent variable in tabulated function
real	yold[ARB]	# i: dependent variable in tabulated function
real	y2[ARB]		# o: spline coefficients (second derivatives)
#--
int	iold
pointer	sp, d, u, w

begin
	# Allocate dynamic memory for temporary arrays

	call smark (sp)
	call salloc (d, numold, TY_REAL)
	call salloc (u, numold, TY_REAL)
	call salloc (w, numold, TY_REAL)

	# Cubic spline algorithm taken from "Algorithms" by
	# Robert Sedgewick, p. 71

	# Compute tridiagonal matrix coefficients

	do iold = 2, numold - 1
	    Memr[d+iold-1] = 2.0 * (xold[iold+1] - xold[iold-1])

	do iold = 1, numold - 1
	    Memr[u+iold-1] = xold[iold+1] - xold[iold]

	do iold = 2, numold - 1 {
	    Memr[w+iold-1] = (yold[iold+1] - yold[iold]) / Memr[u+iold-1] -
			     (yold[iold] - yold[iold-1]) / Memr[u+iold-2]
	}

	# Forward pass of gaussian elimination

	y2[1] = 0.0
	y2[numold] = 0.0

	do iold = 2, numold - 2 {
	    Memr[w+iold] = Memr[w+iold] - Memr[w+iold-1] * 
			   Memr[u+iold-1] / Memr[d+iold-1]
	    Memr[d+iold] = Memr[d+iold] - Memr[u+iold-1] *
			   Memr[u+iold-1] / Memr[d+iold-1]
	}

	# Back pass

	do iold = numold - 1, 2, -1 {
	    y2[iold] = (Memr[w+iold-1] - Memr[u+iold-1] * 
			y2[iold+1]) / Memr[d+iold-1]
	}

	call sfree (sp)
end
